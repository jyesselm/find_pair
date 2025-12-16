# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Modern C++ rewrite of X3DNA v2.4 for nucleic acid structure analysis. The goal is 100% output match with the legacy C code.

**Key principle**: Modern code (`src/`, `include/`) must produce JSON output identical to legacy code (`org/`).

## Important Directives

- **Always use max workers/cores** when running parallel operations
- **Always verify file generation** - check that JSON files are valid after generation
- **Use `$X3DNA` environment variable** - legacy code requires `X3DNA=/Users/jyesselman2/local/installs/x3dna`
- **Use `fp2-validate` CLI systematically** for all validation
- **Floating point differences are NOT a valid explanation** for mismatches

## Build Commands

```bash
# Build modern code (Release)
make release

# Build legacy code (Release)
make org-release

# Build debug versions
make debug
make org-debug

# Run C++ tests
make test

# Clean builds
make clean-all
make org-clean-all
```

## Environment Setup

```bash
# Required for LEGACY code execution only (modern code is self-contained)
export X3DNA=/Users/jyesselman2/local/installs/x3dna

# Install Python package
pip install -e ".[dev]"
```

Note: Modern C++ code uses resources in `resources/config/` and `resources/templates/` and does NOT require the X3DNA environment variable.

## Testing

Primary testing is JSON regression comparison - modern output must match legacy output exactly.

```bash
# Main validation CLI (always use fp2-validate)
# Test sets: 10, 50, 100, 500, 1000, or "fast" (3602 validated PDBs)
fp2-validate validate core --test-set fast    # All fast PDBs, skip H-bonds (recommended)
fp2-validate validate core --test-set 100     # 100-PDB test set, skip H-bonds
fp2-validate validate steps --test-set fast   # Step validation only
fp2-validate validate pairs --test-set 100    # Pair validation only
fp2-validate validate 3 --pdb 1EHZ -v         # Single PDB, verbose

# Stage groups: atoms, frames, pairs, steps, hbonds, core (all except hbonds), all
# Use "core" instead of "all" to skip H-bond validation (known detection differences)

# Run pytest
pytest tests_python/

# Regenerate JSON files (always use rebuild_json.py CLI)
python3 scripts/rebuild_json.py regenerate --test-set fast --modern-only   # Fast PDBs, modern only
python3 scripts/rebuild_json.py regenerate --test-set 100                  # 100-PDB test set, both
python3 scripts/rebuild_json.py regenerate 1EHZ 2BNA                       # Specific PDBs
```

### Current Test Results

**Fast PDB Set (3602 PDBs) - Core Stages (excludes H-bonds):**

| Stage | Pass Rate | Notes |
|-------|-----------|-------|
| Core (all except H-bonds) | 99.0% (3566/3602) | Recommended validation |
| Pairs | 99.9% (3598/3602) | 4 failures |
| Frames | 99.7% (3591/3602) | 11 failures |
| Steps | 99.4% (3579/3602) | 23 failures (legacy JSON indexing issues) |

**100-PDB Test Set:**

| Stage | Pass Rate | Notes |
|-------|-----------|-------|
| 1: Atoms | 100% (100/100) | All pass |
| 2: Residue Indices | 100% (100/100) | All pass |
| 3-5: Frames | 100% (100/100) | All pass |
| 6-7: Pair Validation | 100% (100/100) | All pass |
| 8: H-bond List | 46% (46/100) | Known issue: modern uses stricter H-bond criteria |
| 9: Base Pair | 100% (100/100) | All pass |
| 10: Bestpair Selection | 100% (100/100) | All pass |
| 11-12: Steps | 99% (99/100) | 1 failure: 8RUJ (five2three strand swap algorithmic difference) |

**Note**: Step ordering (`bp_idx`) now matches legacy 100% when helices match. Fixes implemented:
- Fixed `helix_break` cutoff: 7.5 → 7.8 (matches legacy `misc_3dna.par` config)
- Fixed `end_stack_xang` threshold: 110° → 125° (matches legacy `END_STACK_XANG`)
- Fixed backbone extraction to use `legacy_residue_idx` from atoms
- Added second `check_direction` call after `check_strand2` (legacy line 1361)
- Fixed neighbor swapping in `calculate_context` (legacy lines 931-941)
- Added Watson-Crick pair check in `wc_bporien` (legacy `base_pairs[m][3] > 0`)
- Fixed modified base handling in wc_bporien ("Cg" → "CG")
- Added geometric bpid check (dir_x > 0 && dir_y < 0 && dir_z < 0)
- Fixed `calculate_context` to not store neighbor1 when dist > helix_break (matches legacy line 963-965)
- Fixed `has_positive_bpid` to check pair type (WC/wobble only) before geometry check
- Fixed pair type classification to use case-insensitive comparison (matches legacy toupper)
- Fixed helix_organization to swap base_i/base_j display order when strand_swapped=True (matches legacy format)
- Fixed `check_others` to compare SPECIFIC frame pairs instead of sum totals (matches legacy lines 1342-1349)
- **Fixed `has_positive_bpid` to use rotation matrix ROWS not COLUMNS for direction calculation** (legacy uses `&orien[i][0]` which is row 0, not column 0)
- **Fixed atoms JSON output to include insertion codes** (PDB column 27, e.g., "520A" has insertion='A')

**Remaining 1 failure** (8RUJ) has algorithmic differences in five2three strand swap logic for one "outlier" pair:
- bp_idx 148 (helix 11, pos 143): legacy=swapped=false, modern=swapped=true
This is a non-sequential pair (347,348 A-A) inserted in helix 11. The helix ordering matches exactly, but strand assignment differs due to complex cross-strand connectivity. Affects steps 142-143, 143-144 (2 mismatches).

**Tolerance settings**: Step comparison uses 5e-4 tolerance in `step_comparison.py` to account for floating point variations. Verbose reporter uses 1e-4 tolerance for display.

### Validation Stages

| Stage | Name | Test Command |
|-------|------|--------------|
| **core** | All except H-bonds | `fp2-validate validate core --test-set fast` (recommended) |
| 1 | Atoms | `fp2-validate validate 1` |
| 2 | Residue Indices | `fp2-validate validate 2` |
| 3-5 | Frames | `fp2-validate validate frames` |
| 6-7 | Pair Validation | `fp2-validate validate pairs` |
| 8 | H-bond List | `fp2-validate validate 8` (known differences) |
| 9 | Base Pair | `fp2-validate validate 9` |
| 10 | Bestpair Selection | `fp2-validate validate 10` |
| 11-12 | Step/Helical Params | `fp2-validate validate steps` |

## Architecture

### Two-Phase Pipeline

1. **find_pair** (`FindPairProtocol`): Base pair identification
   - Frame calculation via least-squares template fitting
   - Geometric validation (dorg, dNN, plane_angle, d_v)
   - H-bond detection and quality scoring
   - Greedy mutual selection of best pairs

2. **analyze** (`AnalyzeProtocol`): Step parameter calculation
   - Calculates 6 step params: Shift, Slide, Rise, Tilt, Roll, Twist
   - Calculates helical parameters

### Key Directories

- `include/x3dna/`, `src/x3dna/` - Modern C++ library
- `org/` - Legacy C code (reference, don't edit except for JSON output)
- `apps/` - Main executables (find_pair_app, analyze_app)
- `tools/` - Utilities (generate_modern_json)
- `x3dna_json_compare/` - Python comparison library
- `data/json/` - Modern JSON output
- `data/json_legacy/` - Legacy JSON output (reference)
- `resources/templates/` - Base frame templates (Atomic_*.pdb)

### Core Algorithm Classes

- `BaseFrameCalculator` - Reference frame calculation via LS fitting
- `BasePairFinder` - Greedy pair finding with mutual selection
- `BasePairValidator` - Geometric checks, H-bonds, quality scores
- `ParameterCalculator` - Step and helical parameter calculation

## Coding Guidelines (from .cursorrules)

- Modern code in `src/`, `include/`; legacy in `org/`
- Keep functions under 50 lines, max 3 levels of nesting
- Use RAII, avoid raw pointers for ownership
- Prefer pass by const reference or rvalue reference
- Modern code stores legacy indices for direct comparison (1-based vs 0-based)
- Use pytest for Python tests in `tests_python/`
- Use multiple cores when possible
- Floating point differences are NOT a valid explanation for mismatches

## Key Files for Debugging

- `scripts/compare_json.py` - Main JSON comparison tool
- `scripts/rebuild_json.py` - Regenerate JSON files
- `x3dna_json_compare/cli.py` - fp2-validate CLI entry point
- `docs/CODE_FLOW.md` - Detailed algorithm flow
- `docs/TESTING_GUIDE.md` - Complete testing reference

## Common Workflows

### Debug a failing PDB
```bash
# Find the failure
fp2-validate validate 3 --test-set 100 -v -s

# Detailed comparison
fp2-validate compare <PDB_ID> --verbose

# Regenerate JSON
./build/generate_modern_json data/pdb/<PDB>.pdb data/json --stage=frames
```

### Add new comparison logic
1. Add comparison module in `x3dna_json_compare/`
2. Add test in `tests_python/integration/`
3. Verify match with `fp2-validate validate <stage>`
