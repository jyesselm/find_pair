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

# Regenerate JSON files (use fp2-validate rebuild)
fp2-validate rebuild regenerate --test-set fast --modern-only   # Fast PDBs, modern only
fp2-validate rebuild regenerate --test-set 100                  # 100-PDB test set, both
fp2-validate rebuild regenerate 1EHZ 2BNA                       # Specific PDBs
```

**Note**: By default, `fp2-validate validate` automatically regenerates modern JSON before comparison. Use `--skip-regen` to skip regeneration and use existing files.

### Current Test Results (December 2024)

**res_id Matching**: Comparison now uses stable `res_id` identifiers (format: `chain-name-seq[ins]`, e.g., "A-G-1", "A-C-10A") for matching records between legacy and modern JSON. This ensures accurate matching regardless of index ordering differences.

**Fast PDB Set (3602 PDBs) - Core Stages (excludes H-bonds):**

| Stage | Pass Rate | Notes |
|-------|-----------|-------|
| Core (all except H-bonds) | 98.9% (3563/3602) | Recommended validation |
| Pairs | 99.9% (3599/3602) | 3 failures |
| Frames | 99.5% (3585/3602) | 17 failures |
| Steps | 99.4% (3582/3602) | 20 failures |

**100-PDB Test Set:**

| Stage | Pass Rate | Notes |
|-------|-----------|-------|
| 1: Atoms | 100% (100/100) | All pass |
| 2: Residue Indices | 100% (100/100) | All pass |
| 3-5: Frames | 99% (99/100) | 1 failure (6OZK) |
| 6-7: Pair Validation | 100% (100/100) | All pass |
| 8: H-bond List | 46% (46/100) | Known issue: modern uses stricter H-bond criteria |
| 9: Base Pair | 100% (100/100) | All pass |
| 10: Bestpair Selection | 100% (100/100) | All pass |
| 11-12: Steps | 100% (100/100) | All pass |

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

**Tolerance settings**: Step comparison uses 5e-4 tolerance in `step_comparison.py` to account for floating point variations. Verbose reporter uses 1e-4 tolerance for display.

**C++ Unit Tests:**

| Category | Pass Rate | Notes |
|----------|-----------|-------|
| Unit Tests | 383/384 (99.7%) | 1 external GEMMI test not built |
| Integration Tests | ~92 skipped | Require real PDB files |

Run with `make test` or `ctest --test-dir build`.

### Large Structure Validation (December 2024)

**Key Finding**: Step parameter calculations are **correct**. Failures in large structures are due to different helix organization, not calculation errors.

Analysis of 1VQ5 (ribosome, 1535 pairs, 1157 steps):
- Base pairs identified: **Identical** (same 1535 residue combinations)
- Helix organization: **Different** (bp_idx ordering differs significantly)
- Steps with same residues + same frame: **100% parameter match**
- The 19 "mismatches" reported are **false position matches** (coincidental mst_org proximity)

**Slow PDB Set** (521 PDBs > 15 seconds, in `data/slow_pdbs.json`):
- ~50% pass rate (vs 99.4% for fast set)
- Failures are helix organization differences, not calculation errors
- All failures are in complex multi-helix structures (ribosomes, large RNA)

**Conclusion**: When legacy and modern use the same residue pairs AND same frame selection (strand_swapped), parameters match exactly. Differences arise from the `five2three` helix organization algorithm producing different orderings for complex structures.

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

### Comparison Targets

The validation system supports multiple comparison targets. Each target corresponds to a directory under `data/`:

| Target | Directory | Description |
|--------|-----------|-------------|
| `legacy` | `data/json_legacy/` | Legacy C code output (X3DNA v2.4 reference) - default |
| `baseline` | `data/json_baseline/` | Snapshot of modern output for regression testing |

```bash
# Compare against legacy (default)
fp2-validate validate core --test-set 100

# Compare against baseline
fp2-validate validate core --test-set 100 --target baseline

# List available targets
fp2-validate targets

# Create/update baseline from current modern output
fp2-validate baseline create          # Create new baseline
fp2-validate baseline update --force  # Overwrite existing baseline
fp2-validate baseline status          # Show baseline status
```

**Baseline workflow**: Use baseline for regression testing. After making changes, compare against baseline to detect unintended changes. Update baseline when intentional changes are made.

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
- `data/json/` - Modern JSON output (generated)
- `data/json_legacy/` - Legacy JSON output (reference)
- `data/json_baseline/` - Baseline snapshot for regression testing
- `resources/templates/` - Base frame templates (Atomic_*.pdb)

### Core Algorithm Classes

- `BaseFrameCalculator` - Reference frame calculation via LS fitting
- `BasePairFinder` - Greedy pair finding with mutual selection
- `BasePairValidator` - Geometric checks, H-bonds, quality scores
- `ParameterCalculator` - Step and helical parameter calculation
- `HelixOrganizer` - Helix organization and step ordering
- `HydrogenBondFinder` - H-bond detection with donor/acceptor validation

### Recent Changes (December 2024)

**Comparison Targets System**: Generalized comparison to support multiple targets (legacy, baseline, etc.):
- `x3dna_json_compare/targets.py` - ComparisonTarget class and target registry
- `x3dna_json_compare/runner.py` - Updated to accept `--target` parameter
- `x3dna_json_compare/cli.py` - Added `targets` and `baseline` commands
- Supports comparing modern output against any target directory under `data/json_<target>/`

**res_id Comparison Matching**: Updated Python comparison modules to use stable `res_id` identifiers:
- `x3dna_json_compare/step_comparison.py` - Added `build_residue_idx_to_res_id_map()` for legacy res_id construction
- `x3dna_json_compare/base_pair_comparison.py` - Uses res_id for matching base pairs
- `x3dna_json_compare/hbond_comparison.py` - Uses res_id for matching H-bond lists
- `x3dna_json_compare/res_id_utils.py` - Shared utilities for res_id extraction

**Removed ~1,700 LOC of dead/redundant code:**
- HelixDetector (redundant with HelixOrganizer)
- HydrogenBondCounter (consolidated into HydrogenBondFinder)
- Helix submodule files (extracted but never integrated)
- ResidueIndexFixer (unused)

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

- `x3dna_json_compare/cli.py` - fp2-validate CLI entry point
- `x3dna_json_compare/rebuild.py` - JSON regeneration module
- `x3dna_json_compare/runner.py` - Validation runner (auto-regenerates before comparison)
- `scripts/compare_json.py` - Standalone JSON comparison tool
- `docs/CODE_FLOW.md` - Detailed algorithm flow
- `docs/TESTING_GUIDE.md` - Complete testing reference

## H-Bond Prototyping (December 2024)

**Directory**: `prototypes/hbond_optimizer/`

Python prototype for improved H-bond detection with proper saturation tracking. Key features:
- Predicts H atom positions and lone pair directions from known geometry
- Tracks which H slots (for donors) and LP slots (for acceptors) are used
- Greedy best-first selection respecting chemical capacity limits
- Compares against DSSR output for validation

**Key insight**: NH2 groups (N6, N4, N2) can donate 2 H-bonds, but each must use a different hydrogen. sp2 oxygens (O6, O2, O4) can accept 2 H-bonds, but each must use a different lone pair. The prototype computes actual H/LP positions from geometry and tracks which slots are used.

**Design document**: `docs/HBOND_STRATEGIES_SURVEY.md` - Survey of barnaba, ClaRNA, FR3D, HBPLUS approaches

```bash
# Run prototype
cd prototypes/hbond_optimizer
python hbond_optimizer.py data/pdb/1GID.pdb

# Compare with DSSR
python compare_with_dssr.py data/pdb/1GID.pdb
```

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
