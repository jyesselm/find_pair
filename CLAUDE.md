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
# Required for legacy code execution
export X3DNA=/Users/jyesselman2/local/installs/x3dna

# Install Python package
pip install -e ".[dev]"
```

## Testing

Primary testing is JSON regression comparison - modern output must match legacy output exactly.

```bash
# Main validation CLI (always use fp2-validate)
X3DNA=/Users/jyesselman2/local/installs/x3dna python -m x3dna_json_compare.cli validate all --test-set 100
X3DNA=/Users/jyesselman2/local/installs/x3dna python -m x3dna_json_compare.cli validate frames --test-set 100
X3DNA=/Users/jyesselman2/local/installs/x3dna python -m x3dna_json_compare.cli validate 3 --pdb 1EHZ -v

# Run pytest
pytest tests_python/

# Generate JSON for comparison
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=frames

# Regenerate JSON files (always verify after generation)
PYTHONPATH=. X3DNA=/Users/jyesselman2/local/installs/x3dna python3 scripts/rebuild_json.py regenerate --test-set 100
```

### Current Test Results (100-PDB Test Set)

| Stage | Pass Rate | Notes |
|-------|-----------|-------|
| 1: Atoms | 94% (94/100) | 6 failures: 2CV2, 4RQF, 5VOE, 8UPT, 8Z1P, 8ZYC |
| 2: Residue Indices | 100% (100/100) | All pass |
| 3-5: Frames | 99% (99/100) | 1 failure: 6OZK (base_type `g` vs `I` for inosine) |
| 6-7: Pairs | 99% (99/100) | Only 4RQF fails (corrupt legacy JSON) |
| 11-12: Steps | 95% (95/100) | 5 failures: 3UCU, 5CCX, 7YGB, 8RUJ, 8U5Z |

**Note**: Step ordering (`bp_idx`) now matches legacy 100% when helices match. Fixes implemented:
- Fixed `end_stack_xang` threshold: 110° → 125° (matches legacy `END_STACK_XANG`)
- Fixed backbone extraction to use `legacy_residue_idx` from atoms
- Added second `check_direction` call after `check_strand2` (legacy line 1361)
- Fixed neighbor swapping in `calculate_context` (legacy lines 931-941)
- Added Watson-Crick pair check in `wc_bporien` (legacy `base_pairs[m][3] > 0`)
- Fixed modified base handling in wc_bporien ("Cg" → "CG")
- Added geometric bpid check (dir_x > 0 && dir_y < 0 && dir_z < 0)

**Remaining 10 failures** (2EEW, 3UCU, 5FJ1, 5Y85, 6ICZ, 7YGA, 7YGB, 8RUJ, 8U5Z, 8Z1P) have different helix groupings. See `docs/STEP_PARAMETER_INVESTIGATION.md` for details.

### Validation Stages (must pass in order)

| Stage | Name | Test Command |
|-------|------|--------------|
| 1 | Atoms | `fp2-validate validate 1` |
| 2 | Residue Indices | `fp2-validate validate 2` |
| 3-5 | Frames | `fp2-validate validate frames` |
| 6-7 | Pair Validation | `fp2-validate validate pairs` |
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
