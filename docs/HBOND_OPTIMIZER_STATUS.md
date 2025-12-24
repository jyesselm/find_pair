# H-Bond Optimizer Status

## Current Performance (December 2024)

**Fast PDB Set (3602 PDBs):**
- **Recall: 96.68%** (307,060 / 317,592 DSSR H-bonds matched)
- **Precision: 91.50%** (307,060 / 335,572 optimizer H-bonds correct)

## Miss Categories

| Reason | Count | Percentage | Description |
|--------|-------|------------|-------------|
| POOR_ALIGNMENT | 6,241 | 59.3% | Geometric alignment score < 0.3 threshold (for dist >= 3.2A) |
| MISSING_RESIDUE | 2,157 | 20.5% | Residue not found in our PDB parsing |
| MISSING_ATOM | 1,765 | 16.8% | Atom not found in residue |
| SLOT_CONFLICT | 250 | 2.4% | Lost in greedy selection (slot saturated) |
| NO_LP_SLOTS | 62 | 0.6% | Acceptor has no lone pair slots defined |
| NOT_DONOR_ACCEPTOR | 56 | 0.5% | Atom pair not in donor/acceptor tables |
| NO_H_SLOTS | 1 | 0.0% | Donor has no hydrogen slots defined |

## Recent Improvements

### ID Normalization Fix (Dec 2024)
Fixed residue ID normalization to properly handle:
- **DNA prefixes**: "DG38" now correctly normalizes to "G38" (was incorrectly parsed as "G8" via modified residue "DG3")
- **Standard bases with numbers**: "A231" correctly stays as "A231" (was incorrectly parsed as "A1" via modified residue "A23")
- **Priority order**: Standard bases (1 char) checked BEFORE modified residue registry to avoid false matches

This fix improved recall from 90.21% to **96.68%** (+6.47%).

## Key Issues Identified

### 1. POOR_ALIGNMENT (59.3%)
Geometric prediction of H/LP directions doesn't match the actual H-bond geometry for longer-distance bonds. Uses **short distance threshold** (3.2A) - bonds below this distance bypass alignment check.

**Statistics (for dist >= 3.2A):**
- Min alignment: -4.000
- Max alignment: 0.300
- Mean alignment: -0.842
- Median alignment: -0.124

**LP Geometry Model:**
- Carbonyl oxygens (O2, O4, O6): sp2 with two LPs at 120 deg from C=O bond, in base plane
- Ring nitrogens (N1, N3, N7): sp2 with one LP pointing out of ring
- Ribose oxygens (O2', O4'): sp3 with flexible geometry
- Phosphate oxygens (OP1, OP2): Isotropic with 3 orthogonal LP directions

### 2. MISSING_RESIDUE (20.5%)
Residues reported by DSSR that we don't find in our PDB parsing. Remaining cases:
- Non-nucleotides: X.FMN200 (flavin mononucleotide)
- Unusual insertion codes: B.G47^A, B.C47^F, etc.
- Unusual numbering: R.C2E1, T.G2P3000

### 3. Modified Residues Needing Attention

**NOT_DONOR_ACCEPTOR (56 cases):**
- 6AP (6-aminopurine) - N2 donor not in our tables
- XAN (xanthosine) - O2 acceptor not in our tables
- NMN (nicotinamide) - N7 acceptor not in our tables
- AZA (8-azaguanine) - O2 acceptor

**NO_LP_SLOTS (62 cases):**
- 2BA, 6MZ, NMN, EPE, AMP, CFL, 5GP - OP1/OP2 need LP slot definitions

## Architecture

### Files
- `prototypes/hbond_optimizer/optimizer.py` - Main HBondOptimizer class with greedy selection
- `prototypes/hbond_optimizer/geometry.py` - H/LP slot prediction, alignment scoring
- `prototypes/hbond_optimizer/compare_with_dssr.py` - DSSR comparison utilities
- `prototypes/hbond_optimizer/modified_registry.py` - Modified residue -> parent base mapping
- `prototypes/hbond_optimizer/benchmark.py` - Multithreaded benchmark with detailed miss analysis

### Key Concepts
1. **Slot-based model**: Each donor has H slots, each acceptor has LP slots
2. **Bifurcation support**: One H/LP can serve two bonds if angularly separated (>=45 deg)
3. **Greedy selection**: Sort candidates by distance, accept if alignment >= threshold
4. **Short distance bypass**: Bonds < 3.2A skip alignment check (distance-only mode)
5. **Modified residue support**: Map 421+ modified codes to parent base types

## Next Steps

### High Priority
1. **Reduce POOR_ALIGNMENT** - Consider:
   - Lowering alignment threshold further
   - Adding alternative LP directions for carbonyl acceptors
   - Using adaptive thresholds based on distance

### Medium Priority
2. **Add unusual modified residues** - 6AP, XAN, NMN, AZA donor/acceptor entries
3. **Fix NO_LP_SLOTS for modified residues** - 2BA, 6MZ, NMN, EPE, etc.

### Low Priority
4. **Handle non-nucleotide residues** - FMN, EPE ligands
5. **Handle unusual insertion codes** - B.G47^A format

## Running Benchmarks

```bash
# Quick test (100 PDBs)
python3 prototypes/hbond_optimizer/benchmark.py --test-set 100 --workers 8

# Full benchmark (3602 PDBs, ~5 minutes with 10 workers)
python3 prototypes/hbond_optimizer/benchmark.py --test-set fast --workers 10 --output data/hbond_benchmark_fast.json

# Custom parameters
python3 prototypes/hbond_optimizer/benchmark.py --test-set 100 --max-dist 3.5 --min-align 0.2
```

The detailed JSON report includes all missed H-bonds with original DSSR residue IDs for debugging modified residue issues.
