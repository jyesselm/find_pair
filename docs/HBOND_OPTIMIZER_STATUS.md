# H-Bond Optimizer Status

## Current Performance (December 2024)

**Fast PDB Set (3602 PDBs):**
- **Recall: 95.45%** (303,136 / 317,598 DSSR H-bonds matched)
- **Precision: 91.47%** (303,136 / 331,408 optimizer H-bonds correct)

## Miss Categories

| Reason | Count | Percentage | Description |
|--------|-------|------------|-------------|
| POOR_ALIGNMENT | 6,529 | 45.1% | Geometric alignment score < 0.3 threshold (for dist ≥ 3.2Å) |
| MISSING_RESIDUE | 5,388 | 37.3% | Residue not found in our PDB parsing |
| MISSING_ATOM | 1,913 | 13.2% | Atom not found in residue |
| SLOT_CONFLICT | 501 | 3.5% | Lost in greedy selection (slot saturated) |
| NO_LP_SLOTS | 70 | 0.5% | Acceptor has no lone pair slots defined |
| NOT_DONOR_ACCEPTOR | 56 | 0.4% | Atom pair not in donor/acceptor tables |
| NO_H_SLOTS | 5 | 0.0% | Donor has no hydrogen slots defined |

## Top Missed Atom Pairs

| Atom Pair | Count | Notes |
|-----------|-------|-------|
| N2-O2 | 3,458 | G.N2 donor to C/U.O2 acceptor |
| N1-N3 | 2,437 | Watson-Crick edge |
| N4-O6 | 1,864 | C.N4 donor to G.O6 acceptor |
| N6-O2 | 1,492 | A.N6 donor to C/U.O2 acceptor |
| N6-O4 | 1,440 | A.N6 donor to U.O4 acceptor |
| N6-O6 | 1,297 | A-G mismatch type |
| N2-N3 | 1,051 | G-G or G-A type |

## Key Issues Identified

### 1. POOR_ALIGNMENT (45.1%)
Geometric prediction of H/LP directions doesn't match the actual H-bond geometry for longer-distance bonds. Now uses **short distance threshold** (3.2Å) - bonds below this distance bypass alignment check.

**Statistics (for dist ≥ 3.2Å):**
- Min alignment: -4.000
- Max alignment: 0.300
- Mean alignment: -0.843
- Median alignment: -0.139

**Key insight**: For very short distances (< 3.2Å), the geometry is less important because the atoms are so close that H-bonding is highly likely. For longer distances, alignment matters more.

**LP Geometry Model:**
- Carbonyl oxygens (O2, O4, O6): sp2 with two LPs at 120° from C=O bond, in base plane
- Ring nitrogens (N1, N3, N7): sp2 with one LP pointing out of ring
- Ribose oxygens (O2', O4'): sp3 with flexible geometry
- Phosphate oxygens (OP1, OP2): Tetrahedral with 3 LP directions

### 2. MISSING_RESIDUE (37.3%)
Residues reported by DSSR that we don't find in our PDB parsing. Top examples:
- Standard residues: A.A231, A.U238, N.DG31, N.DG33 - likely chain/numbering mismatches
- Unusual chains: X.FMN200 (flavin mononucleotide) - non-nucleotide
- Insertion codes: B.G47^A - unusual format

### 3. Modified Residues Needing Attention

**NOT_DONOR_ACCEPTOR (unusual modified residues):**
- 6AP (6-aminopurine) - N2 donor not in our tables
- XAN (xanthosine) - O2 acceptor not in our tables
- NMN (nicotinamide) - N7 acceptor not in our tables
- AZA (8-azaguanine) - O2 acceptor

**NO_LP_SLOTS (missing phosphate acceptors):**
- 2BA, 6MZ, NMN, EPE, AMP - OP1/OP2 need LP slot definitions

**NO_H_SLOTS:**
- U.N3 - 4 cases where N3 has no H slots (should have 1)
- EPE.N4 - 1 case

## Architecture

### Files
- `prototypes/hbond_optimizer/optimizer.py` - Main HBondOptimizer class with greedy selection
- `prototypes/hbond_optimizer/geometry.py` - H/LP slot prediction, alignment scoring
- `prototypes/hbond_optimizer/compare_with_dssr.py` - DSSR comparison utilities
- `prototypes/hbond_optimizer/modified_registry.py` - Modified residue → parent base mapping
- `prototypes/hbond_optimizer/benchmark.py` - Multithreaded benchmark with detailed miss analysis

### Key Concepts
1. **Slot-based model**: Each donor has H slots, each acceptor has LP slots
2. **Bifurcation support**: One H/LP can serve two bonds if angularly separated (≥45°)
3. **Greedy selection**: Sort candidates by distance, accept if alignment ≥ threshold
4. **Modified residue support**: Map 421+ modified codes to parent base types

## Next Steps

### High Priority
1. **Fix POOR_ALIGNMENT** - The LP prediction for carbonyl oxygens (O2, O4, O6) often points perpendicular to the base plane, but many H-bonds come from in-plane. Consider:
   - Lowering alignment threshold for short distances
   - Adding alternative LP directions for carbonyl acceptors
   - Using distance-only mode as fallback

2. **Add phosphate LP slots** - OP1/OP2 are common acceptors but have no LP slots defined

### Medium Priority
3. **Add unusual modified residues** - 6AP, XAN, NMN, AZA, etc.
4. **Investigate MISSING_RESIDUE** - Many are standard bases; check chain/number matching

### Low Priority
5. **Fix U.N3 H slots** - Should have 1 H slot but sometimes returns 0
6. **Handle non-nucleotide residues** - FMN, EPE, etc.

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
