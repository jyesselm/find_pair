# Pair Identification Prototype Architecture

This document explains the flow, components, and data sources of the pair identification prototype system.

## Overview

The pair identification prototype implements a **template-based Leontis-Westhof (LW) base pair classification** system that:

1. Loads pre-computed reference frames from modern C++ JSON output
2. Validates pairs geometrically (distance, angle, planarity)
3. Classifies pairs using template RMSD + H-bond pattern matching
4. Compares results against DSSR as ground truth

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         DATA SOURCES                                    │
├─────────────────────────────────────────────────────────────────────────┤
│  data/pdb/*.pdb           → Raw PDB structures                          │
│  data/json/ls_fitting/    → Reference frames (from C++ FindPairProtocol)│
│  data/json/pdb_atoms/     → Extracted atom coordinates                  │
│  data/json/all_hbond_list/→ H-bonds detected by modern C++              │
│  data/json_dssr/          → DSSR output (ground truth)                  │
│  basepair-idealized/      → Idealized template coordinates              │
│  basepair-exemplars/      → Exemplar template coordinates               │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                      DATA LOADING LAYER                                 │
├─────────────────────────────────────────────────────────────────────────┤
│  frame_loader.py   → Loads ReferenceFrame objects from JSON             │
│  pair_cache.py     → Orchestrates loading, KDTree neighbor finding      │
│  hbond_loader.py   → Loads H-bonds from modern + DSSR                   │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                     VALIDATION LAYER                                    │
├─────────────────────────────────────────────────────────────────────────┤
│  geometric_validator.py → Frame-based geometric checks                  │
│  cww_validator.py       → cWW-specific H-bond pattern validation        │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                    CLASSIFICATION LAYER                                 │
├─────────────────────────────────────────────────────────────────────────┤
│  template_generator.py → Builds template registry from DSSR + idealized │
│  template_aligner.py   → Kabsch superposition, RMSD calculation         │
│  lw_classifier.py      → Combined RMSD (40%) + H-bond (60%) scoring     │
│  hbond_scorer.py       → Pattern matching against expected H-bonds      │
│  template_overlay.py   → PyMOL visualization generation                 │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                    VALIDATION SCRIPTS                                   │
├─────────────────────────────────────────────────────────────────────────┤
│  validate_classifier.py         → Basic classifier validation           │
│  validate_combined_classifier.py→ Combined RMSD+H-bond validation       │
│  validate_cww_parallel.py       → Parallel cWW validation               │
│  analyze_cww_misses.py          → Root cause analysis of misses         │
│  analyze_standard_wc.py         → Standard WC pair analysis             │
│  compare_cww.py                 → DSSR vs Legacy vs Modern comparison   │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Data Loading Flow

### 1. Frame Loading (`frame_loader.py`)

**Source**: `data/json/ls_fitting/{PDB_ID}.json`

**Key Class**: `ReferenceFrame`
- `origin`: 3D translation vector (numpy array)
- `rotation`: 3x3 rotation matrix (axes as columns)
- `rmsd_fit`: RMSD quality metric from least-squares fitting

**Output**: Dictionary mapping `res_id → ReferenceFrame`

```python
loader = FrameLoader(Path("data/json"))
frames = loader.load_frames("1EHZ")  # Returns Dict[str, ReferenceFrame]
```

### 2. Pair Cache (`pair_cache.py`)

**Sources**:
- Frames: `data/json/ls_fitting/{PDB_ID}.json`
- Atoms: `data/json/pdb_atoms/{PDB_ID}.json`

**Key Classes**:
- `AtomCoords`: Stores C1', N1/N9, ring atoms, H-bond atoms for a residue
- `CachedPair`: Complete data for a validated base pair
- `PairCache`: Main orchestrator

**Workflow**:
```
PairCache.build_cache()
  → Load frames (FrameLoader)
  → Load atoms (pdb_atoms JSON)
  → Build KDTree on frame origins
  → For each residue pair within max_distance (15Å):
      → Validate geometry (GeometricValidator)
      → Create CachedPair
      → Store in self.pairs
```

### 3. H-Bond Loading (`hbond_loader.py`)

**Sources**:
- Modern: `data/json/all_hbond_list/{PDB_ID}.json`
- DSSR: `data/json_dssr/{PDB_ID}.json`

**Key Classes**:
- `HBond`: Single H-bond (donor, acceptor, distance)
- `PairHBonds`: Collection of H-bonds between two residues

**Functions**:
```python
modern_hbonds = load_modern_hbonds(Path("data/json/all_hbond_list/1EHZ.json"))
dssr_hbonds = load_dssr_hbonds(Path("data/json_dssr/1EHZ.json"))
comparison = compare_hbonds(modern_hbonds, dssr_hbonds)
```

---

## Validation Pipeline

### Stage 1: Geometric Validation (`geometric_validator.py`)

Validates pairs based on reference frame geometry:

| Check | Threshold | Description |
|-------|-----------|-------------|
| `dorg` | ≤ 15.0 Å | Origin distance between frames |
| `d_v` | ≤ 2.5 Å | Vertical distance (along helix axis) |
| `plane_angle` | ≤ 65° | Angle between base normals |
| `dNN` | ≥ 4.5 Å | N1/N9 glycosidic nitrogen distance |

**Quality Score**: `dorg + 1.5*d_v + plane_angle/180`

**Result**: `ValidationResult` with all metrics + `is_valid` flag

### Stage 2: cWW Chemistry Validation (`cww_validator.py`)

For Watson-Crick pairs specifically:

| Check | Threshold | Description |
|-------|-----------|-------------|
| N1N9 distance | 7.5-10.5 Å | Expected cWW range |
| H-bond count | ≥ 2 | Minimum H-bonds required |
| H-bond pattern | Match expected | GC: 3 bonds, AU: 2 bonds |

**Expected H-Bond Patterns**:
```
GC: O6↔N4, N1↔N3, N2↔O2 (3 H-bonds)
AU: N6↔O4, N1↔N3         (2 H-bonds)
GU: O6↔N3, N1↔O2         (2 H-bonds, wobble)
```

**Quality Score**: `0.6 * hbond_match_ratio + 0.4 * n1n9_score`

---

## Classification Pipeline

### Template Generation (`template_generator.py`)

**Sources**:
1. DSSR JSON → Statistical properties (mean/std for N1N9, angle, planarity)
2. Idealized PDBs → Atomic coordinates for superposition

**Output**: Template registry `(sequence, lw_class) → IdealizedTemplate`

### Template Alignment (`template_aligner.py`)

**Algorithm**: Kabsch optimal superposition
1. Center both point sets
2. Compute covariance matrix H = P^T @ Q
3. SVD: U, S, Vt = SVD(H)
4. Rotation: R = Vt^T @ U^T
5. Handle reflections (det(R) < 0)

**Atoms Used**: Ring atoms (C2, C4, C5, C6, N1, N3, N7, C8, N9)

### LW Classification (`lw_classifier.py`)

**12 LW Classes Tested**: cWW, tWW, cWH, tWH, cWS, tWS, cHH, tHH, cHS, tHS, cSS, tSS

**Combined Scoring**:
```
rmsd_score = 1.0 - (rmsd / 2.0)  if rmsd < 2.0Å, else 0.0
hbond_score = matched_hbonds / expected_hbonds

combined_score = 0.6 * hbond_score + 0.4 * rmsd_score + 0.05 * hbond_count
confidence = min(1.0, gap_to_second_best / 0.3)
```

### H-Bond Pattern Scoring (`hbond_scorer.py`)

**Pattern Database**: Expected donor-acceptor pairs for each (LW class, sequence)

```python
HBOND_PATTERNS = {
    "cWW": {
        "GC": [("N1", "N3"), ("N2", "O2"), ("O6", "N4")],
        "AU": [("N6", "O4"), ("N1", "N3")],
        "GU": [("N1", "O2")],
    },
    "tWW": {
        "GC": [("N2", "O2"), ("N1", "N3")],
    },
    # ... more patterns
}
```

**Scoring**: `score = matched_count / expected_count`

---

## Validation Scripts

### Quick Reference

| Script | Purpose | Key Flags |
|--------|---------|-----------|
| `validate_classifier.py` | Basic N1N9+RMSD validation | `--max-pdbs`, `--lw-classes` |
| `validate_combined_classifier.py` | RMSD+H-bond combined | `--lw-filter`, `-v` |
| `validate_cww_parallel.py` | Parallel cWW validation | `--workers`, `--pdb-list` |
| `analyze_cww_misses.py` | Miss root cause analysis | `--output` (JSON) |
| `analyze_standard_wc.py` | GC/CG/AU/UA analysis | `--show-misses` |
| `compare_cww.py` | DSSR vs Legacy vs Modern | `--analyze-dssr-only` |

### Example Commands

```bash
# Validate combined classifier on 50 PDBs
python validate_combined_classifier.py \
  --pdb-dir data/pdb \
  --hbond-dir data/json/all_hbond_list \
  --dssr-dir data/json_dssr \
  --max-pdbs 50 -v

# Parallel cWW validation with 10 workers
python validate_cww_parallel.py \
  --pdb-dir data/pdb \
  --hbond-dir data/json/all_hbond_list \
  --dssr-dir data/json_dssr \
  --workers 10 \
  --max-pdbs 500

# Analyze misses and export to JSON
python analyze_cww_misses.py \
  --max-pdbs 500 \
  --output misses_analysis.json

# Compare coverage across systems
python compare_cww.py \
  --dssr-dir data/json_dssr \
  --legacy-dir data/json_legacy/base_pair \
  --modern-dir data/json/base_pair \
  --max-pdbs 100 -v
```

### Miss Reason Categories

When classification fails, misses are categorized:

| Reason | Description |
|--------|-------------|
| `no_hbonds` | No H-bonds detected for the pair |
| `pattern_mismatch` | H-bonds found but wrong pattern |
| `rmsd_prefers_other` | Alternative LW class has better RMSD |
| `partial_match` | Partial H-bond pattern match only |

---

## Example Usage

### Building a Pair Cache

```python
from pair_cache import PairCache
from pathlib import Path

cache = PairCache("1EHZ", Path("data/json"))
cache.build_cache(max_distance=15.0)

# Get valid pairs
valid_pairs = cache.get_valid_pairs()
print(f"Found {len(valid_pairs)} valid pairs")

# Save for later use
cache.save(Path("pair_cache_1EHZ.json"))
```

### Classifying a Pair

```python
from lw_classifier import LWClassifier
from template_aligner import TemplateAligner

aligner = TemplateAligner(
    idealized_dir=Path("basepair-idealized"),
    exemplar_dir=Path("basepair-exemplars")
)

classifier = LWClassifier(aligner)

result = classifier.classify(
    res1_atoms=residue1.ring_atoms,
    res2_atoms=residue2.ring_atoms,
    observed_hbonds=pair_hbonds
)

print(f"Best match: {result.best_lw} (confidence: {result.confidence:.2f})")
```

### Comparing with DSSR

```python
from hbond_loader import load_modern_hbonds, load_dssr_hbonds, compare_hbonds

modern = load_modern_hbonds(Path("data/json/all_hbond_list/1EHZ.json"))
dssr = load_dssr_hbonds(Path("data/json_dssr/1EHZ.json"))

comparison = compare_hbonds(modern, dssr)
print(f"Common pairs: {comparison['common_pairs']}")
print(f"Modern-only: {comparison['modern_only_pairs']}")
print(f"DSSR-only: {comparison['dssr_only_pairs']}")
```

---

## Key Design Decisions

1. **res_id Format**: All modules use `chain-name-seq[ins]` (e.g., "A-G-1", "A-C-10A") for stable residue identification

2. **Two-Stage Validation**: Geometric checks (fast, frame-based) before chemical validation (slow, atom-based)

3. **Combined Scoring**: 60% H-bond pattern + 40% RMSD provides robust classification

4. **Confidence Metric**: Gap between best and second-best scores indicates classification certainty

5. **Parallel Processing**: `validate_cww_parallel.py` uses multiprocessing for large-scale validation

6. **Bidirectional H-Bond Matching**: Handles donor/acceptor labeling inconsistencies

---

## File Dependencies

```
pair_cache.py
  └── frame_loader.py
  └── geometric_validator.py

lw_classifier.py
  └── template_aligner.py
  └── hbond_scorer.py
  └── template_overlay.py (for visualization)

validate_*.py scripts
  └── lw_classifier.py
  └── hbond_loader.py
  └── template_aligner.py
```
