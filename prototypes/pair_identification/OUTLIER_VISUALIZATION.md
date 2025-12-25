# Outlier Visualization Tool

Tool for exploring and visualizing cWW base pair classification discrepancies between our system and DSSR.

## Quick Start

```bash
cd prototypes/pair_identification

# List all outliers
python visualize_outlier.py

# Filter by PDB
python visualize_outlier.py --pdb 1EHZ

# Filter by reason
python visualize_outlier.py --reason wrong_atoms

# Combine filters
python visualize_outlier.py --pdb 1EHZ --reason geometric_outlier

# Render specific outlier (by index from list)
python visualize_outlier.py --render 5 --output-dir viz_output

# Render specific pair directly
python visualize_outlier.py --render-pair 1EHZ A-G-1 A-C-72
```

## Analysis Results Summary

**100 PDB Test Set Results:**
- Total cWW pairs (DSSR): 2,521
- True Positives: 2,196 (87.1%)
- False Negatives: 325 (12.9%)
- **False Positives: 0** (we never incorrectly call cWW)

## Reason Codes

| Reason | Count | Description |
|--------|-------|-------------|
| `distance_issues` | 216 (66.5%) | H-bond distances too long or too short |
| `extra_hbonds` | 157 (48.3%) | Unexpected H-bonds detected |
| `wrong_atoms` | 156 (48.0%) | H-bonds to wrong acceptor atoms (e.g., N2→N3 instead of N2→O2) |
| `missing_hbonds` | 53 (16.3%) | Expected H-bonds not found |
| `geometric_outlier` | 16 (4.9%) | N1N9 distance or angle outlier AND RMSD > 0.5Å |
| `non_canonical` | 13 (4.0%) | DSSR Saenger class is "--" |
| `no_hbonds` | 5 (1.5%) | No H-bonds detected at all |
| `rmsd_prefers_other` | 2 (0.6%) | Geometry fits non-cWW template > 0.5Å better |

## Command Reference

### Listing Outliers

```bash
# Basic list (shows first 50)
python visualize_outlier.py

# Show more
python visualize_outlier.py --limit 100

# Filter by PDB
python visualize_outlier.py --pdb 1GID

# Filter by reason
python visualize_outlier.py --reason wrong_atoms
python visualize_outlier.py --reason geometric_outlier
python visualize_outlier.py --reason distance_issues

# Combine filters
python visualize_outlier.py --pdb 1DDY --reason wrong_atoms
```

### Rendering Visualizations

```bash
# Render by index (from filtered list)
python visualize_outlier.py --pdb 1EHZ --render 0

# Render specific pair directly
python visualize_outlier.py --render-pair 1EHZ A-G-1 A-C-72

# Custom output directory
python visualize_outlier.py --render 5 --output-dir my_visualizations

# Launch PyMOL automatically
python visualize_outlier.py --render 5 --launch
```

### Custom Paths

```bash
python visualize_outlier.py \
  --results analysis_results/cww_analysis_100 \
  --pdb-dir /path/to/pdbs \
  --template-dir /path/to/basepair-idealized \
  --hbond-dir /path/to/slot_hbonds
```

## Output Files

Each visualization generates:
- `target_<PDB>_<SEQ>.pdb` - The actual base pair from the structure
- `aligned_<PDB>_<LW>_<SEQ>.pdb` - Templates aligned to the target
- `view_<PDB>_<RES1>_<RES2>.pml` - PyMOL script

## Individual Report Format

Each PDB has a JSON report in `analysis_results/cww_analysis_100/`:

```json
{
  "pdb_id": "1EHZ",
  "summary": {
    "total_canonical_cww": 17,
    "true_positives": 11,
    "false_negatives": 6,
    "false_positives": 0
  },
  "false_negatives": [
    {
      "res_id1": "A-A-5",
      "res_id2": "A-U-68",
      "sequence": "AU",
      "reasons": ["geometric_outlier"],
      "hbond_diagnostics": {
        "expected_hbonds": [...],
        "found_hbonds": [...],
        "missing_hbonds": [],
        "wrong_atoms": {}
      },
      "geometric_diagnostics": {
        "rmsd_cww": 0.196,
        "best_lw": "cWW",
        "interbase_angle": 17.7,
        "n1n9_distance": 8.84,
        "is_geometric_outlier": true
      }
    }
  ]
}
```

## Interpreting Results

### Geometric Outlier
The pair has unusual geometry:
- **N1N9 distance**: Distance between glycosidic nitrogens. Expected ~8.7-9.2Å for standard WC pairs.
- **Interbase angle**: Angle between base planes. Expected <15° for cWW.

### Wrong Atoms
H-bonds detected to unexpected atoms:
- `N2→N3` instead of `N2→O2` (GC pair)
- `N1→O2` instead of `N1→N3` (GC pair)

### Distance Issues
H-bond distances outside expected range (2.5-3.5Å).

### RMSD Prefers Other
The pair geometry fits a non-cWW template (tWW, cWS, etc.) better than cWW.
