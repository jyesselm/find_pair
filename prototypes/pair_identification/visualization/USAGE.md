# PairVisualizer Usage Guide

Complete guide to using the pair visualization module for base pair analysis.

## Quick Start

### Method 1: Command Line (Easiest)
```bash
cd prototypes/pair_identification/visualization

# Basic usage
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72

# Custom LW classes
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72 --lw cWW,tWW

# Launch PyMOL automatically
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72 --launch-pymol
```

### Method 2: Python API
```python
from pathlib import Path
from visualization import PairVisualizer

# Initialize
viz = PairVisualizer(
    pdb_dir=Path("data/pdb"),
    template_dir=Path("basepair-idealized"),
    output_dir=Path("viz_output"),
)

# Visualize a pair
result = viz.visualize_pair(
    pdb_id="1EHZ",
    res_id1="A-G-1",
    res_id2="A-C-72",
)

# View results
print(f"Best match: {result.best_lw} (RMSD={result.best_rmsd:.3f}Å)")
print(f"PyMOL script: {result.pymol_script}")
```

### Method 3: Example Script
```bash
cd prototypes/pair_identification/visualization
python example_usage.py
```

## Common Use Cases

### 1. Investigate a cWW Miss
Why didn't this pair classify as cWW?

```bash
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72 --lw cWW,tWW,cWS
```

Check:
- RMSD to cWW template (should be < 1.5Å)
- Missing H-bonds (red dashes)
- Structural deviations (overlay with templates)

### 2. Compare Multiple Templates
Which LW class best fits this pair?

```python
result = viz.visualize_pair(
    pdb_id="1EHZ",
    res_id1="A-G-5",
    res_id2="A-U-68",
    lw_classes=["cWW", "tWW", "cWS", "tWH", "tHS"],
)

# Print sorted by RMSD
for lw, rmsd in sorted(result.rmsd_values.items(), key=lambda x: x[1]):
    print(f"{lw}: {rmsd:.4f} Å")
```

### 3. H-Bond Validation
See which expected H-bonds are missing:

```bash
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72
```

In PyMOL:
- Yellow dashes = Observed (detected)
- Green dashes = Expected and found
- Red dashes = Expected but missing

### 4. Batch Processing
Visualize multiple pairs:

```python
pairs = [
    ("1EHZ", "A-G-1", "A-C-72"),
    ("1EHZ", "A-G-2", "A-C-71"),
    ("1EHZ", "A-G-3", "A-C-70"),
]

for pdb_id, res1, res2 in pairs:
    result = viz.visualize_pair(pdb_id, res1, res2)
    if result:
        print(f"{res1}-{res2}: {result.best_lw} ({result.best_rmsd:.3f}Å)")
```

## CLI Options Reference

```bash
python cli.py --help

Required:
  --pdb PDB_ID          PDB identifier (e.g., 1EHZ)
  --res1 RES_ID         First residue (e.g., A-G-1)
  --res2 RES_ID         Second residue (e.g., A-C-72)

Optional:
  --lw LW_CLASSES       Comma-separated LW classes (default: cWW,tWW,cWS,tHS)
  --pdb-dir PATH        PDB directory (default: data/pdb)
  --template-dir PATH   Template directory (default: basepair-idealized)
  --hbond-dir PATH      H-bond JSON directory (default: data/json/slot_hbonds)
  --output-dir PATH     Output directory (default: viz_output)
  --no-hbonds           Skip H-bond loading
  --launch-pymol        Launch PyMOL automatically
```

## Python API Reference

### PairVisualizer

```python
class PairVisualizer:
    def __init__(
        self,
        pdb_dir: Path = Path("data/pdb"),
        template_dir: Path = Path("basepair-idealized"),
        hbond_dir: Path = Path("data/json/slot_hbonds"),
        output_dir: Path = Path("viz_output"),
    )
```

**Methods:**

```python
def visualize_pair(
    self,
    pdb_id: str,
    res_id1: str,
    res_id2: str,
    lw_classes: Optional[List[str]] = None,
    include_hbonds: bool = True,
) -> Optional[VisualizationResult]
```

**Returns:** `VisualizationResult` or `None` if failed

### VisualizationResult

```python
@dataclass
class VisualizationResult:
    pdb_id: str
    res_id1: str
    res_id2: str
    sequence: str              # e.g., "GC"

    target_pdb: Path           # Target pair PDB
    template_pdbs: Dict[str, Path]  # LW class -> aligned PDB
    pymol_script: Path         # Generated .pml script

    rmsd_values: Dict[str, float]   # LW class -> RMSD
    best_lw: str               # Best matching LW class
    best_rmsd: float           # Minimum RMSD
```

## PyMOL Script Features

The generated `.pml` scripts include:

### Visual Elements
- Target pair (cyan, thick sticks)
- Aligned templates (colored, thin sticks)
- Observed H-bonds (yellow dashes)
- Expected H-bonds (green/red dashes)
- Atom labels (N1, N2, N3, etc.)

### Interactive Commands
```pymol
# Toggle individual templates
disable cWW
enable cWW

# Toggle all templates
disable aligned_*

# Toggle H-bonds
disable obs_*
disable exp_*

# Focus on specific atoms
zoom name N1+N3+O6
select hbond_atoms, name N1+N2+N3+N4+N6+N7+O2+O4+O6
```

### Color Scheme
- **Target**: cyan
- **cWW**: green
- **tWW**: orange
- **cWS**: magenta
- **tWH**: yellow
- **tHS**: pink
- **Observed H-bonds**: yellow
- **Expected (found)**: green
- **Expected (missing)**: red

## Interpretation Guide

### RMSD Values

| RMSD (Å) | Interpretation |
|----------|----------------|
| < 0.5    | Excellent match |
| 0.5-1.0  | Good match |
| 1.0-1.5  | Fair match |
| 1.5-2.0  | Poor match |
| > 2.0    | Very poor match |

### Classification Decision

A pair should be classified as LW class X if:
1. RMSD to X template < 1.5Å
2. Expected H-bonds for X are present (>= 50%)
3. RMSD to X is significantly better than other classes (> 0.5Å difference)

### Common Issues

**High RMSD to all templates:**
- Non-canonical pair
- Distorted geometry
- Misidentified pair

**Missing H-bonds:**
- Poor geometry
- Water-mediated interactions
- Modified bases

**Multiple good matches:**
- Intermediate geometry
- Need secondary criteria (H-bonds, context)

## Troubleshooting

### Error: "PDB not found"
Check `--pdb-dir` path and ensure PDB file exists:
```bash
ls data/pdb/1EHZ.pdb
```

### Error: "Residues not found"
Verify res_id format: `chain-base-seq` (e.g., "A-G-1")
```python
from core.pdb_parser import parse_pdb
residues = parse_pdb(Path("data/pdb/1EHZ.pdb"))
print(list(residues.keys())[:10])  # Show first 10 res_ids
```

### Error: "Template not found"
Check template directory structure:
```bash
ls basepair-idealized/cWW/GC.pdb
```

### No H-bonds shown
Check if H-bond JSON exists:
```bash
ls data/json/slot_hbonds/1EHZ.json
```

If missing, generate with:
```bash
python generate_slot_hbonds.py data/pdb/1EHZ.pdb
```

### PyMOL won't launch
Ensure PyMOL is in PATH:
```bash
which pymol
```

Or use manual launch:
```bash
pymol viz_output/view_*.pml
```

## Advanced Usage

### Custom Template Colors
```python
viz = PairVisualizer()
viz.TEMPLATE_COLORS["cWW"] = "blue"
viz.TEMPLATE_COLORS["tWW"] = "red"
```

### Custom Ring Atoms
```python
viz = PairVisualizer()
viz.RING_ATOMS = ["C2", "C4", "C5", "C6", "N1", "N3"]  # Pyrimidine only
```

### Generate Script Only (No Alignment)
```python
from visualization import PyMOLScriptGenerator, TemplateViz

gen = PyMOLScriptGenerator(output_dir=Path("viz_output"))

script = gen.generate_pair_script(
    target_pdb=Path("target.pdb"),
    templates=[
        TemplateViz("cWW", Path("cww.pdb"), 0.5, "green", True),
    ],
    observed_hbonds=[],
    expected_hbonds=[],
    output_name="custom",
    title="Custom Visualization",
)
```

## Integration with Other Tools

### 1. With CWW Annotator
```bash
# Generate annotations
python cww_annotate.py 1EHZ > annotations.txt

# Visualize failures
python cli.py --pdb 1EHZ --res1 <from_annotations> --res2 <from_annotations>
```

### 2. With LW Classifier
```python
from lw_classifier import LWClassifier
from visualization import PairVisualizer

classifier = LWClassifier()
viz = PairVisualizer()

# Classify and visualize
lw_class = classifier.classify(pdb_id, res_id1, res_id2)
result = viz.visualize_pair(pdb_id, res_id1, res_id2)

print(f"Classifier: {lw_class}")
print(f"Best RMSD: {result.best_lw}")
```

### 3. With Validation Tools
```python
from geometric_validator import GeometricValidator
from visualization import PairVisualizer

validator = GeometricValidator()
viz = PairVisualizer()

# Check validation
is_valid, metrics = validator.validate(pdb_id, res_id1, res_id2)

# Visualize if invalid
if not is_valid:
    result = viz.visualize_pair(pdb_id, res_id1, res_id2)
```

## File Organization

Recommended directory structure:
```
project/
├── data/
│   ├── pdb/                    # Input PDBs
│   └── json/
│       └── slot_hbonds/        # H-bond data
├── basepair-idealized/         # Templates
│   ├── cWW/
│   ├── tWW/
│   └── ...
└── viz_output/                 # Generated files
    ├── target_*.pdb
    ├── aligned_*.pdb
    └── view_*.pml
```

## Performance Tips

1. **Skip H-bonds** for faster generation:
   ```python
   result = viz.visualize_pair(..., include_hbonds=False)
   ```

2. **Limit LW classes**:
   ```python
   result = viz.visualize_pair(..., lw_classes=["cWW"])
   ```

3. **Reuse visualizer**:
   ```python
   viz = PairVisualizer()  # Create once
   for pair in pairs:
       result = viz.visualize_pair(*pair)  # Reuse
   ```

## Examples

See `example_usage.py` for complete examples.
