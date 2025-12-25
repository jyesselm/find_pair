# Pair Visualization Module

Comprehensive visualization tools for base pair analysis, combining structural alignment, template comparison, and H-bond validation.

## Components

### 1. `pymol_generator.py`
Low-level PyMOL script generator.

**Classes:**
- `HBondViz`: H-bond visualization data
- `TemplateViz`: Template visualization data
- `PyMOLScriptGenerator`: Generates .pml scripts

**Features:**
- Configurable colors and styles
- Observed vs expected H-bond rendering
- Template overlay with RMSD labels
- Interactive toggle commands

### 2. `pair_visualizer.py`
High-level visualization API.

**Class:** `PairVisualizer`

**Key Methods:**
- `visualize_pair()`: Complete visualization pipeline
  - Loads PDB and extracts residues
  - Aligns templates to target
  - Loads H-bond data
  - Generates PyMOL script

**Returns:** `VisualizationResult` with all generated files and metrics

### 3. `cli.py`
Command-line interface.

**Usage:**
```bash
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72
```

**Options:**
- `--lw`: LW classes to show (default: cWW,tWW,cWS,tHS)
- `--pdb-dir`: PDB directory
- `--template-dir`: Template directory
- `--hbond-dir`: H-bond JSON directory
- `--output-dir`: Output directory
- `--no-hbonds`: Skip H-bond loading
- `--launch-pymol`: Launch PyMOL automatically

## Quick Start

### Example 1: Basic Usage
```python
from pathlib import Path
from pair_visualizer import PairVisualizer

visualizer = PairVisualizer(
    pdb_dir=Path("data/pdb"),
    template_dir=Path("basepair-idealized"),
    output_dir=Path("viz_output"),
)

result = visualizer.visualize_pair(
    pdb_id="1EHZ",
    res_id1="A-G-1",
    res_id2="A-C-72",
)

print(f"PyMOL script: {result.pymol_script}")
print(f"Best match: {result.best_lw} (RMSD={result.best_rmsd:.3f}Å)")
```

### Example 2: Custom LW Classes
```python
result = visualizer.visualize_pair(
    pdb_id="1EHZ",
    res_id1="A-G-5",
    res_id2="A-U-68",
    lw_classes=["cWW", "tWW", "tHS"],  # Only show these
    include_hbonds=True,
)
```

### Example 3: Command Line
```bash
# Basic visualization
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72

# Custom LW classes
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72 --lw cWW,tWW

# Launch PyMOL automatically
python cli.py --pdb 1EHZ --res1 A-G-1 --res2 A-C-72 --launch-pymol
```

## Pipeline Flow

```
1. Load PDB → extract pair residues
2. Load templates for each LW class
3. Align templates to target (ring atoms)
4. Write aligned template PDBs
5. Load H-bond data (observed + expected)
6. Generate PyMOL script with:
   - Target structure (cyan, thick)
   - Templates (colored, thin)
   - Observed H-bonds (yellow dashed)
   - Expected H-bonds (green/red dashed)
   - Atom labels
```

## Output Files

For `visualize_pair(pdb_id="1EHZ", res_id1="A-G-1", res_id2="A-C-72")`:

```
viz_output/
├── target_1EHZ_GC.pdb              # Target pair
├── aligned_1EHZ_cWW_GC.pdb         # Aligned cWW template
├── aligned_1EHZ_tWW_GC.pdb         # Aligned tWW template
├── aligned_1EHZ_cWS_GC.pdb         # Aligned cWS template
├── aligned_1EHZ_tHS_GC.pdb         # Aligned tHS template
└── view_1EHZ_A_G_1_A_C_72.pml      # PyMOL script
```

## PyMOL Script Features

### Visual Elements
- **Target**: Cyan thick sticks with atom labels
- **Templates**: Colored thin sticks (green=cWW, orange=tWW, etc.)
- **Observed H-bonds**: Yellow solid dashes (from slot_hbonds JSON)
- **Expected H-bonds**: Green (found) or red (missing) dotted dashes

### Interactive Commands
```pymol
# Toggle templates
disable cWW
enable cWW

# Toggle H-bonds
disable obs_*
enable exp_*

# Zoom to specific atoms
zoom name N1+N3+O6
```

### Legend (in script comments)
```
# Target = cyan (thick sticks)
# cWW = green (RMSD=0.123Å)
# tWW = orange (RMSD=1.456Å)
# cWS = magenta (RMSD=2.789Å)
#
# Yellow dashed = Observed H-bonds
# Green dashed = Expected H-bonds (found)
# Red dashed = Missing H-bonds
```

## H-Bond Visualization

### Observed H-bonds
Loaded from `data/json/slot_hbonds/<pdb_id>.json`:
- Base-base context only
- Yellow solid dashes
- Distance labels

### Expected H-bonds
Loaded from `hbond.patterns.get_cww_expected()`:
- Canonical Watson-Crick patterns (GC, AU, GU)
- Green if found in observed
- Red if missing
- Dotted dashes

## Alignment Details

### Ring Atoms Used
```python
RING_ATOMS = ["C2", "C4", "C5", "C6", "N1", "N3", "N7", "C8", "N9"]
```

### Alignment Process
1. Collect ring atoms from both residues
2. Kabsch alignment (minimize RMSD)
3. Apply transformation to ALL atoms (not just ring)
4. Compute RMSD on ring atoms
5. Select best match (minimum RMSD)

### Template Colors
```python
TEMPLATE_COLORS = {
    "cWW": "green",
    "tWW": "orange",
    "cWS": "magenta",
    "tWH": "yellow",
    "tHS": "pink",
    "cHS": "salmon",
}
```

## Directory Structure Requirements

```
data/
├── pdb/
│   └── 1EHZ.pdb                    # Input PDB files
└── json/
    └── slot_hbonds/
        └── 1EHZ.json               # H-bond data

basepair-idealized/
├── cWW/
│   ├── GC.pdb
│   ├── AU.pdb
│   └── ...
├── tWW/
│   └── ...
└── cWS/
    └── ...

viz_output/                          # Generated files
```

## Dependencies

- `numpy`: Kabsch alignment
- `core.pdb_parser`: PDB I/O
- `core.alignment`: Kabsch algorithm
- `core.residue`: Residue data structures
- `core.identifiers`: ID conversions
- `hbond.patterns`: Expected H-bond patterns

## Integration

This module is designed to work with:
- **Slot H-bond optimizer**: Uses `data/json/slot_hbonds/` output
- **Template generator**: Uses `basepair-idealized/` templates
- **LW classifier**: Compares RMSD to classify pairs

## Future Enhancements

1. **Multi-pair visualization**: Show all pairs in a structure
2. **Animation**: Morph between templates
3. **Density overlay**: Show electron density if available
4. **Custom H-bond criteria**: Override default patterns
5. **Batch mode**: Generate scripts for many pairs at once
6. **Web export**: Generate 3Dmol.js HTML files
