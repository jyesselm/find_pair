# Modified Nucleotide Residues

This document describes how modified nucleotides are handled in the codebase and documents known classification issues that affect validation.

## Overview

Modified nucleotides are non-standard bases found in RNA/DNA structures. They are identified by:
- **3-letter code**: e.g., `9DG`, `A23`, `5MU`, `PSU`
- **base_type**: Single letter indicating the parent base (lowercase = modified)
  - Uppercase (A, G, C, T, U, I, P): Standard bases
  - Lowercase (a, g, c, t, u, i, p): Modified bases

## Documenting Modified Residues

Use the tool `tools/document_modified_residues.py` to scan all PDBs and generate a comprehensive report:

```bash
# Generate modified residues report
python3 tools/document_modified_residues.py --workers 20

# Output: data/modified_residues.json
```

### Output Format

The generated `data/modified_residues.json` contains:

```json
{
  "summary": {
    "total_pdbs_scanned": 3602,
    "pdbs_with_modified_residues": 1510,
    "total_modified_residues": 5354,
    "unique_modified_types_count": 370
  },
  "modified_residue_types": {
    "9DG": {
      "pdb_count": 4,
      "official_name": "9-DEAZAGUANINE",
      "formula": "C6 H6 N4 O",
      "component_type": "non-polymer"
    },
    "A23": {
      "pdb_count": 39,
      "official_name": "ADENOSINE-5'-PHOSPHATE-2',3'-CYCLIC PHOSPHATE",
      "formula": "C10 H13 N5 O9 P2",
      "component_type": "RNA linking"
    }
  },
  "pdbs": {
    "1Q2R": {
      "count": 4,
      "unique_types": ["9DG"],
      "residues": [
        {
          "residue_name": "9DG",
          "base_type": "u",
          "chain_id": "A",
          "residue_seq": 387,
          "legacy_residue_idx": 1538
        }
      ]
    }
  }
}
```

## Statistics (as of December 2024)

| Metric | Value |
|--------|-------|
| Total PDBs scanned | 3,602 |
| PDBs with modified residues | 1,510 (42%) |
| Total modified residues | 5,354 |
| Unique modified types | 370 |

### Top 20 Modified Residue Types

| Code | Count | Official Name |
|------|-------|---------------|
| GTP | 278 | GUANOSINE-5'-TRIPHOSPHATE |
| CCC | 72 | CYTIDINE-5'-PHOSPHATE-2',3'-CYCLIC PHOSPHATE |
| ATP | 64 | ADENOSINE-5'-TRIPHOSPHATE |
| 5BU | 53 | 5-BROMO-URIDINE-5'-MONOPHOSPHATE |
| LCC | 51 | Locked cytidine derivative |
| OMC | 50 | O2'-METHYLYCYTIDINE-5'-MONOPHOSPHATE |
| SAM | 50 | S-ADENOSYLMETHIONINE |
| 5MU | 46 | 5-METHYLURIDINE 5'-MONOPHOSPHATE |
| CBV | 46 | 5-BROMOCYTIDINE 5'-(DIHYDROGEN PHOSPHATE) |
| LCG | 46 | Locked guanosine derivative |
| ADP | 42 | ADENOSINE-5'-DIPHOSPHATE |
| 5MC | 40 | 5-METHYLCYTIDINE-5'-MONOPHOSPHATE |
| A23 | 39 | ADENOSINE-5'-PHOSPHATE-2',3'-CYCLIC PHOSPHATE |
| GDP | 39 | GUANOSINE-5'-DIPHOSPHATE |
| SAH | 37 | S-ADENOSYL-L-HOMOCYSTEINE |
| H2U | 36 | 5,6-DIHYDROURIDINE-5'-MONOPHOSPHATE |
| AMP | 36 | ADENOSINE MONOPHOSPHATE |
| 1MA | 34 | 6-HYDRO-1-METHYLADENOSINE-5'-MONOPHOSPHATE |
| OMG | 31 | O2'-METHYLGUANOSINE-5'-MONOPHOSPHATE |
| EPE | 31 | HEPES buffer (not a nucleotide!) |

## Known Classification Issues

### Legacy vs Modern Classification Differences

The legacy code has an inconsistency in how it classifies modified nucleotides:

1. **Atom-based classification** (`residue_ident`): Determines purine/pyrimidine based on presence of atoms N7, C8, N9
2. **Template-based classification** (`base_ident`): Assigns base_type based on RMSD fit to templates

These can conflict for modified nucleotides that have unusual atom compositions.

### Problematic Residues

#### 9DG (9-Deazaguanine)

| Property | Value |
|----------|-------|
| Official Name | 9-DEAZAGUANINE |
| True Base Type | Guanine derivative (purine) |
| Assigned base_type | `u` (pyrimidine) |
| Issue | Has N7, C8 atoms but no N9; classified as purine by atoms but pyrimidine by template |

**Impact on dNN calculation:**
- Legacy uses RY classification → looks for N9 → finds C9 instead
- Modern uses base_type → uses N1 correctly for pyrimidine
- Result: ~3Å difference in dNN values

**Affected PDBs:** 1Q2R, 1Q2S, 7NQ4, 8OMR

#### A23 (Modified Adenosine)

| Property | Value |
|----------|-------|
| Official Name | ADENOSINE-5'-PHOSPHATE-2',3'-CYCLIC PHOSPHATE |
| True Base Type | Adenine derivative (purine) |
| Assigned base_type | `a` (modified adenine) |
| Issue | 9-atom RMSD > NT_CUTOFF → falls back to 6-atom fit → classified as pyrimidine |

**Impact on dNN calculation:**
- Legacy: 9-atom fit fails → 6-atom fit → RY=0 (pyrimidine) → uses N1
- Modern: Uses base_type 'a' → correctly uses N9
- Result: ~3Å difference in dNN values

**Affected PDBs:** 2XD0, 2XDD

#### EPE (HEPES Buffer)

| Property | Value |
|----------|-------|
| Official Name | 4-(2-HYDROXYETHYL)-1-PIPERAZINE ETHANESULFONIC ACID |
| True Base Type | Not a nucleotide! |
| Component Type | non-polymer |
| Issue | Gets processed as a modified base due to ring atoms |

**Affected PDBs:** 4E8R and others

## How Classification Works

### Step 1: Atom Detection (`residue_ident`)

```
RA_LIST = [C4, N3, C2, N1, C6, C5, N7, C8, N9]
                                   ^^^^^^^^^^^
                                   Purine-specific atoms

If N7, C8, or N9 found → kr > 0
If RMSD ≤ NT_CUTOFF (0.2618):
  - If kr > 0: Return PURINE (RY=1)
  - Else: Return PYRIMIDINE (RY=0)
If RMSD fails and kr > 0:
  - Try again with only 6 atoms (no N7/C8/N9)
  - If passes: Return PYRIMIDINE (RY=0)  ← BUG: purine atoms ignored!
```

### Step 2: Base Type Assignment (`base_ident`)

Based on residue name matching against base list, or by detecting key atoms:
- O6 present → guanine (g)
- N6 present → adenine (a)
- N4 present → cytosine (c)
- C5M present → thymine (t)
- Default → uridine (u)

### Step 3: N Atom Selection (`glyco_N`)

For dNN calculation, select the glycosidic nitrogen:
```
If RY == 1 (purine):  Find N9, fallback to any atom with '9'
If RY == 0 (pyrimidine): Find N1, fallback to any atom with '1'
```

**The bug:** RY classification and base_type can disagree for modified nucleotides.

## Impact on Validation

Stage 6 (pair_validation) failures due to dNN mismatches are **expected** for:
- PDBs containing 9DG residues
- PDBs containing A23 residues with high RMSD
- PDBs containing non-nucleotide compounds misclassified as bases

These are **not bugs in modern code** - modern code is more consistent by using base_type for all classification decisions.

## Regenerating Modified Residue Data

To update the modified residues documentation:

```bash
cd /path/to/find_pair_2

# Generate fresh data
python3 tools/document_modified_residues.py --workers 20

# View summary
python3 -c "
import json
with open('data/modified_residues.json') as f:
    data = json.load(f)
print(json.dumps(data['summary'], indent=2))
"
```

## References

- RCSB PDB Chemical Component Dictionary: https://www.rcsb.org/docs/general-help/ligand-structure-quality
- Modified nucleotide database: https://modomics.genesilico.pl/
- X3DNA modified base handling: Legacy code in `org/src/cmn_fncs.c`

