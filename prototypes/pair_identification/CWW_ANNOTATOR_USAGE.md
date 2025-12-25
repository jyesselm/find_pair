# cWW Miss Annotator - Usage Guide

## Overview

The cWW miss annotator analyzes Watson-Crick (cWW) base pair classification differences between our system and DSSR. It provides detailed diagnostics explaining why canonical cWW pairs are missed (false negatives) or why non-cWW pairs are incorrectly classified as cWW (false positives).

## Quick Start

```bash
# Analyze single PDB
python prototypes/pair_identification/cww_annotate.py --pdb 1EHZ

# Analyze 100 PDBs with 10 parallel workers
python prototypes/pair_identification/cww_annotate.py \
  --pdb-list data/fast_pdbs.json \
  --max-pdbs 100 \
  --workers 10
```

## Command-Line Options

```
Usage: cww_annotate.py [OPTIONS]

Input Options:
  --pdb PDB_ID              Single PDB ID to analyze
  --pdb-list PATH           JSON file with list of PDB IDs
  --max-pdbs N              Limit processing to first N PDBs

Directory Options:
  --pdb-dir PATH            PDB files directory (default: data/pdb)
  --hbond-dir PATH          Slot H-bond JSON directory (default: data/json/slot_hbonds)
  --dssr-dir PATH           DSSR JSON directory (default: data/json_dssr)
  --idealized-dir PATH      Idealized templates (default: basepair-idealized)
  --exemplar-dir PATH       Exemplar templates (default: basepair-exemplars)

Output Options:
  --output PATH, -o PATH    Output directory (default: results/cww_analysis)
  --verbose, -v             Print detailed per-PDB statistics

Processing Options:
  --workers N, -w N         Number of parallel workers (default: 10)
```

## Output Files

The tool generates two types of JSON files in the output directory:

### 1. Individual PDB Reports (`<PDB_ID>.json`)

Each file contains:
- Summary statistics (TP, FN, FP, sensitivity, precision)
- Detailed false negative annotations
- Detailed false positive annotations

Structure:
```json
{
  "pdb_id": "1EHZ",
  "summary": {
    "total_canonical_cww": 42,
    "true_positives": 38,
    "false_negatives": 4,
    "false_positives": 1,
    "sensitivity": 0.905,
    "precision": 0.974
  },
  "false_negatives": [...],
  "false_positives": [...]
}
```

Each annotation includes:
- Residue IDs and sequence
- Our prediction vs DSSR class
- Reason codes explaining the miss
- H-bond diagnostics (expected vs found)
- Geometric diagnostics (RMSD, angles, distances)

### 2. Aggregate Report (`aggregate.json`)

Combined statistics across all processed PDBs:
- Total counts (PDBs, cWW pairs, TP, FN, FP)
- Overall accuracy, sensitivity, precision
- Reason code distribution (global and per-sequence)
- Wrong atom error patterns
- RMSD gap distribution
- Per-sequence breakdown (GC, CG, AU, UA)

## Diagnostic Reason Codes

The annotator assigns reason codes to explain each miss:

| Code | Meaning |
|------|---------|
| `NO_HBONDS` | No H-bonds detected between bases |
| `MISSING_HBONDS` | Some expected canonical H-bonds missing |
| `WRONG_ATOMS` | H-bonds to wrong acceptor atoms |
| `EXTRA_HBONDS` | Unexpected H-bonds detected |
| `DISTANCE_ISSUES` | H-bond distances out of range (< 2.5 or > 3.5 Å) |
| `OVERLOADED_ACCEPTOR` | Acceptor has > 2 H-bonds |
| `RMSD_PREFERS_OTHER` | Better RMSD fit to non-cWW template |
| `GEOMETRIC_OUTLIER` | Angle > 15° or N1-N9 distance abnormal |
| `NON_CANONICAL` | DSSR Saenger classification is "--" |

## Examples

### Analyze Fast PDB Set

```bash
python prototypes/pair_identification/cww_annotate.py \
  --pdb-list data/fast_pdbs.json \
  --workers 10 \
  --output results/fast_set_cww
```

### Debug Specific PDB

```bash
python prototypes/pair_identification/cww_annotate.py \
  --pdb 1EHZ \
  --verbose \
  --output results/debug_1ehz
```

### Custom Data Directories

```bash
python prototypes/pair_identification/cww_annotate.py \
  --pdb-list my_pdbs.json \
  --pdb-dir /custom/path/to/pdbs \
  --hbond-dir /custom/path/to/hbonds \
  --dssr-dir /custom/path/to/dssr \
  --output /custom/output/path
```

### Analyze Small Test Set

```bash
python prototypes/pair_identification/cww_annotate.py \
  --pdb-list data/fast_pdbs.json \
  --max-pdbs 10 \
  --workers 4 \
  --verbose
```

## Expected Output

Console output shows progress and summary:

```
Processing 100 PDBs with 10 workers...

============================================================
Aggregate Results (100 PDBs)
============================================================
Total canonical cWW:     4,200
True positives:          3,950
False negatives:         250
False positives:         45
Accuracy:                94.05%
Sensitivity:             94.05%
Precision:               98.87%

Top miss reasons:
  MISSING_HBONDS                : 150 ( 60.0%)
  DISTANCE_ISSUES               :  80 ( 32.0%)
  WRONG_ATOMS                   :  45 ( 18.0%)
  GEOMETRIC_OUTLIER             :  30 ( 12.0%)
  RMSD_PREFERS_OTHER            :  25 ( 10.0%)

Per-sequence breakdown:
  GC: 1410/1500 correct (94.0%)
  CG: 1415/1500 correct (94.3%)
  AU: 555/600 correct (92.5%)
  UA: 570/600 correct (95.0%)

Results saved to results/cww_analysis/
  - Individual reports: 100 files
  - Aggregate report: aggregate.json
```

## Architecture

The tool consists of several modules:

```
cww_miss_annotator/
├── cli.py                    # Command-line interface
├── annotator.py              # Main orchestrator
├── loaders.py                # DSSR and H-bond JSON loaders
├── hbond_analyzer.py         # H-bond pattern analysis
├── geometric_analyzer.py     # RMSD and geometric quality
├── diagnostics.py            # Result dataclasses
└── report.py                 # Report generation and aggregation
```

### Key Components

1. **MissAnnotator** (`annotator.py`): Orchestrates the analysis workflow
   - Loads DSSR pairs and slot H-bonds
   - Identifies TP, FN, FP
   - Annotates each miss with diagnostics

2. **HBondAnalyzer** (`hbond_analyzer.py`): Analyzes H-bond patterns
   - Compares found vs expected canonical patterns
   - Identifies missing, extra, and wrong H-bonds
   - Checks distance ranges and acceptor saturation

3. **GeometricAnalyzer** (`geometric_analyzer.py`): Computes geometric metrics
   - RMSD to cWW and other LW templates
   - Interbase angle (plane-plane)
   - N1-N9 glycosidic distance
   - Geometric outlier detection

4. **AggregateReporter** (`report.py`): Aggregates statistics
   - Combines multiple PDB reports
   - Counts reason code distributions
   - Computes per-sequence statistics

## Dependencies

Required files and directories:
- PDB files in `data/pdb/`
- Slot H-bond JSON in `data/json/slot_hbonds/`
- DSSR JSON in `data/json_dssr/`
- Idealized templates in `basepair-idealized/`
- Exemplar templates in `basepair-exemplars/`
- `template_aligner.py` (in parent directory)

## Performance

- Uses multiprocessing for parallel PDB analysis
- Default: 10 workers
- Processing time: ~1-2 seconds per PDB (including RMSD calculations)
- Memory: ~100-200 MB per worker

## Limitations

1. Only analyzes standard Watson-Crick sequences (GC, CG, AU, UA)
2. Requires pre-computed slot H-bond data
3. RMSD calculations require template files
4. Does not handle modified bases beyond standard 4 nucleotides

## File Locations

- **Entry point**: `/Users/jyesselman2/local/code/cpp/find_pair_2/prototypes/pair_identification/cww_annotate.py`
- **Module directory**: `/Users/jyesselman2/local/code/cpp/find_pair_2/prototypes/pair_identification/cww_miss_annotator/`
- **Default output**: `/Users/jyesselman2/local/code/cpp/find_pair_2/results/cww_analysis/`
