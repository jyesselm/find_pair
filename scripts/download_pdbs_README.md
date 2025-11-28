# PDB Download Script

Script to download PDB files from RCSB PDB and place them in `data/pdb/`.

## Quick Start

```bash
# Download PDBs from test set
python3 scripts/download_pdbs.py --test-set 100

# Download specific PDB IDs
python3 scripts/download_pdbs.py 1H4S 2BNA 3DNA

# Download from file
python3 scripts/download_pdbs.py --from-file pdb_list.txt
```

## Features

- ✅ Downloads from RCSB PDB (https://files.rcsb.org/)
- ✅ Validates downloaded files
- ✅ Skips existing files (resume capability)
- ✅ Handles errors gracefully
- ✅ Progress tracking
- ✅ Download log saved

## Usage

### Download from Test Set

```bash
python3 scripts/download_pdbs.py --test-set 100
```

This downloads all PDBs in the test set of size 100.

### Download Specific PDB IDs

```bash
python3 scripts/download_pdbs.py 1H4S 2BNA 3DNA 4ABC
```

### Download from File

Create a file `pdb_list.txt` with one PDB ID per line:
```
1H4S
2BNA
3DNA
4ABC
```

Then download:
```bash
python3 scripts/download_pdbs.py --from-file pdb_list.txt
```

### Skip Existing Files

Resume interrupted downloads by skipping files that already exist:

```bash
python3 scripts/download_pdbs.py --test-set 1000 --skip-existing
```

### Continue on Error

Continue downloading even if some PDBs fail:

```bash
python3 scripts/download_pdbs.py --test-set 500 --continue-on-error
```

### Use Threading for Faster Downloads

Download multiple PDBs in parallel:

```bash
python3 scripts/download_pdbs.py --test-set 100 --threads 10
```

Default is 5 threads. Use more threads for faster downloads, but be respectful of RCSB servers.

### Adjust Download Delay

Be nicer to RCSB servers by increasing delay between downloads:

```bash
python3 scripts/download_pdbs.py --test-set 100 --delay 1.0
```

Note: Delay applies per thread in multi-threaded mode.

## Options

- `--test-set SIZE` - Download PDBs from test set (10, 50, 100, 500, 1000)
- `--from-file FILE` - Load PDB IDs from file (one per line)
- `--output-dir DIR` - Output directory (default: data/pdb/)
- `--skip-existing` - Skip files that already exist
- `--continue-on-error` - Continue even if some downloads fail
- `--delay SECONDS` - Delay between downloads (default: 0.5)
- `--retry N` - Number of retry attempts (default: 3)
- `--validate` - Validate downloaded files (default: True)
- `--threads N` - Number of parallel download threads (default: 5, use 0 for single-threaded)

## Output

Files are saved to `data/pdb/{PDB_ID}.pdb`.

A download log is saved to `data/pdb/.download_log.json` with:
- Successful downloads
- Skipped files
- Failed downloads
- Download statistics

## Examples

### Download test set of 1000 PDBs

```bash
python3 scripts/download_pdbs.py --test-set 1000 --skip-existing
```

### Download specific PDBs with validation

```bash
python3 scripts/download_pdbs.py 1H4S 2BNA --validate
```

### Resume failed downloads

```bash
# Check log
cat data/pdb/.download_log.json

# Retry failed PDBs
python3 scripts/download_pdbs.py $(jq -r '.failed[]' data/pdb/.download_log.json) --continue-on-error
```

## Integration with Cluster Scripts

Before running cluster comparisons, ensure PDBs are downloaded:

```bash
# 1. Download PDBs
python3 scripts/download_pdbs.py --test-set 1000

# 2. Generate test sets (if needed)
python3 scripts/compare_json.py generate-test-sets

# 3. Run cluster comparison
python3 scripts/cluster/run_cluster_comparison.py --test-set 1000
```

## Troubleshooting

### Download fails

- Check internet connection
- Try increasing delay: `--delay 1.0`
- Check if PDB ID is valid (visit https://www.rcsb.org/structure/{PDB_ID})
- Use `--continue-on-error` to continue with other PDBs

### Invalid files

- Files are automatically validated
- Invalid files are deleted
- Re-download with: `python3 scripts/download_pdbs.py {PDB_ID}`

### Rate limiting

If you see 429 errors (rate limited):
- Increase delay: `--delay 2.0`
- Download in smaller batches
- Wait and retry later

## Related Scripts

- `scripts/cluster/run_cluster_comparison.py` - Run comparisons on cluster
- `scripts/compare_json.py` - Compare JSON files

