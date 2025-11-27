# Legacy JSON Generation - Batch Processing

**Status**: ✅ **IN PROGRESS** - Generating legacy JSON for all PDBs

---

## Overview

We're batch-generating legacy JSON files for all PDB files in the repository to expand our test coverage beyond the original 100 PDBs.

---

## Progress

### Current Status
- **Total PDB files**: 4,934
- **PDBs ≤ 10 MB**: 4,923 (target for generation)
- **PDBs with legacy JSON**: 833+ (and growing)
- **Progress**: ~17% complete
- **Remaining**: ~4,090 PDBs

### Size Distribution (Current)
- **< 1 MB**: 792 PDBs with legacy JSON
- **1-5 MB**: 31 PDBs with legacy JSON  
- **5-10 MB**: 10 PDBs with legacy JSON

---

## Tools Created

### 1. `scripts/generate_legacy_json_batch.py`
**Purpose**: Batch generate legacy JSON for multiple PDBs

**Features**:
- Processes PDBs sorted by size (smallest first) for efficiency
- Parallel processing with configurable worker count
- Skips PDBs that already have legacy JSON
- Progress tracking and time estimates
- Handles timeouts and errors gracefully

**Usage**:
```bash
# Generate for all PDBs <= 10 MB (8 parallel workers)
python3 scripts/generate_legacy_json_batch.py --max-size 10 --workers 8

# Generate for smaller PDBs only (faster)
python3 scripts/generate_legacy_json_batch.py --max-size 1 --workers 8

# Limit number of PDBs to process
python3 scripts/generate_legacy_json_batch.py --max-size 5 --limit 100 --workers 4

# Resume from a specific PDB ID
python3 scripts/generate_legacy_json_batch.py --max-size 10 --start-from 1ABC --workers 8
```

**Options**:
- `--max-size MB`: Maximum file size in MB (default: 20)
- `--limit N`: Limit number of PDBs to process
- `--workers N`: Number of parallel workers (default: 4)
- `--start-from PDB_ID`: Start from this PDB ID (alphabetically)
- `--skip-existing`: Skip PDBs that already have legacy JSON (default: True)

### 2. `scripts/check_legacy_json_progress.py`
**Purpose**: Check current progress of legacy JSON generation

**Usage**:
```bash
python3 scripts/check_legacy_json_progress.py
```

**Output**:
- Total PDB files
- PDBs with legacy JSON
- Progress percentage
- Size distribution
- Remaining count

---

## Generation Process

### How It Works

1. **Find Legacy Executable**: Locates `org/build/bin/find_pair_analyze`
2. **Filter PDBs**: Selects PDBs based on size and existing JSON
3. **Sort by Size**: Processes smallest files first for efficiency
4. **Parallel Processing**: Uses ProcessPoolExecutor for concurrent generation
5. **Validation**: Checks for created JSON files after execution
6. **Progress Tracking**: Reports progress and time estimates

### Output Structure

Legacy JSON is written to `data/json_legacy/` in a segmented structure:
```
data/json_legacy/
├── find_bestpair_selection/
│   ├── 1ABC.json
│   ├── 1DEF.json
│   └── ...
├── base_pair/
│   ├── 1ABC.json
│   └── ...
├── base_frame_calc/
│   ├── 1ABC.json
│   └── ...
└── ...
```

### Key Indicator

The presence of `find_bestpair_selection/{PDB_ID}.json` indicates that legacy JSON has been successfully generated for that PDB.

---

## Performance

### Processing Speed
- **Small PDBs (< 1 MB)**: ~0.1-0.5 seconds per PDB
- **Medium PDBs (1-5 MB)**: ~1-5 seconds per PDB
- **Large PDBs (5-10 MB)**: ~5-30 seconds per PDB

### With 8 Workers
- **Estimated time for 4,923 PDBs**: ~2-4 hours (depending on size distribution)
- **Current rate**: ~100-200 PDBs per minute (for small files)

---

## Monitoring

### Check Progress
```bash
# Quick progress check
python3 scripts/check_legacy_json_progress.py

# Count files directly
ls data/json_legacy/find_bestpair_selection/*.json | wc -l
```

### View Logs
If running with `tee`:
```bash
tail -f /tmp/legacy_json_generation.log
```

### Check Running Process
```bash
ps aux | grep generate_legacy_json_batch
```

---

## Error Handling

### Common Issues

1. **Timeout**: Large PDBs may timeout (default: 10 minutes)
   - **Solution**: Process large files separately or increase timeout

2. **Path Issues**: Legacy executable must be run from `org/` directory
   - **Solution**: Script handles this automatically

3. **Missing Executable**: Legacy executable not built
   - **Solution**: Build with `cd org && mkdir -p build && cd build && cmake .. && make`

4. **File System Sync**: Files may take a moment to appear
   - **Solution**: Script includes a small delay before checking

---

## Next Steps

1. ✅ **Continue Generation**: Let the batch script run to completion
2. ✅ **Monitor Progress**: Use `check_legacy_json_progress.py` periodically
3. ✅ **Test Generated JSON**: Once complete, test modern code against all generated legacy JSON
4. ✅ **Expand Coverage**: Consider generating for PDBs > 10 MB if needed

---

## Benefits

### Expanded Test Coverage
- **Before**: 100 PDBs with legacy JSON
- **After**: 4,900+ PDBs with legacy JSON (49x increase)

### Better Validation
- More diverse structures tested
- Edge cases more likely to be discovered
- Confidence in 100% match rate across broader dataset

### Future Testing
- Can test modern code against thousands of PDBs
- Identify any regressions quickly
- Validate fixes across diverse structures

---

## Related Documentation

- `docs/JSON_GENERATION_SYSTEM.md`: Modern JSON generation system
- `docs/FINAL_TESTING_SUMMARY.md`: Current 100% match status
- `scripts/rebuild_json.py`: Alternative JSON regeneration tool

