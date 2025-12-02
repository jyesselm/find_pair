# Index Validation TODO

**Created**: December 2, 2025  
**Status**: Ready to validate remaining PDBs

---

## Current Status

✅ **2,386 PDBs (57.9%)** - Already passed validation  
⏳ **1,737 PDBs (42.1%)** - Need to be tested  
📊 **Total**: 4,123 PDBs in dataset

### Breakdown of Remaining

- **1,536 PDBs**: Not yet tested
- **199 PDBs**: Skipped (no legacy JSON)
- **2 PDBs**: Timeout (need retry)

---

## How to Run Full Validation

### Option 1: Quick Start (Recommended)

```bash
./RUN_FULL_VALIDATION.sh
```

This will:
- Auto-resume from where we left off
- Skip the 2,386 PDBs that already passed
- Test the remaining 1,737 PDBs
- Use 8 threads and batch size 100
- Clean up generated files automatically
- Stop on first failure

### Option 2: Custom Settings

```bash
# More threads for faster processing
python3 scripts/validate_all_indices.py --threads 12 --batch-size 200

# Keep generated files for inspection
python3 scripts/validate_all_indices.py --no-clean

# Revalidate everything (including 2,386 that already passed)
python3 scripts/validate_all_indices.py --revalidate
```

---

## Expected Outcomes

### Best Case (Most Likely)

- All 1,737 remaining PDBs validate successfully
- Final result: 100% of testable PDBs pass
- Generated files cleaned up automatically
- Total: ~2,385 passed, ~199 skipped (no legacy data)

### Worst Case

- Some PDBs fail validation (index mismatch found)
- Script stops on first failure
- Mapping file created at `data/index_mapping/{pdb_id}.json`
- Investigate with: `python3 scripts/analyze_index_mismatches.py`

---

## What Happens During Validation

For each PDB:

1. ✅ Read PDB file and track all residues
2. ✅ Apply modern code filtering
3. ✅ Assign modern indices (0-based)
4. ✅ Load legacy indices from JSON (1-based)
5. ✅ Match residues by chain/seq/insertion
6. ✅ Validate counts and indices match
7. ✅ **Clean up**: Delete generated JSON files if passed
8. ✅ Update status CSV
9. ✅ Continue to next PDB

On PASS:
- Record in CSV
- **Delete generated files** (saves disk space)
- Continue

On FAIL:
- Record in CSV
- **Keep mapping file** for debugging
- **STOP** for investigation

---

## Monitoring Progress

The script prints:
- ✅ for each PASS
- ❌ for each FAIL (stops here)
- ⏭️ for each SKIP
- Batch summaries

Progress is saved continuously to `data/index_validation_status.csv`

You can monitor in another terminal:
```bash
# Watch the CSV update
watch -n 5 "tail -20 data/index_validation_status.csv"

# Count status
awk -F',' 'NR>1 {print $6}' data/index_validation_status.csv | sort | uniq -c
```

---

## After Validation Completes

### Check Results

```bash
python3 scripts/summarize_validation.py
```

### Analyze Any Failures

```bash
python3 scripts/analyze_index_mismatches.py
```

### Update Documentation

If all pass:
- Update INDEX_VALIDATION_REPORT.md with final numbers
- Merge `fix-index-matching` branch to main

---

## Estimated Time

With 8 threads and 100 per batch:

- **Fast PDBs** (small): ~2 seconds each
- **Slow PDBs** (large): ~10 seconds each
- **Average**: ~5 seconds each

**Total estimate**: 1,737 PDBs × 5 sec / 8 threads = **18 minutes**

(May be faster due to batching and parallel processing)

---

## Disk Space

### Before Cleanup (--no-clean)

Each PDB generates ~13 JSON files:
- 1,737 PDBs × 13 files × ~50 KB = **~1.1 GB**

### After Cleanup (default)

Only failed PDBs keep files:
- Assuming 0 failures: **~0 MB**
- Assuming 10 failures: **~7 MB**

**Recommendation**: Use default cleanup (saves disk space)

---

## Troubleshooting

### If it stops with a failure

1. Check the error message
2. Look at mapping file: `data/index_mapping/{pdb_id}.json`
3. Run analysis: `python3 scripts/analyze_index_mismatches.py`
4. Investigate root cause
5. Fix the issue
6. Resume: `python3 scripts/validate_all_indices.py --start-from {pdb_id}`

### If you want to interrupt

- Press Ctrl+C
- Progress is saved in CSV
- Resume with: `python3 scripts/validate_all_indices.py` (auto-resume)

---

## When Complete

You'll have:

✅ Complete validation of all 4,123 PDBs  
✅ CSV with all results  
✅ Confidence that indices match (or knowledge of what doesn't)  
✅ Clean workspace (no unnecessary JSON files)  
✅ Ready to proceed with base pair comparisons

---

## Run It Now!

```bash
./RUN_FULL_VALIDATION.sh
```

Let it run, grab some coffee, and check back in ~20 minutes! ☕

