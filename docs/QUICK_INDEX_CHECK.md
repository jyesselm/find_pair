# Quick Index Validation Reference

**TL;DR**: ✅ All 2,386 PDBs with legacy data have perfect index matching. No issues found.

---

## Quick Commands

### Check summary

```bash
python3 scripts/summarize_validation.py
```

### Analyze any mismatches

```bash
python3 scripts/analyze_index_mismatches.py
```

### Validate all remaining PDBs (auto-resume)

```bash
# Tests only PDBs that haven't passed yet (auto-resumes)
python3 scripts/validate_all_indices.py --threads 8
```

### Validate ALL PDBs (full recheck)

```bash
# Revalidates everything, including those that already passed
python3 scripts/validate_all_indices.py --threads 8 --revalidate
```

### Keep generated files (don't clean up)

```bash
python3 scripts/validate_all_indices.py --no-clean
```

### Check a specific PDB

```bash
./build/generate_modern_json data/pdb/1EHZ.pdb data/json
```

---

## Current Status (Dec 2, 2025)

| Status | Count | Percentage |
|--------|-------|------------|
| ✅ PASS | 2,386 | 92.2% |
| ⏭️ SKIP | 199 | 7.7% |
| ⏱️ TIMEOUT | 2 | 0.1% |
| ❌ FAIL | 0 | 0.0% |

**Result**: 100% of PDBs with legacy data pass validation.

---

## What This Means

1. **Residue indices match perfectly** between legacy and modern code
2. **Comparisons are valid** - we're comparing the same residues
3. **No index-related bugs** - any differences are in calculations, not indexing
4. **Ready for production** - index matching is production-ready

---

## Files

- **Status CSV**: `data/index_validation_status.csv`
- **Full Report**: `docs/INDEX_VALIDATION_REPORT.md`
- **Validation Script**: `scripts/validate_all_indices.py`
- **Analysis Script**: `scripts/analyze_index_mismatches.py`
- **Summary Script**: `scripts/summarize_validation.py`

---

## If You Find a Failure

1. Run: `python3 scripts/analyze_index_mismatches.py --pdb {pdb_id}`
2. Check: `data/index_mapping/{pdb_id}.json`
3. Review: `data/index_mismatch_report.csv`
4. Investigate root cause
5. Fix and revalidate

---

## See Also

- `docs/INDEX_VALIDATION_REPORT.md` - Full detailed report
- `docs/COMPARISON_INDEX_RULES.md` - Rules for using indices in comparisons

