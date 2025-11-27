# File Cleanup Plan

## Summary

Many files were generated during the investigation and comparison process. Most are now obsolete since we've achieved 100% match rate.

## File Categories

### Essential Files (Keep in Root)
- `README.md` - Main project README
- `PROJECT_STATUS.md` - Current project status
- `DOCUMENTATION_INDEX.md` - Documentation index
- `FINAL_100_PDB_COMPARISON_RESULTS.md` - Final comparison results
- `3CME_EDGE_CASE_ANALYSIS.md` - Current edge case analysis
- `TIE_BREAKING_FIX_SUMMARY.md` - Tie-breaking fix documentation
- `TIE_BREAKING_FIX_FINAL_RESULTS.md` - Fix results
- `.cursorrules` - Cursor configuration
- `.gitignore` - Git ignore rules

### Obsolete Analysis Files (Archive to `docs/archive/`)

#### PDB-Specific Analysis (16 files) → `docs/archive/pdb_analysis/`
These document PDBs that had differences - now mostly fixed:
- `100D_mismatched_pairs_analysis.md`
- `1T0K_mismatched_pairs_analysis.md`
- `3DIZ_mismatched_pairs_analysis.md`
- `3G8T_mismatched_pairs_analysis.md`
- `3IWN_mismatched_pairs_analysis.md`
- `4JV5_mismatched_pairs_analysis.md`
- `5IWA_mismatched_pairs_analysis.md`
- `6CAP_mismatched_pairs_analysis.md`
- `6CAQ_mismatched_pairs_analysis.md`
- `6J6G_mismatched_pairs_analysis.md`
- `6ZXH_mismatched_pairs_analysis.md`
- `7XHT_mismatched_pairs_analysis.md`
- `7XUE_mismatched_pairs_analysis.md`
- `8T2T_mismatched_pairs_analysis.md`
- `deep_analysis_3CF5.md`
- `deep_analysis_top3.md`

#### Investigation Documents (8 files) → `docs/archive/investigation/`
- `REMAINING_0.5_PERCENT_ANALYSIS.md`
- `REMAINING_TASKS.md`
- `FINAL_COMPREHENSIVE_SUMMARY.md`
- `FINAL_INVESTIGATION_STATUS.md`
- `FINAL_STATUS_SUMMARY.md`
- `INVESTIGATION_SUMMARY_AND_NEXT_STEPS.md`
- `TIE_BREAKING_ANALYSIS.md`
- `find_bestpair_differences_analysis.md`
- `remaining_differences_analysis.md`
- `detailed_report.md`

#### Comparison Reports (7 files) → `docs/archive/comparison_reports/`
- `COMPLETE_100_PDB_COMPARISON_SUMMARY.md`
- `comparison_complete_100_pdbs.md`
- `comparison_final_100_pdbs.md`
- `comparison_report_100_pdbs.md`
- `comparison_report_100_pdbs_detailed.md`
- `comparison_results_summary.md`
- `comparison_summary_100_pdbs.md`
- `cleanup_summary.md`

### Temporary Data Files (Can be Removed)
- `bestpairs.pdb` - Temporary PDB output
- `bp_order.dat` - Temporary data file
- `hel_regions.pdb` - Temporary PDB output
- `ref_frames.dat` - Temporary data file
- `col_chains.scr` - Temporary script output
- `col_helices.scr` - Temporary script output
- `remaining_differences_analysis_detailed.json` - Obsolete analysis data
- `comparison_results_summary.json` - Can be regenerated

### Directories
- `base_frame_calc/` - Check if temporary or needed
- `data/` - Keep (contains PDB files and JSON)
- `docs/` - Keep (contains organized documentation)
- `scripts/` - Keep (contains Python scripts)

## Action Plan

1. **Archive obsolete files** to `docs/archive/` subdirectories
2. **Remove temporary data files** (PDB, DAT, SCR, obsolete JSON)
3. **Keep essential files** in root
4. **Update DOCUMENTATION_INDEX.md** with archive locations

## Commands to Execute

```bash
# Archive PDB-specific analysis
mv *_mismatched_pairs_analysis.md docs/archive/pdb_analysis/
mv deep_analysis_*.md docs/archive/pdb_analysis/

# Archive investigation documents
mv REMAINING_*.md docs/archive/investigation/
mv FINAL_*.md docs/archive/investigation/  # (except FINAL_100_PDB_COMPARISON_RESULTS.md)
mv INVESTIGATION_*.md docs/archive/investigation/
mv TIE_BREAKING_ANALYSIS.md docs/archive/investigation/
mv find_bestpair_differences_analysis.md docs/archive/investigation/
mv remaining_differences_analysis.md docs/archive/investigation/
mv detailed_report.md docs/archive/investigation/

# Archive comparison reports
mv COMPLETE_100_PDB_COMPARISON_SUMMARY.md docs/archive/comparison_reports/
mv comparison_*.md docs/archive/comparison_reports/
mv cleanup_summary.md docs/archive/comparison_reports/

# Remove temporary files
rm -f bestpairs.pdb bp_order.dat hel_regions.pdb ref_frames.dat
rm -f col_chains.scr col_helices.scr
rm -f remaining_differences_analysis_detailed.json
# Keep comparison_results_summary.json (might be useful)
```

## Result

After cleanup:
- **~9 essential files** in root
- **~31 files archived** for reference
- **~8 temporary files removed**
- Clean, organized root directory

