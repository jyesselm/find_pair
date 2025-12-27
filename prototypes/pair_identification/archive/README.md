# Archived Analysis Scripts

This directory contains one-off analysis scripts that were used during development and debugging of the pair identification prototype. These scripts are no longer actively maintained but are preserved for reference.

## Script Categories

### DSSR Comparison Scripts
- `compare_cww.py` - Compare CWW pair detection against DSSR
- `compare_frames_dssr.py` - Compare reference frame calculations with DSSR
- `compare_validation_dssr.py` - Compare validation metrics with DSSR

### Failure Analysis Scripts
- `analyze_cww_misses.py` - Analyze cases where CWW pairs were missed
- `analyze_greedy_conflicts.py` - Analyze conflicts in greedy pair selection
- `analyze_standard_wc.py` - Analyze Watson-Crick pair detection
- `analyze_validation_failures.py` - Analyze validation failures across PDBs

### Validation Scripts
- `validate_classifier.py` - Validate template classification logic
- `validate_combined_classifier.py` - Validate combined classifier approach
- `validate_cww_parallel.py` - Parallel validation of CWW pairs

### Pair Finding Scripts
- `find_cww_pairs.py` - Find and analyze CWW pairs in structures
- `run_all_pdbs_sequential.py` - Run pair finding sequentially on test set

### Visualization Scripts
- `generate_alignment_viz.py` - Generate alignment visualizations
- `visualize_outlier.py` - Visualize outlier pairs
- `show_rmsd_comparison.py` - Show RMSD comparison plots

### Metric Calculation
- `calculate_adjusted_metrics.py` - Calculate adjusted validation metrics

### Template Generation
- `template_overlay.py` - Overlay templates for visualization
- `template_generator.py` - Generate reference templates
- `generate_slot_hbonds.py` - Generate H-bond slot assignments

## Usage

These scripts are not guaranteed to work with the current codebase and may have hardcoded paths or dependencies. They are kept for reference only.

For current analysis tools, use:
- Main prototype modules in `core/` and `hbond/`
- CLI tools in `cli/`
- Current test suite in `tests/`

## Migration Notes

If you need to resurrect any of these scripts:
1. Check dependencies and imports
2. Update paths to use current config system
3. Consider integrating functionality into main modules rather than standalone scripts
