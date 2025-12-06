#!/usr/bin/env python3
"""
Unified CLI for X3DNA JSON validation and comparison.

Single entry point that replaces all the individual validate_*.py and test_*.py scripts.

Usage:
    fp2-validate                         # All stages, all fast PDBs
    fp2-validate frames                  # Frames only
    fp2-validate hbonds atoms            # Multiple stages
    fp2-validate compare 1EHZ --verbose  # Detailed comparison
    fp2-validate --stop-on-first         # Stop at first failure
    fp2-validate --diff                  # Document differences
    fp2-validate --quiet                 # Exit code only (for CI)
"""

import sys
from pathlib import Path
import multiprocessing

import click

from .runner import ValidationRunner
from .pdb_list import get_pdb_list, load_test_set
from .json_comparison import JsonComparator
from .config import load_config


@click.group()
@click.version_option(version='1.0.0')
def main():
    """X3DNA JSON validation and comparison tool.
    
    Validates that modern JSON output matches legacy output across all stages.
    """
    pass


@main.command()
@click.argument('stages', nargs=-1)
@click.option('--pdb', '-p', multiple=True, help='Specific PDB(s) to validate')
@click.option('--max', '-n', 'max_count', type=int, help='Maximum number of PDBs to process')
@click.option('--test-set', type=click.Choice(['10', '50', '100', '500', '1000']),
              help='Use a saved test set of specified size')
@click.option('--workers', '-w', type=int, help='Number of parallel workers (default: 10)')
@click.option('--quiet', '-q', is_flag=True, help='Suppress output, exit code only')
@click.option('--verbose', '-v', is_flag=True, help='Show per-PDB results')
@click.option('--stop-on-first', '-s', is_flag=True, help='Stop at first failure for debugging')
@click.option('--diff', is_flag=True, help='Document all differences to file')
@click.option('--diff-file', type=click.Path(path_type=Path), help='Custom differences output file')
@click.option('--checkpoint', type=click.Path(path_type=Path), 
              help='Save progress to checkpoint file for resume')
@click.option('--resume', is_flag=True, 
              help='Resume from checkpoint, skipping already-passed PDBs')
@click.option('--clean-on-match', is_flag=True,
              help='Delete modern JSON files that match legacy (save disk space)')
@click.option('--project-root', type=click.Path(path_type=Path, exists=True),
              help='Project root directory (default: current directory)')
def validate(stages, pdb, max_count, test_set, workers, quiet, verbose, 
             stop_on_first, diff, diff_file, checkpoint, resume, clean_on_match,
             project_root):
    """Validate legacy vs modern JSON outputs.
    
    \b
    STAGES (12 stages in legacy execution order):
      1  pdb_atoms            - Atom parsing
      2  residue_indices      - Residue to atom mapping  
      3  base_frame_calc      - Base frame calculation
      4  ls_fitting           - Least squares fitting
      5  frame_calc           - Reference frame calculation
      6  pair_validation      - Pair validation checks
      7  distance_checks      - Distance measurements
      8  hbond_list           - Hydrogen bond list
      9  base_pair            - Base pair records
      10 find_bestpair_selection - Final pair selection
      11 bpstep_params        - Step parameters
      12 helical_params       - Helical parameters
    
    \b
    STAGE GROUPS (for convenience):
      atoms     = 1
      frames    = 3,4,5
      pairs     = 6,7,9,10
      hbonds    = 8
      steps     = 11,12
      all       = all stages
    
    \b
    Examples:
        fp2-validate 1                    # Stage 1 (atoms) only
        fp2-validate 3 4 5                # Stages 3-5 (frames)
        fp2-validate frames               # Same as above (group)
        fp2-validate pdb_atoms            # By stage name
        fp2-validate --pdb 1EHZ -v        # Single PDB, verbose
        fp2-validate --stop-on-first      # Debug mode
        fp2-validate --diff               # Document differences
        fp2-validate --test-set 100       # Use saved test set
    
    \b
    Checkpoint/Resume:
        fp2-validate --checkpoint run.json     # Save progress
        fp2-validate --checkpoint run.json --resume  # Resume from checkpoint
    
    \b
    Clean up matched files:
        fp2-validate --clean-on-match     # Delete modern JSON that matches
    """
    # Determine project root
    if project_root is None:
        # Try to find project root by looking for data/valid_pdbs_fast.json
        cwd = Path.cwd()
        if (cwd / 'data' / 'valid_pdbs_fast.json').exists():
            project_root = cwd
        elif (cwd.parent / 'data' / 'valid_pdbs_fast.json').exists():
            project_root = cwd.parent
        else:
            project_root = cwd
    
    # Convert stages to list
    stages_list = list(stages) if stages else None
    
    # Create runner
    runner = ValidationRunner(
        project_root=project_root,
        stages=stages_list,
        workers=workers,
        quiet=quiet,
        verbose=verbose,
        stop_on_first=stop_on_first,
        document_differences=diff,
        diff_file=diff_file,
        checkpoint_file=checkpoint,
        resume=resume,
        clean_on_match=clean_on_match,
    )
    
    # Get PDB list
    try:
        pdb_ids = runner.get_pdb_list(
            specific=list(pdb) if pdb else None,
            max_count=max_count,
            test_set=int(test_set) if test_set else None
        )
    except FileNotFoundError as e:
        if not quiet:
            click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    
    if not pdb_ids:
        if not quiet:
            click.echo("No PDBs to validate!", err=True)
        sys.exit(1)
    
    # Run validation
    summary = runner.validate(pdb_ids)
    
    # Exit code: 0 if all pass, 1 if any fail
    sys.exit(0 if summary.all_passed else 1)


@main.command('list-pdbs')
@click.option('--test-set', type=click.Choice(['10', '50', '100', '500', '1000']),
              help='List PDBs in a test set')
@click.option('--max', '-n', 'max_count', type=int, help='Maximum number to show')
@click.option('--project-root', type=click.Path(path_type=Path, exists=True))
def list_pdbs(test_set, max_count, project_root):
    """List available PDBs for validation."""
    if project_root is None:
        project_root = Path.cwd()
    
    try:
        pdb_ids = get_pdb_list(
            project_root,
            test_set=int(test_set) if test_set else None,
            max_count=max_count
        )
    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    
    click.echo(f"Total PDBs: {len(pdb_ids)}")
    if test_set:
        click.echo(f"Test set: {test_set}")
    
    for pdb_id in pdb_ids:
        click.echo(f"  {pdb_id}")


@main.command('info')
@click.option('--project-root', type=click.Path(path_type=Path, exists=True))
def info(project_root):
    """Show validation environment info."""
    if project_root is None:
        project_root = Path.cwd()
    
    click.echo("X3DNA Validation Environment")
    click.echo("=" * 40)
    click.echo(f"Project root: {project_root}")
    
    # Check for executables
    legacy_exe = project_root / 'org' / 'build' / 'bin' / 'find_pair_analyze'
    modern_exe = project_root / 'build' / 'generate_modern_json'
    
    click.echo(f"Legacy executable: {'✅' if legacy_exe.exists() else '❌'} {legacy_exe}")
    click.echo(f"Modern executable: {'✅' if modern_exe.exists() else '❌'} {modern_exe}")
    
    # Check for data files
    fast_pdbs = project_root / 'data' / 'valid_pdbs_fast.json'
    if fast_pdbs.exists():
        try:
            pdb_ids = get_pdb_list(project_root)
            click.echo(f"Fast PDBs: {len(pdb_ids)} available")
        except Exception:
            click.echo("Fast PDBs: ❌ Could not load")
    else:
        click.echo(f"Fast PDBs: ❌ {fast_pdbs} not found")
    
    # Check for test sets
    click.echo("\nTest sets:")
    for size in [10, 50, 100, 500, 1000]:
        test_pdbs = load_test_set(project_root, size)
        if test_pdbs:
            click.echo(f"  {size}: ✅ {len(test_pdbs)} PDBs")
        else:
            click.echo(f"  {size}: ❌ Not found")


# Common options for alias commands
def alias_options(f):
    """Decorator to add common options to alias commands."""
    f = click.option('--pdb', '-p', multiple=True, help='Specific PDB(s) to validate')(f)
    f = click.option('--max', '-n', 'max_count', type=int, help='Maximum number of PDBs')(f)
    f = click.option('--test-set', type=click.Choice(['10', '50', '100', '500', '1000']),
                     help='Use a saved test set')(f)
    f = click.option('--workers', '-w', type=int, help='Number of parallel workers')(f)
    f = click.option('--quiet', '-q', is_flag=True, help='Suppress output')(f)
    f = click.option('--verbose', '-v', is_flag=True, help='Show per-PDB results')(f)
    f = click.option('--stop-on-first', '-s', is_flag=True, help='Stop at first failure')(f)
    f = click.option('--diff', is_flag=True, help='Document differences to file')(f)
    f = click.option('--checkpoint', type=click.Path(path_type=Path), 
                     help='Save progress to checkpoint file')(f)
    f = click.option('--resume', is_flag=True, help='Resume from checkpoint')(f)
    f = click.option('--clean-on-match', is_flag=True, 
                     help='Delete modern JSON that matches legacy')(f)
    return f


@main.command('atoms')
@alias_options
@click.pass_context
def atoms(ctx, **kwargs):
    """Validate atom records (alias for: validate atoms)."""
    ctx.invoke(validate, stages=['atoms'], **kwargs)


@main.command('frames')
@alias_options
@click.pass_context
def frames(ctx, **kwargs):
    """Validate frame calculations (alias for: validate frames)."""
    ctx.invoke(validate, stages=['frames'], **kwargs)


@main.command('hbonds')
@alias_options
@click.pass_context
def hbonds(ctx, **kwargs):
    """Validate H-bond lists (alias for: validate hbonds)."""
    ctx.invoke(validate, stages=['hbonds'], **kwargs)


@main.command('pairs')
@alias_options
@click.pass_context
def pairs(ctx, **kwargs):
    """Validate base pairs (alias for: validate pairs)."""
    ctx.invoke(validate, stages=['pairs'], **kwargs)


@main.command('steps')
@alias_options
@click.pass_context
def steps(ctx, **kwargs):
    """Validate step parameters (alias for: validate steps)."""
    ctx.invoke(validate, stages=['steps'], **kwargs)


@main.command('compare')
@click.argument('pdb_ids', nargs=-1)
@click.option('--verbose', '-v', is_flag=True, help='Detailed field-by-field comparison (single PDB only)')
@click.option('--output', '-o', type=click.Path(path_type=Path), help='Save report to file')
@click.option('--test-set', type=click.Choice(['10', '50', '100', '500', '1000']),
              help='Use a saved test set')
@click.option('--workers', '-w', type=int, help='Number of parallel workers')
@click.option('--project-root', type=click.Path(path_type=Path, exists=True),
              help='Project root (default: current directory)')
def compare(pdb_ids, verbose, output, test_set, workers, project_root):
    """Generate detailed comparison reports with verbose field-by-field output.
    
    \b
    Examples:
        fp2-validate compare 1EHZ --verbose          # Detailed comparison
        fp2-validate compare --test-set 10           # Compare 10 PDBs
        fp2-validate compare 1EHZ -o report.txt      # Save to file
    """
    from .verbose_reporter import VerboseReporter, create_record_comparison_from_dicts
    
    # Determine project root
    if not project_root:
        project_root = Path.cwd()
    
    # Get PDB IDs
    if test_set:
        pdb_ids = load_test_set(project_root, int(test_set))
        if not pdb_ids:
            click.echo(f"Error: Test set {test_set} not found.", err=True)
            sys.exit(1)
    elif not pdb_ids:
        click.echo("Error: Specify PDB ID(s) or use --test-set", err=True)
        sys.exit(1)
    
    # Load config and create comparator
    config = load_config(project_root / "comparison_config.yaml")
    comparator = JsonComparator(
        tolerance=config.get("tolerance", 1e-6),
        compare_atoms=True, compare_frames=True, compare_steps=True,
        compare_pairs=True, compare_hbond_list=True, compare_residue_indices=True
    )
    
    # Run comparisons
    click.echo(f"Comparing {len(pdb_ids)} PDB(s)...")
    results = {}
    for pdb_id in pdb_ids:
        # Note: These files may not exist - JsonComparator handles split files
        legacy_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
        modern_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
        pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
            
        result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)
        
        # Check for errors
        if result.status == 'error':
            click.echo(f"Error comparing {pdb_id}: {', '.join(result.errors)}", err=True)
            continue
            
        results[pdb_id] = result
    
    if not results:
        click.echo("Error: No valid comparisons completed", err=True)
        sys.exit(1)
    
    # Generate report
    if verbose and len(pdb_ids) == 1:
        # Full verbose mode with field-by-field comparison
        from .verbose_reporter import VerboseReporter, create_record_comparison_from_dicts
        
        pdb_id = list(pdb_ids)[0]
        result = results.get(pdb_id)
        if not result:
            click.echo(f"Error: No result for {pdb_id}", err=True)
            sys.exit(1)
        
        # Use VerboseReporter for full field-by-field output
        reporter = VerboseReporter(
            tolerance=comparator.tolerance,
            show_provenance=False,
            show_related=True,
            diff_only=False
        )
        
        # Build stages list
        stages = []
        if result.atom_comparison:
            stages.append("atoms")
        if result.frame_comparison:
            stages.append("frames")
        if result.step_comparison:
            stages.append("steps")
        if result.helical_comparison:
            stages.append("helical")
        if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
            stages.append("distance_checks")
        if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
            stages.append("hbond_list")
        if hasattr(result, 'pair_validation_comparison') and result.pair_validation_comparison:
            stages.append("pair_validation")
        if hasattr(result, 'find_bestpair_comparison') and result.find_bestpair_comparison:
            stages.append("find_bestpair_selection")
        if hasattr(result, 'base_pair_comparison') and result.base_pair_comparison:
            stages.append("base_pair")
        
        reporter.add_header(pdb_id, stages)
        
        # Add distance checks
        if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
            dcc = result.distance_checks_comparison
            reporter.add_stage_header("distance_checks", 3)
            reporter.add_stage_summary(
                dcc.total_legacy, dcc.total_modern,
                dcc.total_legacy - len(dcc.missing_in_modern),
                len(dcc.missing_in_modern), len(dcc.extra_in_modern),
                len(dcc.mismatched_checks)
            )
            for mismatch in dcc.mismatched_checks[:reporter.max_mismatches_per_stage]:
                pair_key = (mismatch.get('base_i'), mismatch.get('base_j'))
                rec_comp = create_record_comparison_from_dicts(
                    record_key=pair_key,
                    record_type="distance_checks",
                    legacy_record=mismatch.get('legacy_record', {}),
                    modern_record=mismatch.get('modern_record', {}),
                    fields_to_compare=['dorg', 'dNN', 'plane_angle', 'd_v', 'overlap_area'],
                    tolerance=comparator.tolerance,
                    legacy_source=f"data/json_legacy/distance_checks/{pdb_id}.json",
                    modern_source=f"data/json/distance_checks/{pdb_id}.json"
                )
                reporter.add_record_comparison(rec_comp)
        
        # Add hbond list
        if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
            hlc = result.hbond_list_comparison
            reporter.add_stage_header("hbond_list", 4)
            reporter.add_stage_summary(
                hlc.total_legacy, hlc.total_modern,
                hlc.total_legacy - len(hlc.missing_in_modern),
                len(hlc.missing_in_modern), len(hlc.extra_in_modern),
                len(hlc.mismatched_pairs)
            )
            for mismatch in hlc.mismatched_pairs[:reporter.max_mismatches_per_stage]:
                pair_key = (mismatch.get('base_i'), mismatch.get('base_j'))
                rec_comp = create_record_comparison_from_dicts(
                    record_key=pair_key,
                    record_type="hbond_list",
                    legacy_record=mismatch.get('legacy_record', {}),
                    modern_record=mismatch.get('modern_record', {}),
                    fields_to_compare=['num_hbonds'],
                    tolerance=comparator.tolerance,
                    legacy_source=f"data/json_legacy/hbond_list/{pdb_id}.json",
                    modern_source=f"data/json/hbond_list/{pdb_id}.json"
                )
                reporter.add_record_comparison(rec_comp)
        
        # Add frames
        if result.frame_comparison:
            fc = result.frame_comparison
            reporter.add_stage_header("frames", 2)
            reporter.add_stage_summary(
                fc.total_legacy, fc.total_modern,
                fc.total_legacy - len(fc.missing_residues),
                len(fc.missing_residues),
                0,
                len(fc.mismatched_calculations)
            )
            for mismatch in fc.mismatched_calculations[:reporter.max_mismatches_per_stage]:
                residue_key = mismatch.residue_key
                record_type = "unknown"
                mismatches_dict = mismatch.mismatches
                if 'base_frame_calc' in mismatches_dict:
                    record_type = "base_frame_calc"
                    fields = ['rms_fit', 'num_matched_atoms', 'base_type']
                elif 'ls_fitting' in mismatches_dict:
                    record_type = "ls_fitting"
                    fields = ['rms_fit', 'num_points']
                elif 'frame_calc' in mismatches_dict:
                    record_type = "frame_calc"
                    fields = ['rms_fit', 'num_matched_atoms']
                else:
                    fields = ['rms_fit']
                
                rec_comp = create_record_comparison_from_dicts(
                    record_key=residue_key,
                    record_type=record_type,
                    legacy_record=mismatch.legacy_record,
                    modern_record=mismatch.modern_record,
                    fields_to_compare=fields,
                    tolerance=comparator.tolerance,
                    legacy_source=f"data/json_legacy/{record_type}/{pdb_id}.json",
                    modern_source=f"data/json/{record_type}/{pdb_id}.json"
                )
                reporter.add_record_comparison(rec_comp)
        
        # Add summary
        stages_compared = len(stages)
        stages_with_diffs = 0
        diff_details = []
        
        if result.frame_comparison and len(result.frame_comparison.mismatched_calculations) > 0:
            stages_with_diffs += 1
            diff_details.append(f"frames: {len(result.frame_comparison.mismatched_calculations)} mismatches")
        
        if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
            if len(result.distance_checks_comparison.mismatched_checks) > 0:
                stages_with_diffs += 1
                diff_details.append(f"distance_checks: {len(result.distance_checks_comparison.mismatched_checks)} mismatches")
        
        if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
            if len(result.hbond_list_comparison.mismatched_pairs) > 0:
                stages_with_diffs += 1
                diff_details.append(f"hbond_list: {len(result.hbond_list_comparison.mismatched_pairs)} mismatches")
        
        reporter.add_summary(
            stages_compared=stages_compared,
            perfect_matches=stages_compared - stages_with_diffs,
            stages_with_diffs=stages_with_diffs,
            diff_details=diff_details
        )
        
        report = reporter.generate_report()
    else:
        # Summary mode
        if verbose: click.echo("Warning: Verbose only works for single PDB", err=True)
        total = len(results)
        matches = sum(1 for r in results.values() if not r.has_differences())
        lines = []
        lines.append(f"Comparison Summary: {total} PDB(s)")
        lines.append(f"  ✅ Matches: {matches}")
        lines.append(f"  ❌ Differences: {total-matches}")
        if total - matches > 0:
            lines.append("")
            lines.append("PDBs with differences:")
            for pdb_id, r in results.items():
                if r.has_differences():
                    lines.append(f"  - {pdb_id}")
        lines.append("")
        lines.append("TIP: Use --verbose with single PDB for details")
        lines.append("  Example: fp2-validate compare 1EHZ --verbose")
        report = "\n".join(lines)
    
    # Output
    if output:
        output.write_text(report)
        click.echo(f"Report saved to: {output}")
    else:
        click.echo(report)
    
    sys.exit(1 if any(r.has_differences() for r in results.values()) else 0)


if __name__ == '__main__':
    main()

