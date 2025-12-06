#!/usr/bin/env python3
"""
Unified CLI for X3DNA JSON validation and comparison.

Single entry point that replaces all the individual validate_*.py and test_*.py scripts.

Usage:
    fp2-validate                         # All stages, all fast PDBs
    fp2-validate frames                  # Frames only
    fp2-validate hbonds atoms            # Multiple stages
    fp2-validate --pdb 1EHZ --verbose    # Single PDB, detailed
    fp2-validate --stop-on-first         # Stop at first failure
    fp2-validate --diff                  # Document differences
    fp2-validate --quiet                 # Exit code only (for CI)
"""

import sys
from pathlib import Path

import click

from .runner import ValidationRunner
from .pdb_list import get_pdb_list, load_test_set


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
@click.option('--project-root', type=click.Path(path_type=Path, exists=True),
              help='Project root directory (default: current directory)')
def validate(stages, pdb, max_count, test_set, workers, quiet, verbose, 
             stop_on_first, diff, diff_file, project_root):
    """Validate legacy vs modern JSON outputs.
    
    STAGES: atoms, frames, hbonds, pairs, steps (default: all)
    
    \b
    Examples:
        fp2-validate                      # All stages, all fast PDBs
        fp2-validate frames               # Frames only
        fp2-validate hbonds atoms         # Multiple stages
        fp2-validate --pdb 1EHZ -v        # Single PDB, verbose
        fp2-validate --stop-on-first      # Debug mode
        fp2-validate --diff               # Document differences
        fp2-validate --max 100            # First 100 PDBs only
        fp2-validate --test-set 100       # Use saved test set
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


# Aliases for common operations
@main.command('atoms')
@click.option('--pdb', '-p', multiple=True)
@click.option('--max', '-n', 'max_count', type=int)
@click.option('--quiet', '-q', is_flag=True)
@click.option('--verbose', '-v', is_flag=True)
@click.option('--stop-on-first', '-s', is_flag=True)
@click.pass_context
def atoms(ctx, **kwargs):
    """Validate atom records (alias for: validate atoms)."""
    ctx.invoke(validate, stages=['atoms'], **kwargs)


@main.command('frames')
@click.option('--pdb', '-p', multiple=True)
@click.option('--max', '-n', 'max_count', type=int)
@click.option('--quiet', '-q', is_flag=True)
@click.option('--verbose', '-v', is_flag=True)
@click.option('--stop-on-first', '-s', is_flag=True)
@click.pass_context
def frames(ctx, **kwargs):
    """Validate frame calculations (alias for: validate frames)."""
    ctx.invoke(validate, stages=['frames'], **kwargs)


@main.command('hbonds')
@click.option('--pdb', '-p', multiple=True)
@click.option('--max', '-n', 'max_count', type=int)
@click.option('--quiet', '-q', is_flag=True)
@click.option('--verbose', '-v', is_flag=True)
@click.option('--stop-on-first', '-s', is_flag=True)
@click.pass_context
def hbonds(ctx, **kwargs):
    """Validate H-bond lists (alias for: validate hbonds)."""
    ctx.invoke(validate, stages=['hbonds'], **kwargs)


@main.command('pairs')
@click.option('--pdb', '-p', multiple=True)
@click.option('--max', '-n', 'max_count', type=int)
@click.option('--quiet', '-q', is_flag=True)
@click.option('--verbose', '-v', is_flag=True)
@click.option('--stop-on-first', '-s', is_flag=True)
@click.pass_context
def pairs(ctx, **kwargs):
    """Validate base pairs (alias for: validate pairs)."""
    ctx.invoke(validate, stages=['pairs'], **kwargs)


@main.command('steps')
@click.option('--pdb', '-p', multiple=True)
@click.option('--max', '-n', 'max_count', type=int)
@click.option('--quiet', '-q', is_flag=True)
@click.option('--verbose', '-v', is_flag=True)
@click.option('--stop-on-first', '-s', is_flag=True)
@click.pass_context
def steps(ctx, **kwargs):
    """Validate step parameters (alias for: validate steps)."""
    ctx.invoke(validate, stages=['steps'], **kwargs)


if __name__ == '__main__':
    main()

