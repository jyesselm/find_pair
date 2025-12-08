"""
Command-line interface for validation runner using Click.

Usage:
    python -m validation.cli 3 4 5           # Run stages 3, 4, 5
    python -m validation.cli frames          # Run frame stages
    python -m validation.cli all --max 50    # All stages, 50 PDBs
"""

import time
from pathlib import Path
from typing import List, Optional

import click

from .config import STAGES, STAGE_GROUPS, resolve_stages
from .runner import validate_stage, print_stage_result, StageResult


def _get_project_root() -> Path:
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


def _get_pdb_ids(json_dir: Path, max_pdbs: Optional[int] = None) -> List[str]:
    """
    Get list of PDB IDs from JSON directory.
    
    Scans any subdirectory (e.g., pdb_atoms/) for JSON files to find PDB IDs.
    """
    # Try to find PDBs from any stage subdirectory
    for stage_dir in ["pdb_atoms", "residue_indices", "frame_calc"]:
        subdir = json_dir / stage_dir
        if subdir.exists():
            pdbs = sorted(
                f.stem for f in subdir.iterdir()
                if f.suffix == ".json" and len(f.stem) == 4
            )
            if pdbs:
                return pdbs[:max_pdbs] if max_pdbs else pdbs
    
    return []


def _run_stages(
    stages: List[int],
    pdb_ids: List[str],
    json_dir: Path,
    verbose: bool,
    stop_on_failure: bool,
    num_workers: int = 1
) -> List[StageResult]:
    """Run validation for all stages."""
    results = []
    
    for stage_num in stages:
        start = time.time()
        result = validate_stage(
            stage_num, pdb_ids, json_dir,
            verbose, stop_on_failure, num_workers
        )
        elapsed = time.time() - start
        
        print_stage_result(result)
        click.echo(f"  Time: {elapsed:.2f}s")
        results.append(result)
    
    return results


def _print_summary(results: List[StageResult]) -> None:
    """Print overall summary."""
    click.echo(f"\n{'='*60}")
    click.echo("OVERALL SUMMARY")
    click.echo(f"{'='*60}")
    
    for result in results:
        status = "✅" if result.failed == 0 else "❌"
        click.echo(f"  {status} Stage {result.stage_num} ({result.stage_name}): "
                   f"{result.passed}/{result.total} passed")


@click.command()
@click.argument("stages", nargs=-1)
@click.option(
    "--json-dir", "-d",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="JSON directory (default: data/json)"
)
@click.option(
    "--max-pdbs", "-m",
    type=int,
    default=None,
    help="Maximum PDBs to process"
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Verbose output"
)
@click.option(
    "--stop-on-failure", "-s",
    is_flag=True,
    help="Stop on first failure"
)
@click.option(
    "--workers", "-w",
    type=int,
    default=1,
    help="Number of parallel workers (default: 1)"
)
def main(
    stages: tuple,
    json_dir: Optional[Path],
    max_pdbs: Optional[int],
    verbose: bool,
    stop_on_failure: bool,
    workers: int
) -> None:
    """
    Validate legacy vs modern JSON outputs.
    
    STAGES can be numbers (1-12), names (frame_calc), or groups (frames, pairs, all).
    
    Examples:
    
        python -m validation.cli 3 4 5
        
        python -m validation.cli frames --max-pdbs 50
        
        python -m validation.cli all -v
    """
    # Default to all stages if none specified
    stage_list = list(stages) if stages else ["all"]
    resolved = resolve_stages(stage_list)
    
    if not resolved:
        click.echo("No valid stages specified", err=True)
        raise SystemExit(1)
    
    # Default json directory
    if json_dir is None:
        json_dir = _get_project_root() / "data" / "json"
    
    pdb_ids = _get_pdb_ids(json_dir, max_pdbs)
    
    if not pdb_ids:
        click.echo(f"No PDBs found in {json_dir}", err=True)
        raise SystemExit(1)
    
    click.echo(f"Running validation for stages: {resolved}")
    click.echo(f"PDBs: {len(pdb_ids)}")
    if workers > 1:
        click.echo(f"Workers: {workers}")
    
    results = _run_stages(resolved, pdb_ids, json_dir, verbose, stop_on_failure, workers)
    _print_summary(results)
    
    all_passed = all(r.failed == 0 for r in results)
    raise SystemExit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
