"""
Output formatting utilities.

Controls all output in quiet/verbose/normal modes to avoid extra output.
"""

import sys
from dataclasses import dataclass
from typing import List, Optional, TextIO
from pathlib import Path
from .verbose_reporter import VerboseReporter


@dataclass
class ValidationSummary:
    """Summary of validation results."""
    total: int
    passed: int
    failed: int
    skipped: int
    stages_tested: List[str]
    differences_file: Optional[Path] = None
    first_failure_pdb: Optional[str] = None
    checkpoint_file: Optional[Path] = None
    
    @property
    def all_passed(self) -> bool:
        return self.failed == 0
    
    @property 
    def success_rate(self) -> float:
        testable = self.total - self.skipped
        return self.passed / testable if testable > 0 else 0.0


class OutputFormatter:
    """Unified output control - quiet/verbose/normal modes.
    
    - quiet: No output, only exit code
    - verbose: Per-PDB results
    - normal (default): Summary only, progress every 100
    """
    
    def __init__(self, quiet: bool = False, verbose: bool = False, 
                 stream: TextIO = None):
        self.quiet = quiet
        self.verbose = verbose
        self.stream = stream or sys.stdout
        self._current_failures = []
    
    def _print(self, msg: str = "", flush: bool = False):
        """Print if not quiet."""
        if not self.quiet:
            print(msg, file=self.stream, flush=flush)
    
    def start(self, total: int, stages: Optional[List[str]] = None):
        """Print start message."""
        if self.quiet:
            return
        self._print(f"Validating {total} PDBs...", flush=True)
        if stages and stages != ['all']:
            self._print(f"Stages: {', '.join(stages)}")
    
    def progress(self, current: int, total: int, pdb_id: str, result):
        """Print progress update."""
        if self.quiet:
            return
        
        # Handle both dict and object results
        status_val = result.get('status') if isinstance(result, dict) else result.status
            
        if self.verbose:
            # Show detailed comparison for each PDB
            if status_val == 'match':
                status = '✅'
            elif status_val == 'skip':
                status = '⏭️ '
            else:
                status = '❌'
                self._current_failures.append(pdb_id)
            
            self._print(f"\n{status} {pdb_id}")
            
            # Show detailed verbose report if available
            full_result = result.get('full_result') if isinstance(result, dict) else None
            
            if full_result:
                # Generate verbose report for detailed comparison
                try:
                    reporter = VerboseReporter(pdb_id, full_result)
                    verbose_output = reporter.generate_report()
                    
                    # VerboseReporter returns empty string for single-stage comparisons
                    # TODO: Enhance VerboseReporter to handle single-stage comparisons
                    if verbose_output:
                        self._print(verbose_output)
                    else:
                        # Fallback: show basic comparison info
                        self._print(f"  For detailed comparison: fp2-validate compare {pdb_id} --verbose")
                except Exception as e:
                    # Fallback if VerboseReporter fails
                    self._print(f"  For detailed comparison: fp2-validate compare {pdb_id} --verbose")
            elif status_val == 'match':
                self._print(f"  All comparisons passed")
            else:
                # Show basic info for non-verbose mode
                self._print(f"  Status: {status_val}")
        else:
            # Summary mode - progress every 100 or at end
            if current % 100 == 0 or current == total:
                self._print(f"Progress: {current}/{total}", flush=True)
    
    def first_failure(self, pdb_id: str, result):
        """Print details of first failure when --stop-on-first."""
        if self.quiet:
            return
            
        self._print()
        self._print("=" * 60)
        self._print(f"FIRST FAILURE: {pdb_id}")
        self._print("=" * 60)
        
        # Handle both dict and object results
        if isinstance(result, dict):
            errors = result.get('errors', [])
            atom_comp = result.get('atom_comparison')
            frame_comp = result.get('frame_comparison')
            hbond_comp = result.get('hbond_comparison')
        else:
            errors = result.errors if hasattr(result, 'errors') else []
            atom_comp = getattr(result, 'atom_comparison', None)
            frame_comp = getattr(result, 'frame_comparison', None)
            hbond_comp = getattr(result, 'hbond_list_comparison', None)
        
        if errors:
            self._print("Errors:")
            for err in errors:
                self._print(f"  - {err}")
        
        # Show comparison details (for dict results, these are summaries)
        if atom_comp:
            if isinstance(atom_comp, dict):
                if atom_comp.get('missing_in_modern') or atom_comp.get('extra_in_modern'):
                    self._print(f"\nAtom differences:")
                    self._print(f"  Missing in modern: {atom_comp.get('missing_in_modern', 0)}")
                    self._print(f"  Extra in modern: {atom_comp.get('extra_in_modern', 0)}")
            elif hasattr(atom_comp, 'has_differences') and atom_comp.has_differences():
                self._print(f"\nAtom differences:")
                self._print(f"  Missing in modern: {len(atom_comp.missing_in_modern)}")
                self._print(f"  Extra in modern: {len(atom_comp.extra_in_modern)}")
                self._print(f"  Mismatched fields: {len(atom_comp.mismatched_fields)}")
        
        if frame_comp:
            if isinstance(frame_comp, dict):
                if frame_comp.get('missing_residues') or frame_comp.get('mismatched_calculations'):
                    self._print(f"\nFrame differences:")
                    self._print(f"  Missing residues: {frame_comp.get('missing_residues', 0)}")
                    self._print(f"  Mismatched: {frame_comp.get('mismatched_calculations', 0)}")
            elif frame_comp.mismatched_calculations or frame_comp.missing_residues:
                self._print(f"\nFrame differences:")
                self._print(f"  Missing residues: {len(frame_comp.missing_residues)}")
                self._print(f"  Mismatched: {len(frame_comp.mismatched_calculations)}")
        
        if hbond_comp:
            if isinstance(hbond_comp, dict):
                if hbond_comp.get('missing_in_modern') or hbond_comp.get('extra_in_modern'):
                    self._print(f"\nH-bond differences:")
                    self._print(f"  Missing in modern: {hbond_comp.get('missing_in_modern', 0)}")
                    self._print(f"  Extra in modern: {hbond_comp.get('extra_in_modern', 0)}")
            elif hbond_comp.missing_in_modern or hbond_comp.extra_in_modern:
                self._print(f"\nH-bond differences:")
                self._print(f"  Missing in modern: {len(hbond_comp.missing_in_modern)}")
                self._print(f"  Extra in modern: {len(hbond_comp.extra_in_modern)}")
        
        self._print()
        self._print("Use --verbose for more details or investigate manually:")
        self._print(f"  fp2-validate --pdb {pdb_id} --verbose")
    
    def summary(self, summary: ValidationSummary):
        """Print final summary."""
        if self.quiet:
            return
            
        self._print()
        self._print("=" * 60)
        self._print("VALIDATION SUMMARY")
        self._print("=" * 60)
        self._print(f"Total: {summary.total}")
        self._print(f"Passed: {summary.passed} ({summary.success_rate:.1%})")
        self._print(f"Failed: {summary.failed}")
        if summary.skipped > 0:
            self._print(f"Skipped: {summary.skipped}")
        
        if summary.differences_file:
            self._print(f"\nDifferences saved to: {summary.differences_file}")
        
        if summary.checkpoint_file:
            self._print(f"Checkpoint saved to: {summary.checkpoint_file}")
        
        if summary.first_failure_pdb:
            self._print(f"\nStopped at first failure: {summary.first_failure_pdb}")
        
        # Show failures in non-verbose mode
        if not self.verbose and self._current_failures:
            self._print(f"\nFailed PDBs ({len(self._current_failures)}):")
            for pdb_id in self._current_failures[:20]:
                self._print(f"  - {pdb_id}")
            if len(self._current_failures) > 20:
                self._print(f"  ... and {len(self._current_failures) - 20} more")
        
        if summary.all_passed:
            self._print("\n✅ All validations passed!")
        else:
            self._print(f"\n❌ {summary.failed} validation(s) failed")

