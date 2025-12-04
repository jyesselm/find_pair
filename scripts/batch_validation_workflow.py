#!/usr/bin/env python3
"""
Batch validation workflow for staged X3DNA validation.

This script orchestrates:
1. Batch validation of specific stages (e.g., distance_checks, hbonds, etc.)
2. Result logging to data/validation_results/
3. Cleanup of passing JSON (optional)
4. Detailed mismatch reporting

Usage:
    python scripts/batch_validation_workflow.py --stage distance_checks
    python scripts/batch_validation_workflow.py --stage all --cleanup
    python scripts/batch_validation_workflow.py --list-pdbs
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
from dataclasses import dataclass, asdict

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from x3dna_json_compare.frame_comparison import (
    compare_frames,
    compare_ls_fitting,
)


@dataclass
class StageConfig:
    """Configuration for a validation stage."""
    name: str
    json_types: List[str]
    description: str
    comparison_function: Optional[str] = None


# Define all stages
STAGES = {
    "residue_indices": StageConfig(
        name="residue_indices",
        json_types=["residue_indices"],
        description="Residue index mapping",
        comparison_function="compare_residue_indices"
    ),
    "ls_fitting": StageConfig(
        name="ls_fitting",
        json_types=["ls_fitting"],
        description="Least-squares fitting data",
        comparison_function="compare_ls_fitting"
    ),
    "frames": StageConfig(
        name="frames",
        json_types=["base_frame_calc", "frame_calc"],
        description="Reference frame calculations",
        comparison_function="compare_frames"
    ),
    "distance_checks": StageConfig(
        name="distance_checks",
        json_types=["distance_checks"],
        description="Geometric distance measurements",
        comparison_function="compare_distance_checks"
    ),
    "hbonds": StageConfig(
        name="hbonds",
        json_types=["hbond_list"],
        description="Hydrogen bond detection",
        comparison_function="compare_hbonds"
    ),
    "pair_validation": StageConfig(
        name="pair_validation",
        json_types=["pair_validation"],
        description="Base pair validation decisions",
        comparison_function="compare_pair_validation"
    ),
    "pair_selection": StageConfig(
        name="pair_selection",
        json_types=["find_bestpair_selection", "base_pair"],
        description="Final base pair selection (PRIMARY OUTPUT)",
        comparison_function="compare_pair_selection"
    ),
}


@dataclass
class ValidationResult:
    """Result of validating a single PDB for a stage."""
    pdb_id: str
    stage: str
    matches: bool
    total_comparisons: int = 0
    mismatches: List[Dict] = None
    error: Optional[str] = None
    timestamp: str = None

    def __post_init__(self):
        if self.mismatches is None:
            self.mismatches = []
        if self.timestamp is None:
            self.timestamp = datetime.now().isoformat()


class BatchValidator:
    """Batch validation orchestrator."""

    def __init__(self, root_dir: Path, cleanup_on_pass: bool = False):
        self.root_dir = root_dir
        self.json_dir = root_dir / "data" / "json"
        self.legacy_json_dir = root_dir / "data" / "json_legacy"
        self.results_dir = root_dir / "data" / "validation_results"
        self.cleanup_on_pass = cleanup_on_pass

        # Create results directory
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # Setup logging
        self.setup_logging()

    def setup_logging(self):
        """Setup logging to both file and console."""
        log_file = self.results_dir / f"batch_validation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Batch validation started. Log file: {log_file}")

    def find_pdbs_for_stage(self, stage: StageConfig) -> List[str]:
        """Find all PDBs that have JSON for this stage."""
        pdbs = set()

        for json_type in stage.json_types:
            json_type_dir = self.json_dir / json_type
            if not json_type_dir.exists():
                continue

            for json_file in json_type_dir.glob("*.json"):
                pdbs.add(json_file.stem)

        return sorted(list(pdbs))

    def validate_pdb(self, pdb_id: str, stage: StageConfig) -> ValidationResult:
        """Validate a single PDB for a given stage."""
        self.logger.info(f"Validating {pdb_id} for stage {stage.name}")

        try:
            # Check if all required JSON files exist
            for json_type in stage.json_types:
                modern_file = self.json_dir / json_type / f"{pdb_id}.json"
                legacy_file = self.legacy_json_dir / json_type / f"{pdb_id}.json"

                if not modern_file.exists():
                    return ValidationResult(
                        pdb_id=pdb_id,
                        stage=stage.name,
                        matches=False,
                        error=f"Modern JSON missing: {modern_file}"
                    )

                if not legacy_file.exists():
                    return ValidationResult(
                        pdb_id=pdb_id,
                        stage=stage.name,
                        matches=False,
                        error=f"Legacy JSON missing: {legacy_file}"
                    )

            # Perform comparison based on stage type
            if stage.name == "ls_fitting":
                result = compare_ls_fitting(pdb_id)
            elif stage.name == "frames":
                result = compare_frames(pdb_id)
            else:
                # For now, just do a simple JSON comparison
                result = self._simple_json_compare(pdb_id, stage)

            if result.matches:
                self.logger.info(f"✅ {pdb_id}: PASS")
                if self.cleanup_on_pass:
                    self._cleanup_passing_json(pdb_id, stage)
            else:
                self.logger.warning(f"❌ {pdb_id}: FAIL - {len(result.mismatches)} mismatches")

            return ValidationResult(
                pdb_id=pdb_id,
                stage=stage.name,
                matches=result.matches,
                total_comparisons=result.total_comparisons,
                mismatches=result.mismatches[:10] if hasattr(result, 'mismatches') else []  # Limit to 10
            )

        except Exception as e:
            self.logger.error(f"Error validating {pdb_id}: {e}")
            return ValidationResult(
                pdb_id=pdb_id,
                stage=stage.name,
                matches=False,
                error=str(e)
            )

    def _simple_json_compare(self, pdb_id: str, stage: StageConfig) -> ValidationResult:
        """Simple JSON comparison for stages without specialized comparison functions."""
        all_match = True
        total_comparisons = 0
        mismatches = []

        for json_type in stage.json_types:
            modern_file = self.json_dir / json_type / f"{pdb_id}.json"
            legacy_file = self.legacy_json_dir / json_type / f"{pdb_id}.json"

            with open(modern_file) as f:
                modern = json.load(f)
            with open(legacy_file) as f:
                legacy = json.load(f)

            # Simple list length comparison
            if len(modern) != len(legacy):
                all_match = False
                mismatches.append({
                    "type": json_type,
                    "issue": f"Count mismatch: modern={len(modern)}, legacy={len(legacy)}"
                })
                continue

            total_comparisons += len(modern)

            # TODO: Implement detailed comparison for each JSON type
            # For now, just check counts

        result = ValidationResult(
            pdb_id=pdb_id,
            stage=stage.name,
            matches=all_match,
            total_comparisons=total_comparisons,
            mismatches=mismatches
        )
        return result

    def _cleanup_passing_json(self, pdb_id: str, stage: StageConfig):
        """Delete JSON files for a passing PDB."""
        for json_type in stage.json_types:
            modern_file = self.json_dir / json_type / f"{pdb_id}.json"
            legacy_file = self.legacy_json_dir / json_type / f"{pdb_id}.json"

            if modern_file.exists():
                modern_file.unlink()
                self.logger.info(f"  Deleted modern JSON: {modern_file}")

            if legacy_file.exists():
                legacy_file.unlink()
                self.logger.info(f"  Deleted legacy JSON: {legacy_file}")

    def validate_batch(self, pdbs: List[str], stage: StageConfig) -> List[ValidationResult]:
        """Validate a batch of PDBs for a given stage."""
        results = []

        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"Validating {len(pdbs)} PDBs for stage: {stage.name}")
        self.logger.info(f"Description: {stage.description}")
        self.logger.info(f"JSON types: {', '.join(stage.json_types)}")
        self.logger.info(f"{'='*60}\n")

        for i, pdb_id in enumerate(pdbs, 1):
            self.logger.info(f"[{i}/{len(pdbs)}] Processing {pdb_id}")
            result = self.validate_pdb(pdb_id, stage)
            results.append(result)

        return results

    def save_results(self, stage_name: str, results: List[ValidationResult]):
        """Save validation results to JSON file."""
        result_file = self.results_dir / f"{stage_name}_results.json"

        # Load existing results if any
        existing_results = []
        if result_file.exists():
            with open(result_file) as f:
                existing_results = json.load(f)

        # Merge new results (replace existing entries for same PDB)
        pdb_to_result = {r["pdb_id"]: r for r in existing_results}
        for result in results:
            pdb_to_result[result.pdb_id] = asdict(result)

        # Save merged results
        with open(result_file, 'w') as f:
            json.dump(list(pdb_to_result.values()), f, indent=2)

        self.logger.info(f"\n✅ Results saved to: {result_file}")

    def print_summary(self, stage_name: str, results: List[ValidationResult]):
        """Print validation summary."""
        total = len(results)
        passed = sum(1 for r in results if r.matches)
        failed = sum(1 for r in results if not r.matches)
        errors = sum(1 for r in results if r.error is not None)

        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"SUMMARY: {stage_name}")
        self.logger.info(f"{'='*60}")
        self.logger.info(f"Total PDFs: {total}")
        self.logger.info(f"✅ Passed: {passed} ({100*passed/total:.1f}%)")
        self.logger.info(f"❌ Failed: {failed} ({100*failed/total:.1f}%)")
        if errors > 0:
            self.logger.info(f"⚠️  Errors: {errors}")

        # List failures
        if failed > 0:
            self.logger.info(f"\nFailed PDBs:")
            for r in results:
                if not r.matches:
                    error_msg = r.error if r.error else f"{len(r.mismatches)} mismatches"
                    self.logger.info(f"  - {r.pdb_id}: {error_msg}")

        self.logger.info(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Batch validation workflow for staged X3DNA validation"
    )
    parser.add_argument(
        "--stage",
        choices=list(STAGES.keys()) + ["all"],
        default="distance_checks",
        help="Stage to validate (default: distance_checks)"
    )
    parser.add_argument(
        "--pdbs",
        nargs="+",
        help="Specific PDBs to validate (default: all available for stage)"
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Delete passing JSON files after validation"
    )
    parser.add_argument(
        "--list-pdbs",
        action="store_true",
        help="List available PDBs for each stage and exit"
    )

    args = parser.parse_args()

    # Find root directory
    root_dir = Path(__file__).parent.parent
    validator = BatchValidator(root_dir, cleanup_on_pass=args.cleanup)

    # List PDFs mode
    if args.list_pdbs:
        print("\nAvailable PDBs by stage:")
        print("=" * 60)
        for stage_name, stage_config in STAGES.items():
            pdbs = validator.find_pdbs_for_stage(stage_config)
            print(f"{stage_name:20s} {len(pdbs):4d} PDBs")
            if len(pdbs) <= 10:
                print(f"  {', '.join(pdbs)}")
        print("=" * 60)
        return 0

    # Determine which stages to validate
    if args.stage == "all":
        stages_to_validate = list(STAGES.values())
    else:
        stages_to_validate = [STAGES[args.stage]]

    # Run validation for each stage
    for stage in stages_to_validate:
        # Find PDBs to validate
        if args.pdbs:
            pdbs = args.pdbs
        else:
            pdbs = validator.find_pdbs_for_stage(stage)

        if not pdbs:
            validator.logger.warning(f"No PDBs found for stage: {stage.name}")
            continue

        # Validate batch
        results = validator.validate_batch(pdbs, stage)

        # Save and summarize results
        validator.save_results(stage.name, results)
        validator.print_summary(stage.name, results)

    return 0


if __name__ == "__main__":
    sys.exit(main())

