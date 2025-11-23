#!/usr/bin/env python3
"""
Generate detailed human-readable report with PDB content examples.

This script runs all Phase 1 analyses and creates a comprehensive report
with actual PDB atom lines for examples.

Usage:
    python3 scripts/generate_detailed_report.py
    python3 scripts/generate_detailed_report.py --output detailed_report.md
"""

import json
import sys
import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional

from x3dna_json_compare import JsonComparator, ComparisonResult


class PdbReader:
    """Simple PDB file reader to extract atom lines."""

    def __init__(self, pdb_file: Path):
        self.pdb_file = pdb_file
        self.atoms = {}
        self._load_atoms()

    def _load_atoms(self):
        """Load all ATOM records from PDB file."""
        if not self.pdb_file.exists():
            return

        try:
            with open(self.pdb_file, "r") as f:
                for line in f:
                    if line.startswith("ATOM  ") or line.startswith("HETATM"):
                        atom_name = line[12:16].strip()
                        residue_name = line[17:20].strip()
                        chain_id = line[21:22]
                        residue_seq = int(line[22:26].strip())
                        insertion = line[26:27] if len(line) > 26 else " "

                        key = (chain_id, residue_seq, insertion, atom_name)
                        self.atoms[key] = {
                            "line": line.rstrip(),
                            "atom_name": atom_name,
                            "residue_name": residue_name,
                            "chain_id": chain_id,
                            "residue_seq": residue_seq,
                            "insertion": insertion,
                            "x": float(line[30:38].strip()) if len(line) > 38 else 0.0,
                            "y": float(line[38:46].strip()) if len(line) > 46 else 0.0,
                            "z": float(line[46:54].strip()) if len(line) > 54 else 0.0,
                        }
        except Exception as e:
            print(f"Warning: Could not read PDB file {self.pdb_file}: {e}")

    def get_atom_line(
        self, chain_id: str, residue_seq: int, insertion: str, atom_name: str
    ) -> Optional[str]:
        """Get PDB line for a specific atom."""
        # Normalize atom name (add spaces if needed)
        normalized_name = f" {atom_name.strip()} "
        if len(normalized_name) < 4:
            normalized_name = normalized_name[:4]

        key = (chain_id, residue_seq, insertion, normalized_name)
        if key in self.atoms:
            return self.atoms[key]["line"]

        # Try without spaces
        key2 = (chain_id, residue_seq, insertion, atom_name.strip())
        if key2 in self.atoms:
            return self.atoms[key2]["line"]

        # Try with different spacing
        for stored_key, atom_data in self.atoms.items():
            if (
                stored_key[0] == chain_id
                and stored_key[1] == residue_seq
                and stored_key[2] == insertion
                and stored_key[3].strip() == atom_name.strip()
            ):
                return atom_data["line"]

        return None

    def get_residue_atoms(
        self, chain_id: str, residue_seq: int, insertion: str
    ) -> List[Dict]:
        """Get all atoms for a residue."""
        atoms = []
        for key, atom_data in self.atoms.items():
            if key[0] == chain_id and key[1] == residue_seq and key[2] == insertion:
                atoms.append(atom_data)
        return sorted(atoms, key=lambda x: x["atom_name"])


def extract_residue_info(residue_str: str) -> tuple:
    """Extract chain_id, residue_seq, insertion from residue string like 'A:123 '."""
    parts = residue_str.split(":")
    if len(parts) != 2:
        return None, None, None

    chain_id = parts[0]
    seq_part = parts[1]

    # Extract sequence number and insertion code
    residue_seq = int(seq_part.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ "))
    insertion = seq_part[-1] if seq_part[-1].isalpha() else " "

    return chain_id, residue_seq, insertion


def generate_report(project_root: Path, output_file: Path) -> None:
    """Generate comprehensive detailed report."""
    print("=" * 80)
    print("GENERATING DETAILED ANALYSIS REPORT")
    print("=" * 80)
    print()

    # Run analyses and collect data
    print("Step 1: Analyzing base type differences...")
    from analyze_base_type_differences import analyze_base_type
    import tempfile
    import subprocess

    # We'll collect data by running comparisons directly
    json_dir = project_root / "data/json"
    legacy_dir = project_root / "data/json_legacy"

    modern_files = set(
        f.stem for f in json_dir.glob("*.json") if not f.stem.endswith("_legacy")
    )
    legacy_files = set(f.stem for f in legacy_dir.glob("*.json"))
    common = sorted(modern_files & legacy_files)

    comparator = JsonComparator(enable_cache=False)

    # Collect examples for each base type
    examples_by_base = defaultdict(list)
    base_stats = defaultdict(
        lambda: {
            "total": 0,
            "exact_matches": 0,
            "real_diffs": 0,
            "patterns": defaultdict(int),
        }
    )

    print(f"Processing {len(common)} files to collect examples...")

    import multiprocessing
    from concurrent.futures import ThreadPoolExecutor, as_completed

    max_workers = multiprocessing.cpu_count()

    def compare_pdb(pdb_id):
        """Compare and collect examples."""
        legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
        modern_file = project_root / f"data/json/{pdb_id}.json"
        pdb_file = project_root / f"data/pdb/{pdb_id}.pdb"

        if not all([legacy_file.exists(), modern_file.exists()]):
            return []

        try:
            result = comparator.compare_files(
                legacy_file, modern_file, pdb_file, pdb_id
            )

            pdb_examples = []

            if result.frame_comparison:
                fc = result.frame_comparison

                for mismatch in fc.mismatched_calculations:
                    base_type = mismatch.legacy_record.get("base_type", "?")
                    leg_atoms = set(mismatch.legacy_matched_atoms)
                    mod_atoms = set(mismatch.modern_matched_atoms)
                    chain_id, residue_seq, insertion = mismatch.residue_key

                    # Collect all mismatches, not just atom set differences
                    # This includes order-only differences too
                    only_legacy = sorted(leg_atoms - mod_atoms)
                    only_modern = sorted(mod_atoms - leg_atoms)

                    # Only include if there are real atom set differences
                    # or if RMS difference is significant
                    rms_diff = abs(
                        mismatch.legacy_record.get("rms_fit", 0)
                        - mismatch.modern_record.get("rms_fit", 0)
                    )

                    if leg_atoms != mod_atoms or rms_diff > 0.001:
                        example = {
                            "pdb_id": pdb_id,
                            "residue": f"{chain_id}:{residue_seq}{insertion}",
                            "chain_id": chain_id,
                            "residue_seq": residue_seq,
                            "insertion": insertion,
                            "residue_name": mismatch.legacy_record.get(
                                "residue_name", ""
                            ),
                            "base_type": base_type,
                            "legacy_atoms": sorted(leg_atoms),
                            "modern_atoms": sorted(mod_atoms),
                            "only_legacy": only_legacy,
                            "only_modern": only_modern,
                            "rms_legacy": mismatch.legacy_record.get("rms_fit", 0),
                            "rms_modern": mismatch.modern_record.get("rms_fit", 0),
                            "rms_diff": abs(
                                mismatch.legacy_record.get("rms_fit", 0)
                                - mismatch.modern_record.get("rms_fit", 0)
                            ),
                            "pdb_file": pdb_file,
                        }

                        pdb_examples.append(example)

            return pdb_examples
        except Exception as e:
            return []

    # Process in parallel
    completed = 0
    all_examples = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {
            executor.submit(compare_pdb, pdb_id): pdb_id for pdb_id in common
        }

        for future in as_completed(future_to_pdb):
            pdb_examples = future.result()
            all_examples.extend(pdb_examples)
            completed += 1

            if completed % 100 == 0:
                print(
                    f"  Progress: {completed}/{len(common)} ({completed*100//len(common)}%)",
                    flush=True,
                )

    # Group by base type and limit examples per base type
    max_examples_per_base = 50  # Limit to keep report manageable

    for example in all_examples:
        base_type = example["base_type"]

        # Only count real differences (atom set differences)
        if example["only_legacy"] or example["only_modern"]:
            base_stats[base_type]["real_diffs"] += 1

            # Track pattern
            pattern = (tuple(example["only_legacy"]), tuple(example["only_modern"]))
            base_stats[base_type]["patterns"][pattern] += 1

        # Add to examples (limit per base type)
        if len(examples_by_base[base_type]) < max_examples_per_base:
            examples_by_base[base_type].append(example)

        base_stats[base_type]["total"] += 1

    print(
        f"\nCollected {len(all_examples)} examples across {len(examples_by_base)} base types"
    )
    print()

    # Generate report
    print("Step 2: Generating detailed report with PDB content...")

    report_lines = []
    report_lines.append("# Detailed Analysis Report with PDB Content")
    report_lines.append("")
    report_lines.append(
        f"**Generated:** {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    )
    report_lines.append(f"**Total Examples Analyzed:** {len(all_examples)}")
    report_lines.append("")

    # Priority order: A, a, G, g, then others
    priority_bases = ["A", "a", "G", "g", "C", "c", "T", "t", "U", "u"]
    other_bases = sorted(set(examples_by_base.keys()) - set(priority_bases))
    ordered_bases = [b for b in priority_bases if b in examples_by_base] + other_bases

    for base_type in ordered_bases:
        examples = examples_by_base[base_type]
        stats = base_stats[base_type]

        report_lines.append("=" * 80)
        report_lines.append(f"BASE TYPE: {base_type}")
        report_lines.append("=" * 80)
        report_lines.append("")
        report_lines.append(f"**Total Differences:** {stats['real_diffs']}")
        report_lines.append(f"**Unique Patterns:** {len(stats['patterns'])}")
        report_lines.append("")

        # Show top patterns
        sorted_patterns = sorted(
            stats["patterns"].items(), key=lambda x: x[1], reverse=True
        )

        report_lines.append("### Top 5 Difference Patterns")
        report_lines.append("")

        for i, (pattern, count) in enumerate(sorted_patterns[:5], 1):
            only_legacy, only_modern = pattern
            report_lines.append(
                f"**Pattern {i}:** {count} occurrences ({count*100/stats['real_diffs']:.1f}%)"
            )
            if only_legacy:
                report_lines.append(f"- Only in legacy: `{list(only_legacy)}`")
            if only_modern:
                report_lines.append(f"- Only in modern: `{list(only_modern)}`")
            report_lines.append("")

        # Show detailed examples with PDB content
        report_lines.append("### Detailed Examples with PDB Content")
        report_lines.append("")

        # Get top 10 examples by RMS difference or pattern frequency
        examples_to_show = sorted(
            examples,
            key=lambda x: (
                x["rms_diff"],
                -len(x["only_legacy"]) - len(x["only_modern"]),
            ),
            reverse=True,
        )[:10]

        for i, example in enumerate(examples_to_show, 1):
            pdb_file = example["pdb_file"]
            chain_id = example["chain_id"]
            residue_seq = example["residue_seq"]
            insertion = example["insertion"]

            report_lines.append(
                f"#### Example {i}: {example['pdb_id']} {example['residue']}"
            )
            report_lines.append("")
            report_lines.append(
                f"- **Residue:** {example['residue_name']} {example['residue']}"
            )
            report_lines.append(f"- **Base Type:** {example['base_type']}")
            report_lines.append(f"- **RMS Difference:** {example['rms_diff']:.6f} Å")
            report_lines.append("")

            # Legacy atoms
            report_lines.append("**Legacy Matched Atoms:**")
            report_lines.append(f"```")
            report_lines.append(f"{', '.join(example['legacy_atoms'])}")
            report_lines.append(f"```")
            report_lines.append("")

            # Modern atoms
            report_lines.append("**Modern Matched Atoms:**")
            report_lines.append(f"```")
            report_lines.append(f"{', '.join(example['modern_atoms'])}")
            report_lines.append(f"```")
            report_lines.append("")

            # Differences
            if example["only_legacy"]:
                report_lines.append(f"**Only in Legacy:** `{example['only_legacy']}`")
                report_lines.append("")
            if example["only_modern"]:
                report_lines.append(f"**Only in Modern:** `{example['only_modern']}`")
                report_lines.append("")

            # PDB content
            if pdb_file.exists():
                pdb_reader = PdbReader(pdb_file)
                residue_atoms = pdb_reader.get_residue_atoms(
                    chain_id, residue_seq, insertion
                )

                if residue_atoms:
                    report_lines.append("**PDB Atom Lines for this Residue:**")
                    report_lines.append("```")

                    # Show all atoms in residue
                    for atom in residue_atoms:
                        atom_name = atom["atom_name"]
                        in_legacy = any(atom_name in a for a in example["legacy_atoms"])
                        in_modern = any(atom_name in a for a in example["modern_atoms"])

                        marker = ""
                        if atom_name in [a.strip() for a in example["only_legacy"]]:
                            marker = " ← ONLY LEGACY"
                        elif atom_name in [a.strip() for a in example["only_modern"]]:
                            marker = " ← ONLY MODERN"
                        elif in_legacy and in_modern:
                            marker = " ← IN BOTH"
                        elif in_legacy:
                            marker = " ← IN LEGACY"
                        elif in_modern:
                            marker = " ← IN MODERN"

                        report_lines.append(f"{atom['line']}{marker}")

                    report_lines.append("```")
                    report_lines.append("")

                    # Show specific atom lines for differences
                    if example["only_legacy"] or example["only_modern"]:
                        report_lines.append("**PDB Lines for Differing Atoms:**")
                        report_lines.append("```")

                        for atom_name in example["only_legacy"]:
                            atom_line = pdb_reader.get_atom_line(
                                chain_id, residue_seq, insertion, atom_name
                            )
                            if atom_line:
                                report_lines.append(f"{atom_line} ← ONLY IN LEGACY")

                        for atom_name in example["only_modern"]:
                            atom_line = pdb_reader.get_atom_line(
                                chain_id, residue_seq, insertion, atom_name
                            )
                            if atom_line:
                                report_lines.append(f"{atom_line} ← ONLY IN MODERN")

                        report_lines.append("```")
                        report_lines.append("")
                else:
                    report_lines.append("*PDB file exists but residue not found*")
                    report_lines.append("")
            else:
                report_lines.append(f"*PDB file not found: {pdb_file}*")
                report_lines.append("")

            report_lines.append("---")
            report_lines.append("")

    # Write report
    report_content = "\n".join(report_lines)
    output_file.write_text(report_content)

    print(f"\n✅ Detailed report saved to: {output_file}")
    print(f"   Report contains {len(all_examples)} examples with PDB content")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Generate detailed report with PDB content"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="docs/DETAILED_ANALYSIS_REPORT.md",
        help="Output file path (default: docs/DETAILED_ANALYSIS_REPORT.md)",
    )

    args = parser.parse_args()

    project_root = Path(__file__).parent.parent
    output_file = project_root / args.output

    # Create docs directory if needed
    output_file.parent.mkdir(exist_ok=True)

    generate_report(project_root, output_file)


if __name__ == "__main__":
    main()
