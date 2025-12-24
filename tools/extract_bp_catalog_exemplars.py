#!/usr/bin/env python3
"""
Extract and align base pair exemplars from the BP database catalog.

Parses data/bp_database.txt, fetches PDB files, applies symmetry operators,
calculates base pair reference frames, and outputs aligned PDB files.

Usage:
    python tools/extract_bp_catalog_exemplars.py
    python tools/extract_bp_catalog_exemplars.py --skip-modeled -v
    python tools/extract_bp_catalog_exemplars.py --lw-class cWW --verbose
"""

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional

# Add tools directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from bp_catalog_parser import load_database, ExemplarEntry
from bp_catalog_extractor import BasePairExemplarExtractor


def extract_entry(
    extractor: BasePairExemplarExtractor,
    entry: ExemplarEntry,
    verbose: bool
) -> Optional[Dict]:
    """Extract a single entry, returning metadata dict or None."""
    try:
        result = extractor.extract(entry, verbose=verbose)
        if result:
            return {
                'lw_class': entry.lw_class,
                'sequence': entry.sequence,
                'bp_key': entry.bp_key,
                'iso_group': entry.iso_group,
                'pdb_id': entry.pdb_id,
                'chain1': entry.chain1,
                'chain2': entry.chain2,
                'resnum1': entry.resnum1,
                'resnum2': entry.resnum2,
                'resolution': entry.resolution,
                'needs_symmetry': entry.needs_symmetry,
                'output_file': str(result.relative_to(result.parent.parent))
            }
    except Exception as e:
        if verbose:
            print(f"  Error extracting {entry}: {e}")
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Extract and align base pair exemplars from BP database catalog"
    )
    parser.add_argument(
        "--database", "-d",
        default="data/bp_database.txt",
        help="Path to bp_database.txt (default: data/bp_database.txt)"
    )
    parser.add_argument(
        "--pdb-dir",
        default="data/pdb",
        help="Directory containing PDB files (default: data/pdb)"
    )
    parser.add_argument(
        "--template-dir",
        default="resources/templates",
        help="Directory containing base templates (default: resources/templates)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="basepair-catalog-exemplars",
        help="Output directory for exemplar PDBs (default: basepair-catalog-exemplars)"
    )
    parser.add_argument(
        "--lw-class",
        help="Only process specific LW class (e.g., cWW, tHS)"
    )
    parser.add_argument(
        "--skip-modeled",
        action="store_true",
        help="Skip modeled (non-experimental) entries"
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Don't download missing PDB files"
    )
    parser.add_argument(
        "--workers", "-w",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1, sequential)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Just print what would be done, don't extract"
    )

    args = parser.parse_args()

    # Resolve paths relative to project root
    project_root = Path(__file__).parent.parent
    db_path = project_root / args.database
    pdb_dir = project_root / args.pdb_dir
    template_dir = project_root / args.template_dir
    output_dir = project_root / args.output_dir

    # Load database
    print(f"Loading database from {db_path}...")
    entries = load_database(db_path)
    total_entries = sum(len(v) for v in entries.values())
    print(f"Loaded {total_entries} entries across {len(entries)} LW classes")

    # Filter by LW class if specified
    if args.lw_class:
        if args.lw_class not in entries:
            print(f"Error: LW class '{args.lw_class}' not found in database")
            print(f"Available: {', '.join(sorted(entries.keys()))}")
            return 1
        entries = {args.lw_class: entries[args.lw_class]}

    # Collect entries to process
    to_process: List[ExemplarEntry] = []
    for lw_class, entry_list in sorted(entries.items()):
        for entry in entry_list:
            if args.skip_modeled and entry.is_modeled:
                continue
            to_process.append(entry)

    print(f"\nProcessing {len(to_process)} entries...")

    if args.dry_run:
        print("\nDry run - entries to process:")
        for entry in to_process:
            sym_note = " [symmetry]" if entry.needs_symmetry else ""
            print(f"  {entry.lw_class}/{entry.sequence}: {entry.pdb_id} "
                  f"{entry.chain1}.{entry.resnum1}-{entry.chain2}.{entry.resnum2}{sym_note}")
        return 0

    # Create extractor
    extractor = BasePairExemplarExtractor(
        pdb_dir=pdb_dir,
        template_dir=template_dir,
        output_dir=output_dir,
        download_missing=not args.no_download
    )

    # Process entries
    output_dir.mkdir(parents=True, exist_ok=True)
    results: List[Dict] = []
    extracted = 0
    failed = 0

    if args.workers > 1:
        # Parallel extraction (note: file downloads are sequential internally)
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(extract_entry, extractor, entry, args.verbose): entry
                for entry in to_process
            }

            for i, future in enumerate(as_completed(futures), 1):
                entry = futures[future]
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                        extracted += 1
                        if args.verbose:
                            print(f"  [{i}/{len(to_process)}] Extracted: {entry.lw_class}/{entry.sequence}")
                    else:
                        failed += 1
                except Exception as e:
                    failed += 1
                    if args.verbose:
                        print(f"  [{i}/{len(to_process)}] Failed: {entry} - {e}")

                if i % 20 == 0 and not args.verbose:
                    print(f"  Processed {i}/{len(to_process)}...")
    else:
        # Sequential extraction
        for i, entry in enumerate(to_process, 1):
            result = extract_entry(extractor, entry, args.verbose)
            if result:
                results.append(result)
                extracted += 1
                if args.verbose:
                    print(f"  [{i}/{len(to_process)}] Extracted: {entry.lw_class}/{entry.sequence}")
            else:
                failed += 1

            if i % 20 == 0 and not args.verbose:
                print(f"  Processed {i}/{len(to_process)}...")

    # Write summary JSON
    summary_file = output_dir / "summary.json"
    with open(summary_file, 'w') as f:
        json.dump({
            'total_entries': len(to_process),
            'extracted': extracted,
            'failed': failed,
            'exemplars': results
        }, f, indent=2)

    # Print summary by LW class
    print(f"\n{'='*50}")
    print("Summary by LW class:")
    print(f"{'LW Class':<10} {'Extracted':<12} {'Total':<10}")
    print("-" * 35)

    lw_counts: Dict[str, Dict[str, int]] = {}
    for entry in to_process:
        if entry.lw_class not in lw_counts:
            lw_counts[entry.lw_class] = {'total': 0, 'extracted': 0}
        lw_counts[entry.lw_class]['total'] += 1

    for result in results:
        lw_counts[result['lw_class']]['extracted'] += 1

    for lw_class in sorted(lw_counts.keys()):
        counts = lw_counts[lw_class]
        print(f"{lw_class:<10} {counts['extracted']:<12} {counts['total']:<10}")

    print("-" * 35)
    print(f"{'TOTAL':<10} {extracted:<12} {len(to_process):<10}")

    print(f"\nOutput directory: {output_dir}")
    print(f"Summary file: {summary_file}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
