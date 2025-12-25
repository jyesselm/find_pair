"""Example usage of the cWW miss annotator data loaders.

Demonstrates how to load DSSR pairs and slot H-bonds for comparison.
"""

from pathlib import Path

from loaders import (
    DSSRPair,
    SlotHBond,
    dssr_to_slot_id,
    load_dssr_pairs,
    load_slot_hbonds,
    slot_to_dssr_id,
)


def main():
    """Example workflow for loading and comparing data."""

    # Example PDB ID
    pdb_id = "1A9N"

    # Get absolute path to data directory
    base_dir = Path(__file__).parent.parent.parent.parent

    # Load DSSR pairs
    dssr_path = base_dir / "data" / "json_dssr" / f"{pdb_id}.json"
    print(f"Loading DSSR pairs from {dssr_path}")
    dssr_pairs = load_dssr_pairs(dssr_path, lw_filter="cWW")
    print(f"Found {len(dssr_pairs)} cWW pairs\n")

    # Load slot H-bonds
    hbond_path = base_dir / "data" / "json" / "slot_hbonds" / f"{pdb_id}.json"
    print(f"Loading slot H-bonds from {hbond_path}")
    slot_hbonds = load_slot_hbonds(hbond_path)
    unique_pairs = len(slot_hbonds) // 2  # Bidirectional storage
    print(f"Found {unique_pairs} unique residue pairs with H-bonds\n")

    # Example: Check H-bonds for each DSSR pair
    print("Checking H-bonds for DSSR pairs:")
    print("-" * 80)

    for (nt1, nt2), dssr_pair in list(dssr_pairs.items())[:5]:
        print(f"\nDSSR Pair: {nt1} - {nt2}")
        print(f"  Base pair: {dssr_pair.bp}")
        print(f"  Saenger: {dssr_pair.saenger}")
        print(f"  DSSR H-bonds: {dssr_pair.hbonds_desc}")
        print(f"  Interbase angle: {dssr_pair.interbase_angle:.1f}°")
        print(f"  N1-N9 distance: {dssr_pair.n1n9_distance:.2f}Å")

        # Check slot H-bonds
        if (nt1, nt2) in slot_hbonds:
            hbond_list = slot_hbonds[(nt1, nt2)]
            print(f"  Slot H-bonds ({len(hbond_list)} total):")

            # Filter to base-base H-bonds only
            base_hbonds = [hb for hb in hbond_list if hb.context == "base_base"]
            for hb in base_hbonds:
                print(f"    {hb.donor_atom} -> {hb.acceptor_atom}: "
                      f"{hb.distance:.2f}Å (h_slot={hb.h_slot}, lp_slot={hb.lp_slot})")
        else:
            print("  No slot H-bonds found")

    # Example: Find residue pairs with H-bonds but not in DSSR
    print("\n" + "=" * 80)
    print("Residue pairs with H-bonds NOT in DSSR cWW pairs:")
    print("-" * 80)

    shown = 0
    for (res_i, res_j) in slot_hbonds.keys():
        if shown >= 5:
            break
        # Only check one direction
        if res_i < res_j:
            if (res_i, res_j) not in dssr_pairs:
                hbond_list = slot_hbonds[(res_i, res_j)]
                base_hbonds = [hb for hb in hbond_list if hb.context == "base_base"]
                if base_hbonds:
                    print(f"\n{res_i} - {res_j}:")
                    print(f"  Base-base H-bonds: {len(base_hbonds)}")
                    for hb in base_hbonds[:3]:
                        print(f"    {hb.donor_atom} -> {hb.acceptor_atom}: {hb.distance:.2f}Å")
                    shown += 1


if __name__ == "__main__":
    main()
