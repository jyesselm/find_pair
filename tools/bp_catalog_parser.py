#!/usr/bin/env python3
"""
Parser for bp_database.txt basepair catalog file.

Parses the Leontis-Westhof classified base pair exemplar database into
structured data classes for extraction and alignment.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass
class SymmetryOperator:
    """
    Crystallographic symmetry operator.

    Format in file: "n_xyz" where:
    - n = symmetry operation number (1-48 typically)
    - x,y,z = translations in a,b,c directions (5 = 0, 4 = -1, 6 = +1, etc.)

    Example: "1_555" = identity (no transformation)
             "4_555" = apply symmetry operation 4, no translation
             "2_755" = apply operation 2, translate +2 in a direction
    """

    operation_id: int  # 1-48 typically
    t1: int  # translation in a direction (0 = no translation)
    t2: int  # translation in b direction
    t3: int  # translation in c direction

    @classmethod
    def from_string(cls, s: str) -> "SymmetryOperator":
        """Parse symmetry operator from string like '1_555' or '4_755'."""
        s = s.strip()
        if "_" not in s:
            return cls(1, 0, 0, 0)  # Default to identity

        parts = s.split("_")
        op_id = int(parts[0])

        # Parse translations: 5 = 0, 4 = -1, 6 = +1, etc.
        trans_str = parts[1] if len(parts) > 1 else "555"
        t1 = int(trans_str[0]) - 5 if len(trans_str) > 0 else 0
        t2 = int(trans_str[1]) - 5 if len(trans_str) > 1 else 0
        t3 = int(trans_str[2]) - 5 if len(trans_str) > 2 else 0

        return cls(op_id, t1, t2, t3)

    @property
    def is_identity(self) -> bool:
        """Check if this is the identity operator (no transformation)."""
        return self.operation_id == 1 and self.t1 == 0 and self.t2 == 0 and self.t3 == 0

    def __str__(self) -> str:
        return f"{self.operation_id}_{self.t1+5}{self.t2+5}{self.t3+5}"


@dataclass
class ExemplarEntry:
    """A single base pair exemplar entry from the database."""

    sequence: str  # "GC", "Cc" (lowercase = modified base)
    iso_group: str  # "1.1", "[5.3]", "i1.2", "9.1/9.2"
    count: int  # Number of occurrences in database
    pdb_id: str  # "3R1E", "Modeled"
    resolution: Optional[float]  # None for "NA" or "Modeled"
    model_num: int  # 1, 2, etc.
    chain1: str  # "A", "0" (0 typically means first/only chain)
    chain2: str
    resnum1: int
    resnum2: int
    symop1: SymmetryOperator
    symop2: SymmetryOperator
    lw_class: str  # "cWW", "tHS", etc.

    @property
    def is_modeled(self) -> bool:
        """Check if this is a modeled (not experimental) structure."""
        return self.pdb_id.lower() == "modeled"

    @property
    def needs_symmetry(self) -> bool:
        """Check if either residue requires symmetry transformation."""
        return not self.symop1.is_identity or not self.symop2.is_identity

    @property
    def base1(self) -> str:
        """Get first base letter (uppercase)."""
        return self.sequence[0].upper() if self.sequence else ""

    @property
    def base2(self) -> str:
        """Get second base letter (uppercase)."""
        return self.sequence[1].upper() if len(self.sequence) > 1 else ""

    @property
    def bp_key(self) -> str:
        """Get base pair key like 'G-C-cWW'."""
        return f"{self.base1}-{self.base2}-{self.lw_class}"

    def __str__(self) -> str:
        return f"{self.sequence}-{self.lw_class} ({self.pdb_id}:{self.chain1}.{self.resnum1}-{self.chain2}.{self.resnum2})"


class BpDatabaseParser:
    """Parser for bp_database.txt file."""

    # Known LW classifications
    KNOWN_LW_CLASSES = {
        "cWW",
        "tWW",
        "cWH",
        "tWH",
        "cWS",
        "tWS",
        "cHW",
        "tHW",
        "cHH",
        "tHH",
        "cHS",
        "tHS",
        "cSW",
        "tSW",
        "cSH",
        "tSH",
        "cSS",
        "tSS",
    }

    def parse(self, path: Path) -> Dict[str, List[ExemplarEntry]]:
        """
        Parse bp_database.txt into dict of LW class -> entries.

        Returns:
            Dictionary mapping LW class names to lists of ExemplarEntry objects.
        """
        result: Dict[str, List[ExemplarEntry]] = {}
        current_lw: Optional[str] = None

        with open(path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip("\n\r")

                # Skip empty lines
                if not line.strip():
                    continue

                # Check for LW class header
                stripped = line.strip()
                if stripped in self.KNOWN_LW_CLASSES:
                    current_lw = stripped
                    if current_lw not in result:
                        result[current_lw] = []
                    continue

                # Skip column header lines
                if self._is_header_line(line):
                    continue

                # Parse data line
                if current_lw and "\t" in line:
                    try:
                        entry = self._parse_data_line(line, current_lw)
                        if entry:
                            result[current_lw].append(entry)
                    except Exception as e:
                        print(f"Warning: Failed to parse line {line_num}: {e}")
                        continue

        return result

    def _is_header_line(self, line: str) -> bool:
        """Check if line is the column header."""
        return line.startswith("Sequence\t")

    def _parse_data_line(self, line: str, lw_class: str) -> Optional[ExemplarEntry]:
        """Parse a tab-separated data line."""
        parts = line.split("\t")
        if len(parts) < 9:
            return None

        sequence = parts[0].strip()
        iso_group = parts[1].strip()

        # Parse count
        try:
            count = int(parts[2].strip())
        except ValueError:
            count = 0

        pdb_id = parts[3].strip()

        # Parse resolution (may be "NA" for modeled structures)
        resolution_str = parts[4].strip()
        if resolution_str.upper() == "NA" or not resolution_str:
            resolution = None
        else:
            try:
                resolution = float(resolution_str)
            except ValueError:
                resolution = None

        # Parse model number
        try:
            model_num = int(parts[5].strip())
        except ValueError:
            model_num = 1

        # Parse chains (format: "A, B" or "0, 0")
        chains = self._parse_pair_field(parts[6])
        chain1 = chains[0] if chains else "A"
        chain2 = chains[1] if len(chains) > 1 else chain1

        # Parse residue numbers (format: "7, 2")
        resnums = self._parse_pair_field(parts[7])
        try:
            resnum1 = int(resnums[0]) if resnums else 1
            resnum2 = int(resnums[1]) if len(resnums) > 1 else 1
        except ValueError:
            resnum1 = resnum2 = 1

        # Parse symmetry operators (format: "1_555, 1_555")
        symops = self._parse_pair_field(parts[8])
        symop1 = SymmetryOperator.from_string(symops[0]) if symops else SymmetryOperator(1, 0, 0, 0)
        symop2 = SymmetryOperator.from_string(symops[1]) if len(symops) > 1 else SymmetryOperator(1, 0, 0, 0)

        return ExemplarEntry(
            sequence=sequence,
            iso_group=iso_group,
            count=count,
            pdb_id=pdb_id,
            resolution=resolution,
            model_num=model_num,
            chain1=chain1,
            chain2=chain2,
            resnum1=resnum1,
            resnum2=resnum2,
            symop1=symop1,
            symop2=symop2,
            lw_class=lw_class,
        )

    def _parse_pair_field(self, field: str) -> List[str]:
        """Parse a comma-separated pair field like 'A, B' or '7, 2'."""
        parts = [p.strip() for p in field.split(",")]
        return [p for p in parts if p]


def load_database(path: Path) -> Dict[str, List[ExemplarEntry]]:
    """Convenience function to load the database."""
    parser = BpDatabaseParser()
    return parser.parse(path)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python bp_catalog_parser.py <bp_database.txt>")
        sys.exit(1)

    db_path = Path(sys.argv[1])
    entries = load_database(db_path)

    print(f"Parsed {sum(len(v) for v in entries.values())} entries across {len(entries)} LW classes:\n")

    for lw_class in sorted(entries.keys()):
        entry_list = entries[lw_class]
        needs_sym = sum(1 for e in entry_list if e.needs_symmetry)
        modeled = sum(1 for e in entry_list if e.is_modeled)
        print(f"  {lw_class}: {len(entry_list)} entries ({needs_sym} need symmetry, {modeled} modeled)")
