#!/usr/bin/env python3
"""
Base pair exemplar extractor with symmetry handling and frame alignment.

Extracts residues from PDB files, applies crystallographic symmetry operators,
calculates base pair reference frames, and transforms to identity frame at origin.
"""

import urllib.request
import gzip
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
from scipy.spatial.transform import Rotation

try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False
    print("Warning: gemmi not available, symmetry operations will be limited")

from bp_catalog_parser import ExemplarEntry, SymmetryOperator


@dataclass
class BasePairFrame:
    """Reference frame for a base pair."""
    origin: Tuple[float, float, float]
    x_axis: Tuple[float, float, float]
    y_axis: Tuple[float, float, float]
    z_axis: Tuple[float, float, float]


@dataclass
class Atom:
    """Simple atom representation."""
    name: str
    x: float
    y: float
    z: float
    element: str = ""
    resname: str = ""
    chain: str = ""
    resnum: int = 0

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])

    def copy_with_coords(self, x: float, y: float, z: float) -> 'Atom':
        return Atom(self.name, x, y, z, self.element, self.resname, self.chain, self.resnum)


# Ring atoms for least-squares fitting (in order around the ring)
# Purines: 9-membered ring (fused 5+6)
PURINE_RING_ATOMS = ["N9", "C8", "N7", "C5", "C4", "N3", "C2", "N1", "C6"]
# Pyrimidines: 6-membered ring
PYRIMIDINE_RING_ATOMS = ["N1", "C2", "N3", "C4", "C5", "C6"]

PURINES = {"A", "G", "I"}
PYRIMIDINES = {"C", "U", "T"}


class BaseTemplates:
    """Loads and stores standard base templates."""

    def __init__(self, template_dir: Path):
        self.templates: Dict[str, Dict[str, np.ndarray]] = {}
        self._load_templates(template_dir)

    def _load_templates(self, template_dir: Path):
        """Load all Atomic_*.pdb templates."""
        for base in ["A", "G", "C", "U", "T", "I"]:
            template_path = template_dir / f"Atomic_{base}.pdb"
            if template_path.exists():
                self.templates[base] = self._parse_template(template_path)

    def _parse_template(self, path: Path) -> Dict[str, np.ndarray]:
        """Parse template PDB file into atom name -> coords dict."""
        atoms = {}
        with open(path) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms[name] = np.array([x, y, z])
        return atoms

    def get_ring_atoms(self, base: str) -> Tuple[List[str], Dict[str, np.ndarray]]:
        """Get ring atom names and template coords for a base."""
        base = base.upper()
        if base not in self.templates:
            # Try to map modified bases to parent
            base = self._get_parent_base(base)

        if base in PURINES:
            ring_names = PURINE_RING_ATOMS
        else:
            ring_names = PYRIMIDINE_RING_ATOMS

        template = self.templates.get(base, {})
        return ring_names, template

    def _get_parent_base(self, base: str) -> str:
        """Map modified base to parent base type."""
        # Common modifications
        if base in {"1MA", "MIA", "M1A", "6MA", "M6A"}:
            return "A"
        if base in {"1MG", "2MG", "7MG", "M2G", "OMG"}:
            return "G"
        if base in {"5MC", "OMC", "5CM"}:
            return "C"
        if base in {"PSU", "H2U", "5MU", "OMU"}:
            return "U"
        # Default to A for unknown
        return "A"


class BaseFrameCalculator:
    """Calculates reference frames for nucleotide bases."""

    def __init__(self, templates: BaseTemplates):
        self.templates = templates

    def calculate_frame(self, atoms: List[Atom], base_type: str) -> Optional[BasePairFrame]:
        """
        Calculate reference frame by least-squares fitting ring atoms to template.

        Returns frame with origin at ring centroid, axes from rotation matrix.
        """
        ring_names, template = self.templates.get_ring_atoms(base_type)

        # Extract matching ring atoms
        exp_coords = []
        ref_coords = []

        atom_dict = {a.name: a.coords for a in atoms}

        for name in ring_names:
            if name in atom_dict and name in template:
                exp_coords.append(atom_dict[name])
                ref_coords.append(template[name])

        if len(exp_coords) < 4:
            return None

        exp_coords = np.array(exp_coords)
        ref_coords = np.array(ref_coords)

        # Calculate centroids
        exp_centroid = np.mean(exp_coords, axis=0)
        ref_centroid = np.mean(ref_coords, axis=0)

        # Center the coordinates
        exp_centered = exp_coords - exp_centroid
        ref_centered = ref_coords - ref_centroid

        # Calculate rotation using Kabsch algorithm (SVD)
        H = ref_centered.T @ exp_centered
        U, S, Vt = np.linalg.svd(H)

        # Handle reflection case
        d = np.linalg.det(Vt.T @ U.T)
        if d < 0:
            Vt[-1, :] *= -1

        R = Vt.T @ U.T

        # Extract axes from rotation matrix (columns)
        x_axis = tuple(R[:, 0])
        y_axis = tuple(R[:, 1])
        z_axis = tuple(R[:, 2])

        return BasePairFrame(
            origin=tuple(exp_centroid),
            x_axis=x_axis,
            y_axis=y_axis,
            z_axis=z_axis
        )


class SymmetryApplier:
    """Applies crystallographic symmetry operators to coordinates."""

    def __init__(self):
        self.structure: Optional[gemmi.Structure] = None
        self.cell: Optional[gemmi.UnitCell] = None
        self.spacegroup: Optional[gemmi.SpaceGroup] = None

    def load_structure(self, pdb_path: Path) -> bool:
        """Load structure to get cell and space group info."""
        if not HAS_GEMMI:
            return False

        try:
            self.structure = gemmi.read_structure(str(pdb_path))
            self.cell = self.structure.cell
            self.spacegroup = self.structure.spacegroup_hm
            return True
        except Exception as e:
            print(f"Warning: Could not load structure for symmetry: {e}")
            return False

    def apply(self, coords: np.ndarray, symop: SymmetryOperator) -> np.ndarray:
        """
        Apply symmetry operator to coordinates.

        Args:
            coords: Nx3 array of coordinates
            symop: Symmetry operator to apply

        Returns:
            Transformed coordinates
        """
        if symop.is_identity:
            return coords.copy()

        if not HAS_GEMMI or self.cell is None:
            raise RuntimeError("Cannot apply symmetry without gemmi and cell info")

        # Get symmetry operations from space group
        try:
            sg = gemmi.SpaceGroup(self.spacegroup) if self.spacegroup else gemmi.SpaceGroup("P 1")
            operations = sg.operations().sym_ops

            if symop.operation_id > len(operations):
                print(f"Warning: Symmetry operation {symop.operation_id} not found, using identity")
                return coords.copy()

            # Get the symmetry operation (0-indexed internally)
            op = operations[symop.operation_id - 1]

            result = np.zeros_like(coords)
            for i, coord in enumerate(coords):
                # Convert to fractional coordinates
                frac = self.cell.fractionalize(gemmi.Position(*coord))

                # Apply symmetry operation
                frac_new = op.apply_to_xyz([frac.x, frac.y, frac.z])

                # Apply translation
                frac_new[0] += symop.t1
                frac_new[1] += symop.t2
                frac_new[2] += symop.t3

                # Convert back to Cartesian
                cart = self.cell.orthogonalize(gemmi.Fractional(*frac_new))
                result[i] = [cart.x, cart.y, cart.z]

            return result

        except Exception as e:
            print(f"Warning: Symmetry application failed: {e}")
            return coords.copy()


class PdbFetcher:
    """Downloads PDB files from RCSB."""

    RCSB_URL = "https://files.rcsb.org/download/{}.pdb.gz"

    def fetch(self, pdb_id: str, output_dir: Path) -> Optional[Path]:
        """
        Download PDB file if not already present.

        Returns path to PDB file or None if download failed.
        """
        output_path = output_dir / f"{pdb_id}.pdb"

        if output_path.exists():
            return output_path

        # Try to download
        url = self.RCSB_URL.format(pdb_id.upper())
        try:
            print(f"  Downloading {pdb_id}...")
            with urllib.request.urlopen(url, timeout=30) as response:
                compressed = response.read()

            # Decompress
            with gzip.open(io.BytesIO(compressed), 'rt') as f:
                content = f.read()

            # Write to file
            output_dir.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(content)

            return output_path

        except Exception as e:
            print(f"  Failed to download {pdb_id}: {e}")
            return None


class ResidueExtractor:
    """Extracts residue atoms from PDB files."""

    def __init__(self, symmetry_applier: Optional[SymmetryApplier] = None):
        self.symmetry_applier = symmetry_applier or SymmetryApplier()

    def extract(
        self,
        pdb_path: Path,
        chain_id: str,
        resnum: int,
        model_num: int = 1,
        symop: Optional[SymmetryOperator] = None
    ) -> List[Atom]:
        """
        Extract atoms for a specific residue.

        Args:
            pdb_path: Path to PDB file
            chain_id: Chain identifier (or "0" for first chain)
            resnum: Residue number
            model_num: Model number (for NMR structures)
            symop: Symmetry operator to apply

        Returns:
            List of Atom objects
        """
        atoms = []
        target_chain = chain_id if chain_id != "0" else None
        current_model = 0
        found_chain = None

        with open(pdb_path) as f:
            for line in f:
                # Track model number
                if line.startswith("MODEL"):
                    try:
                        current_model = int(line[10:14].strip())
                    except:
                        current_model += 1
                    continue

                if line.startswith("ENDMDL"):
                    if current_model == model_num:
                        break  # Done with target model
                    continue

                if not line.startswith(("ATOM", "HETATM")):
                    continue

                # For single-model files, current_model stays 0
                if current_model != 0 and current_model != model_num:
                    continue

                # Parse chain
                chain = line[21].strip()

                # Handle chain "0" - use first nucleotide chain
                if target_chain is None:
                    resname = line[17:20].strip()
                    if resname in {"A", "G", "C", "U", "T", "DA", "DG", "DC", "DT",
                                   "ADE", "GUA", "CYT", "URA", "THY"}:
                        target_chain = chain
                        found_chain = chain

                if chain != target_chain:
                    continue

                # Parse residue number
                try:
                    rnum = int(line[22:26].strip())
                except ValueError:
                    continue

                if rnum != resnum:
                    continue

                # Skip alternate conformations (take first one)
                altloc = line[16]
                if altloc not in (" ", "A"):
                    continue

                # Parse atom
                name = line[12:16].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                resname = line[17:20].strip()
                element = line[76:78].strip() if len(line) > 76 else name[0]

                atoms.append(Atom(
                    name=name,
                    x=x, y=y, z=z,
                    element=element,
                    resname=resname,
                    chain=chain,
                    resnum=rnum
                ))

        # Apply symmetry if needed
        if symop and not symop.is_identity and atoms:
            if self.symmetry_applier.structure is None:
                self.symmetry_applier.load_structure(pdb_path)

            coords = np.array([[a.x, a.y, a.z] for a in atoms])
            new_coords = self.symmetry_applier.apply(coords, symop)

            atoms = [a.copy_with_coords(new_coords[i, 0], new_coords[i, 1], new_coords[i, 2])
                     for i, a in enumerate(atoms)]

        return atoms


def transform_point(point: np.ndarray, frame: BasePairFrame) -> np.ndarray:
    """
    Transform a point from world coordinates to the base pair frame.

    Aligns the base pair so its reference frame becomes identity at origin.
    The transformation is: p' = R^T * (p - origin)
    """
    origin = np.array(frame.origin)
    R = np.array([frame.x_axis, frame.y_axis, frame.z_axis])

    translated = point - origin
    return R @ translated


def calculate_bp_frame(frame1: BasePairFrame, frame2: BasePairFrame) -> BasePairFrame:
    """
    Calculate combined base pair reference frame as average of two base frames.

    Origin is midpoint, axes are averaged and re-orthogonalized.
    """
    # Midpoint origin
    origin = tuple((np.array(frame1.origin) + np.array(frame2.origin)) / 2)

    # Average rotation matrices
    R1 = np.array([frame1.x_axis, frame1.y_axis, frame1.z_axis]).T
    R2 = np.array([frame2.x_axis, frame2.y_axis, frame2.z_axis]).T

    # Simple average and re-orthogonalize using SVD
    R_avg = (R1 + R2) / 2
    U, _, Vt = np.linalg.svd(R_avg)
    R_ortho = U @ Vt

    # Handle reflection
    if np.linalg.det(R_ortho) < 0:
        R_ortho[:, -1] *= -1

    return BasePairFrame(
        origin=origin,
        x_axis=tuple(R_ortho[:, 0]),
        y_axis=tuple(R_ortho[:, 1]),
        z_axis=tuple(R_ortho[:, 2])
    )


def transform_atoms(atoms: List[Atom], frame: BasePairFrame) -> List[Atom]:
    """Transform all atoms to align base pair frame to identity at origin."""
    result = []
    for atom in atoms:
        new_coords = transform_point(atom.coords, frame)
        result.append(atom.copy_with_coords(new_coords[0], new_coords[1], new_coords[2]))
    return result


def format_pdb_coords(x: float, y: float, z: float) -> str:
    """Format coordinates for PDB file (8.3f format)."""
    return f"{x:8.3f}{y:8.3f}{z:8.3f}"


def write_exemplar_pdb(
    output_path: Path,
    atoms1: List[Atom],
    atoms2: List[Atom],
    entry: ExemplarEntry,
    bp_frame: Optional[BasePairFrame] = None
) -> None:
    """
    Write exemplar base pair to PDB file with REMARK headers.

    Optionally transforms coordinates to align base pair frame to identity.
    """
    with open(output_path, 'w') as f:
        f.write(f"REMARK   1 Base pair exemplar: {entry.bp_key}\n")
        f.write(f"REMARK   2 Source PDB: {entry.pdb_id}\n")
        f.write(f"REMARK   3 Residue 1: {entry.chain1}.{entry.resnum1}\n")
        f.write(f"REMARK   4 Residue 2: {entry.chain2}.{entry.resnum2}\n")
        f.write(f"REMARK   5 LW classification: {entry.lw_class}\n")
        f.write(f"REMARK   6 Sequence: {entry.sequence}\n")
        f.write(f"REMARK   7 Isostericity group: {entry.iso_group}\n")
        if entry.resolution:
            f.write(f"REMARK   8 Resolution: {entry.resolution:.2f}\n")
        else:
            f.write(f"REMARK   8 Resolution: N/A\n")
        if entry.needs_symmetry:
            f.write(f"REMARK   9 Symmetry: {entry.symop1} / {entry.symop2}\n")
        f.write(f"REMARK  10 Aligned to origin with identity frame\n")

        # Transform atoms if frame provided
        if bp_frame:
            atoms1 = transform_atoms(atoms1, bp_frame)
            atoms2 = transform_atoms(atoms2, bp_frame)

        # Write atoms
        atom_num = 1
        for atom in atoms1 + atoms2:
            # PDB ATOM format
            line = f"ATOM  {atom_num:5d}  {atom.name:<3s} {atom.resname:>3s} {atom.chain:1s}{atom.resnum:4d}    "
            line += format_pdb_coords(atom.x, atom.y, atom.z)
            line += f"  1.00  0.00          {atom.element:>2s}\n"
            f.write(line)
            atom_num += 1

        f.write("END\n")


def encode_sequence_for_filename(sequence: str) -> str:
    """
    Encode sequence for case-insensitive filesystem safety.

    Lowercase letters are prefixed with underscore to distinguish them.
    E.g., "Cc" -> "C_c", "cC" -> "_cC", "GC" -> "GC"
    """
    result = []
    for char in sequence:
        if char.islower():
            result.append(f"_{char}")
        else:
            result.append(char)
    return "".join(result)


class BasePairExemplarExtractor:
    """Main orchestrator for extracting and aligning base pair exemplars."""

    def __init__(
        self,
        pdb_dir: Path,
        template_dir: Path,
        output_dir: Path,
        download_missing: bool = True
    ):
        self.pdb_dir = pdb_dir
        self.output_dir = output_dir
        self.download_missing = download_missing

        self.templates = BaseTemplates(template_dir)
        self.frame_calculator = BaseFrameCalculator(self.templates)
        self.symmetry_applier = SymmetryApplier()
        self.residue_extractor = ResidueExtractor(self.symmetry_applier)
        self.fetcher = PdbFetcher()

    def extract(self, entry: ExemplarEntry, verbose: bool = False) -> Optional[Path]:
        """
        Extract and align a single base pair exemplar.

        Returns path to output PDB file or None if extraction failed.
        """
        if entry.is_modeled:
            if verbose:
                print(f"  Skipping modeled entry: {entry}")
            return None

        # Get PDB file
        pdb_path = self.pdb_dir / f"{entry.pdb_id}.pdb"
        if not pdb_path.exists():
            if self.download_missing:
                pdb_path = self.fetcher.fetch(entry.pdb_id, self.pdb_dir)
                if pdb_path is None:
                    return None
            else:
                if verbose:
                    print(f"  PDB not found: {entry.pdb_id}")
                return None

        # Load structure for symmetry if needed
        if entry.needs_symmetry:
            self.symmetry_applier.load_structure(pdb_path)

        # Extract residue atoms
        atoms1 = self.residue_extractor.extract(
            pdb_path, entry.chain1, entry.resnum1,
            entry.model_num, entry.symop1
        )
        atoms2 = self.residue_extractor.extract(
            pdb_path, entry.chain2, entry.resnum2,
            entry.model_num, entry.symop2
        )

        if not atoms1 or not atoms2:
            if verbose:
                print(f"  Could not extract atoms: {entry}")
            return None

        # Calculate base frames
        frame1 = self.frame_calculator.calculate_frame(atoms1, entry.base1)
        frame2 = self.frame_calculator.calculate_frame(atoms2, entry.base2)

        if frame1 is None or frame2 is None:
            if verbose:
                print(f"  Could not calculate frames: {entry}")
            return None

        # Calculate combined base pair frame
        bp_frame = calculate_bp_frame(frame1, frame2)

        # Create output directory structure
        lw_dir = self.output_dir / entry.lw_class
        lw_dir.mkdir(parents=True, exist_ok=True)

        # Write output PDB (encode sequence for case-insensitive filesystems)
        safe_sequence = encode_sequence_for_filename(entry.sequence)
        output_file = lw_dir / f"{safe_sequence}.pdb"
        write_exemplar_pdb(output_file, atoms1, atoms2, entry, bp_frame)

        return output_file


if __name__ == "__main__":
    # Quick test
    from bp_catalog_parser import load_database

    project_root = Path(__file__).parent.parent
    db_path = project_root / "data" / "bp_database.txt"
    pdb_dir = project_root / "data" / "pdb"
    template_dir = project_root / "resources" / "templates"
    output_dir = project_root / "basepair-catalog-exemplars"

    entries = load_database(db_path)

    extractor = BasePairExemplarExtractor(pdb_dir, template_dir, output_dir)

    # Test with first cWW entry
    if "cWW" in entries and entries["cWW"]:
        entry = entries["cWW"][0]
        print(f"Testing with: {entry}")
        result = extractor.extract(entry, verbose=True)
        if result:
            print(f"Created: {result}")
