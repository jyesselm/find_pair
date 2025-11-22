"""
PDB file reading utilities.

Provides efficient access to PDB file lines by line number and atom lookup.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict


class PdbFileReader:
    """Efficient reader for PDB files with line number access."""
    
    def __init__(self, pdb_file: Path):
        """
        Initialize PDB file reader.
        
        Args:
            pdb_file: Path to PDB file
        """
        self.pdb_file = Path(pdb_file)
        self._lines: Optional[List[str]] = None
        self._atom_index: Optional[Dict[Tuple[str, int, str, str], Tuple[int, str]]] = None
        
        if not self.pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    @property
    def lines(self) -> List[str]:
        """Lazy load PDB lines."""
        if self._lines is None:
            with open(self.pdb_file, 'r') as f:
                self._lines = f.readlines()
        return self._lines
    
    def get_line(self, line_number: int) -> Optional[str]:
        """
        Get specific line by number (1-based).
        
        Args:
            line_number: Line number (1-based)
            
        Returns:
            Line content (without newline) or None if out of range
        """
        if line_number < 1 or line_number > len(self.lines):
            return None
        return self.lines[line_number - 1].rstrip('\n\r')
    
    def _build_atom_index(self):
        """Build index of atoms for fast lookup."""
        if self._atom_index is not None:
            return
        
        self._atom_index = {}
        for line_num, line in enumerate(self.lines, 1):
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 16:
                chain_id = line[21] if len(line) > 21 else ' '
                try:
                    residue_seq = int(line[22:26].strip())
                except (ValueError, IndexError):
                    continue
                insertion = line[26] if len(line) > 26 else ' '
                atom_name = line[12:16] if len(line) >= 16 else ''
                
                if atom_name:
                    key = (chain_id, residue_seq, insertion, atom_name)
                    self._atom_index[key] = (line_num, line.rstrip('\n\r'))
    
    def find_atom_line(self, chain_id: str, residue_seq: int,
                      insertion: str, atom_name: str) -> Optional[Tuple[int, str]]:
        """
        Find PDB line for specific atom.
        
        Args:
            chain_id: Chain identifier
            residue_seq: Residue sequence number
            insertion: Insertion code (space if none)
            atom_name: Atom name (4 characters with spaces)
            
        Returns:
            Tuple of (line_number, line_content) or None if not found
        """
        self._build_atom_index()
        key = (chain_id, residue_seq, insertion, atom_name)
        return self._atom_index.get(key)
    
    def get_residue_atoms(self, chain_id: str, residue_seq: int,
                         insertion: str = ' ') -> List[Tuple[str, int, str]]:
        """
        Get all atoms for a residue.
        
        Args:
            chain_id: Chain identifier
            residue_seq: Residue sequence number
            insertion: Insertion code (space if none)
            
        Returns:
            List of (atom_name, line_number, line_content) tuples
        """
        self._build_atom_index()
        result = []
        for (ch, seq, ins, atom), (line_num, line) in self._atom_index.items():
            if ch == chain_id and seq == residue_seq and ins == insertion:
                result.append((atom, line_num, line))
        return sorted(result, key=lambda x: x[0])  # Sort by atom name
    
    def get_atom_lines_by_names(self, chain_id: str, residue_seq: int,
                                insertion: str, atom_names: List[str]) -> Dict[str, Tuple[int, str]]:
        """
        Get PDB lines for multiple atoms.
        
        Args:
            chain_id: Chain identifier
            residue_seq: Residue sequence number
            insertion: Insertion code (space if none)
            atom_names: List of atom names to find
            
        Returns:
            Dictionary mapping atom_name -> (line_number, line_content)
        """
        result = {}
        for atom_name in atom_names:
            atom_line = self.find_atom_line(chain_id, residue_seq, insertion, atom_name)
            if atom_line:
                result[atom_name] = atom_line
        return result
    
    def reload(self):
        """Reload PDB file (invalidate cache)."""
        self._lines = None
        self._atom_index = None


def get_pdb_line(pdb_file: Path, line_number: int) -> Optional[str]:
    """
    Convenience function to get a single line from a PDB file.
    
    Args:
        pdb_file: Path to PDB file
        line_number: Line number (1-based)
        
    Returns:
        Line content or None
    """
    try:
        reader = PdbFileReader(pdb_file)
        return reader.get_line(line_number)
    except Exception:
        return None

