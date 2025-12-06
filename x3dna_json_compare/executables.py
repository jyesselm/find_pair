"""
Executable management utilities.

Handles finding and running legacy/modern executables for JSON generation.
"""

import subprocess
import shutil
from pathlib import Path
from typing import Tuple, Optional


def find_executables(project_root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find legacy and modern executables.
    
    Args:
        project_root: Project root directory
        
    Returns:
        Tuple of (legacy_exe, modern_exe), either may be None if not found
    """
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    modern_exe = project_root / "build" / "generate_modern_json"
    
    if not legacy_exe.exists():
        legacy_exe = None
    if not modern_exe.exists():
        modern_exe = None
    
    return legacy_exe, modern_exe


def generate_modern_json(
    pdb_id: str, 
    pdb_file: Path, 
    output_dir: Path,
    modern_exe: Path, 
    stage: str = "all"
) -> Tuple[bool, str]:
    """Generate modern JSON for a specific stage.
    
    Args:
        pdb_id: PDB identifier
        pdb_file: Path to PDB file
        output_dir: Output directory for generated files
        modern_exe: Path to modern executable
        stage: Stage to generate (atoms, frames, hbonds, pairs, steps, all)
    
    Returns:
        (success, message) tuple
    """
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        
        cmd = [str(modern_exe), str(pdb_file), str(output_dir)]
        if stage != "all":
            cmd.append(f"--stage={stage}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout
        )
        
        if result.returncode != 0:
            return False, f"Generation failed: {result.stderr[:200]}"
        
        return True, "OK"
        
    except subprocess.TimeoutExpired:
        return False, "Generation timeout (exceeded 2 minutes)"
    except Exception as e:
        return False, f"Error: {str(e)}"


def generate_legacy_json(
    pdb_id: str,
    pdb_file: Path,
    output_dir: Path,
    legacy_exe: Path,
    project_root: Path,
    record_types: list
) -> Tuple[bool, str]:
    """Generate legacy JSON for specified record types.
    
    Args:
        pdb_id: PDB identifier
        pdb_file: Path to PDB file
        output_dir: Output directory for generated files
        legacy_exe: Path to legacy executable
        project_root: Project root directory
        record_types: List of record types to check (e.g., ["ls_fitting", "base_frame_calc"])
    
    Returns:
        (success, message) tuple
    """
    try:
        legacy_output = output_dir / "legacy"
        legacy_output.mkdir(parents=True, exist_ok=True)
        
        # Check if legacy JSON files already exist
        legacy_base = project_root / "data" / "json_legacy"
        existing_files = {}
        all_exist = True
        
        for record_type in record_types:
            json_file = legacy_base / record_type / f"{pdb_id}.json"
            # Check for ref_frame as alternative to frame_calc
            if record_type == "frame_calc" and not json_file.exists():
                ref_frame_file = legacy_base / "ref_frame" / f"{pdb_id}.json"
                if ref_frame_file.exists():
                    json_file = ref_frame_file
            
            existing_files[record_type] = json_file
            if not json_file.exists():
                all_exist = False
        
        if all_exist:
            # Copy existing files
            for record_type, json_file in existing_files.items():
                (legacy_output / record_type).mkdir(exist_ok=True)
                shutil.copy2(json_file, legacy_output / record_type / f"{pdb_id}.json")
            return True, "OK (existing files)"
        
        # Files don't exist - regenerate using legacy executable
        if pdb_file.is_absolute():
            try:
                pdb_file_rel = pdb_file.relative_to(project_root)
            except ValueError:
                pdb_file_rel = pdb_file
        else:
            pdb_file_rel = pdb_file
        
        # Build path for org/ directory
        if str(pdb_file_rel).startswith("data/"):
            pdb_path_for_org = "../" + str(pdb_file_rel)
        else:
            pdb_path_for_org = str(pdb_file_rel)
        
        cmd = [str(legacy_exe.resolve()), pdb_path_for_org]
        
        # Run with timeout
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                timeout=300,  # 5 minute timeout
                cwd=str(project_root / "org")
            )
            
            if result.returncode != 0:
                return False, f"Legacy generation failed (return code {result.returncode})"
        except subprocess.TimeoutExpired:
            return False, "Legacy generation timeout (exceeded 5 minutes)"
        
        # Check if files were created and copy them
        all_created = True
        for record_type, json_file in existing_files.items():
            if record_type == "frame_calc" and not json_file.exists():
                ref_frame_file = legacy_base / "ref_frame" / f"{pdb_id}.json"
                if ref_frame_file.exists():
                    json_file = ref_frame_file
            
            if json_file.exists():
                (legacy_output / record_type).mkdir(exist_ok=True)
                shutil.copy2(json_file, legacy_output / record_type / f"{pdb_id}.json")
            else:
                all_created = False
        
        if all_created:
            return True, "OK"
        else:
            missing = [rt for rt, f in existing_files.items() if not f.exists()]
            return False, f"Legacy JSON files not created: {', '.join(missing)}"
        
    except Exception as e:
        return False, f"Legacy error: {str(e)}"


def cleanup_output_dir(output_dir: Path):
    """Clean up output directory after test."""
    if output_dir.exists():
        shutil.rmtree(output_dir)

