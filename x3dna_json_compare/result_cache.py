"""
Result caching system for comparison results.

Provides caching of comparison results to avoid recalculation when
source files haven't changed.
"""

import hashlib
import json
import os
from pathlib import Path
from typing import Dict, Optional, Any
from datetime import datetime
from .models import ComparisonResult


class ComparisonCache:
    """Manages cached comparison results."""
    
    def __init__(self, cache_dir: Path = Path(".comparison_cache")):
        """
        Initialize cache directory.
        
        Args:
            cache_dir: Directory to store cache files
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        self.metadata_file = self.cache_dir / "metadata.json"
        self.results_dir = self.cache_dir / "results"
        self.results_dir.mkdir(exist_ok=True)
        
        self._metadata: Optional[Dict[str, Any]] = None
    
    @property
    def metadata(self) -> Dict[str, Any]:
        """Lazy load metadata."""
        if self._metadata is None:
            if self.metadata_file.exists():
                try:
                    with open(self.metadata_file, 'r') as f:
                        self._metadata = json.load(f)
                except Exception:
                    self._metadata = {}
            else:
                self._metadata = {}
        return self._metadata
    
    def _save_metadata(self):
        """Save metadata to file."""
        self.metadata['version'] = '1.0'
        self.metadata['last_updated'] = datetime.now().isoformat()
        with open(self.metadata_file, 'w') as f:
            json.dump(self._metadata, f, indent=2)
    
    def _get_file_hash(self, file_path: Path) -> Optional[str]:
        """Compute SHA256 hash of file."""
        if not file_path.exists():
            return None
        
        sha256_hash = hashlib.sha256()
        try:
            with open(file_path, "rb") as f:
                # Read file in chunks to handle large files
                for byte_block in iter(lambda: f.read(4096), b""):
                    sha256_hash.update(byte_block)
            return sha256_hash.hexdigest()
        except Exception:
            return None
    
    def _get_file_mtime(self, file_path: Path) -> Optional[float]:
        """Get file modification time."""
        if not file_path.exists():
            return None
        try:
            return os.path.getmtime(file_path)
        except Exception:
            return None
    
    def get_cache_key(self, legacy_file: Path, modern_file: Path,
                     pdb_file: Path) -> str:
        """
        Generate cache key from file modification times (faster than hashing).
        
        Uses mtime + file size as a proxy for file identity.
        Much faster than computing SHA256 hashes for every file.
        
        Args:
            legacy_file: Path to legacy JSON file
            modern_file: Path to modern JSON file
            pdb_file: Path to PDB file
            
        Returns:
            Cache key string
        """
        # Use mtime + size for fast key generation (much faster than hashing)
        def get_file_id(file_path: Path) -> str:
            if not file_path.exists():
                return 'none'
            try:
                stat = os.stat(file_path)
                return f"{stat.st_mtime}:{stat.st_size}"
            except Exception:
                return 'none'
        
        legacy_id = get_file_id(legacy_file)
        modern_id = get_file_id(modern_file)
        pdb_id = get_file_id(pdb_file)
        
        # Combine IDs to create unique key
        combined = f"{legacy_id}:{modern_id}:{pdb_id}"
        return hashlib.md5(combined.encode()).hexdigest()
    
    def is_cache_valid(self, cache_key: str, legacy_file: Path,
                      modern_file: Path, pdb_file: Path) -> bool:
        """
        Check if cached result is still valid.
        
        Uses file modification time for fast validation instead of hashing.
        Only computes hashes when actually needed.
        
        Args:
            cache_key: Cache key for this comparison
            legacy_file: Path to legacy JSON file
            modern_file: Path to modern JSON file
            pdb_file: Path to PDB file
            
        Returns:
            True if cache is valid, False otherwise
        """
        if cache_key not in self.metadata.get('results', {}):
            return False
        
        cached_info = self.metadata['results'][cache_key]
        stored_mtimes = cached_info.get('file_mtimes', {})
        
        # Quick check using modification times (much faster than hashing)
        # Only fall back to hashes if mtimes don't match
        legacy_mtime = self._get_file_mtime(legacy_file)
        modern_mtime = self._get_file_mtime(modern_file)
        pdb_mtime = self._get_file_mtime(pdb_file)
        
        # Check mtimes first (fast)
        if stored_mtimes.get('legacy') != legacy_mtime:
            return False
        if stored_mtimes.get('modern') != modern_mtime:
            return False
        if stored_mtimes.get('pdb') != pdb_mtime:
            return False
        
        # If mtimes match, cache is valid (don't need to hash)
        return True
    
    def get_cached_result(self, pdb_id: str) -> Optional[ComparisonResult]:
        """
        Retrieve cached result if valid.
        
        Args:
            pdb_id: PDB identifier
            
        Returns:
            ComparisonResult if found and valid, None otherwise
        """
        cache_file = self.results_dir / f"{pdb_id}_full.json"
        if not cache_file.exists():
            return None
        
        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)
            
            cache_key = data.get('cache_key')
            if not cache_key or cache_key not in self.metadata.get('results', {}):
                return None
            
            cached_info = self.metadata['results'][cache_key]
            legacy_file = Path(cached_info['files']['legacy'])
            modern_file = Path(cached_info['files']['modern'])
            pdb_file = Path(cached_info['files']['pdb'])
            
            if not self.is_cache_valid(cache_key, legacy_file, modern_file, pdb_file):
                return None
            
            # Return result (simplified - full deserialization would restore all fields)
            return ComparisonResult.from_dict(data)
        except Exception:
            return None
    
    def cache_result(self, result: ComparisonResult, legacy_file: Path,
                    modern_file: Path, pdb_file: Path):
        """
        Store comparison result in cache.
        
        Args:
            result: ComparisonResult to cache
            legacy_file: Path to legacy JSON file
            modern_file: Path to modern JSON file
            pdb_file: Path to PDB file
        """
        if not result.cache_key:
            result.cache_key = self.get_cache_key(legacy_file, modern_file, pdb_file)
        
        cache_key = result.cache_key
        result_file = self.results_dir / f"{result.pdb_id}_full.json"
        
        # Store file hashes and paths in metadata
        if 'results' not in self.metadata:
            self.metadata['results'] = {}
        
        # Store mtimes for fast validation (compute hashes only when needed)
        legacy_mtime = self._get_file_mtime(legacy_file)
        modern_mtime = self._get_file_mtime(modern_file)
        pdb_mtime = self._get_file_mtime(pdb_file)
        
        self.metadata['results'][cache_key] = {
            'pdb_id': result.pdb_id,
            'timestamp': result.timestamp or datetime.now().timestamp(),
            'file_mtimes': {
                'legacy': legacy_mtime,
                'modern': modern_mtime,
                'pdb': pdb_mtime,
            },
            'files': {
                'legacy': str(legacy_file),
                'modern': str(modern_file),
                'pdb': str(pdb_file),
            }
        }
        
        # Save result data
        result_data = result.to_dict()
        result_data['cache_key'] = cache_key
        result_data['timestamp'] = result.timestamp or datetime.now().timestamp()
        
        with open(result_file, 'w') as f:
            json.dump(result_data, f, indent=2)
        
        # Update metadata
        self._save_metadata()
    
    def invalidate_cache(self, pdb_id: Optional[str] = None):
        """
        Invalidate cache for specific PDB or all.
        
        Args:
            pdb_id: PDB ID to invalidate, or None for all
        """
        if pdb_id:
            # Remove specific PDB cache
            result_file = self.results_dir / f"{pdb_id}_full.json"
            if result_file.exists():
                result_file.unlink()
            
            # Remove from metadata
            if 'results' in self.metadata:
                keys_to_remove = [
                    k for k, v in self.metadata['results'].items()
                    if v.get('pdb_id') == pdb_id
                ]
                for key in keys_to_remove:
                    del self.metadata['results'][key]
                self._save_metadata()
        else:
            # Remove all cache
            for cache_file in self.results_dir.glob("*_full.json"):
                cache_file.unlink()
            self.metadata['results'] = {}
            self._save_metadata()
    
    def clear_all(self):
        """Clear all cached results."""
        self.invalidate_cache()

