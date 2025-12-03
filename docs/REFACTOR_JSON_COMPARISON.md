# Refactoring Plan: x3dna_json_compare Module

## Current Problems

### 1. Code Duplication
- Comparison logic duplicated in test scripts (`test_ls_fitting.py`, `test_frames_batch.py`, `test_atoms_batch.py`)
- JSON loading/normalization logic repeated everywhere
- File finding logic scattered across scripts
- Each test script reimplements comparison workflow

### 2. Inconsistent Data Models
- Some comparisons return `FrameComparison` objects
- Some return dictionaries
- Some return tuples `(bool, Dict)`
- No unified result model
- Hard to extend or compose comparisons

### 3. No Regeneration Support
- Can't regenerate specific JSON types (e.g., just `ls_fitting` for a PDB)
- No way to check what's missing and regenerate only that
- Regeneration logic scattered in test scripts

### 4. Threading/Parallelization Issues
- `ParallelExecutor` exists but not consistently used
- Test scripts implement their own threading (`ProcessPoolExecutor`)
- No unified parallel processing API
- No progress tracking or result aggregation

### 5. No Batching Infrastructure
- Each test script implements its own batching
- No shared batching utilities
- No way to resume interrupted batches
- No batch result aggregation

### 6. File Finding Scattered
- `json_file_finder.py` exists but not used consistently
- Test scripts implement their own file finding
- No unified API for finding legacy/modern JSON files

## Target Architecture

### Core Principles
1. **Single Source of Truth**: All comparison logic in `x3dna_json_compare`
2. **Unified Data Models**: Consistent result types across all comparisons
3. **Composable**: Easy to combine multiple comparison types
4. **Regeneratable**: Can regenerate specific JSON types on demand
5. **Parallelizable**: Built-in threading/process support
6. **Batchable**: Built-in batching with resume support

## Proposed Structure

```
x3dna_json_compare/
├── __init__.py                 # Public API
├── core/                       # Core comparison engine
│   ├── __init__.py
│   ├── comparator.py          # Main JsonComparator (refactored)
│   ├── batch_processor.py     # Batch processing with resume
│   ├── parallel_executor.py    # Unified parallel execution
│   └── result_aggregator.py   # Aggregate batch results
├── models/                     # Data models
│   ├── __init__.py
│   ├── json_models.py         # JSON data models (AtomRecord, FrameRecord, etc.)
│   ├── comparison_models.py  # Comparison result models
│   └── batch_models.py        # Batch processing models
├── comparisons/                # Comparison implementations
│   ├── __init__.py
│   ├── atom_comparison.py     # Atom comparison
│   ├── frame_comparison.py    # Frame comparison
│   ├── step_comparison.py     # Step parameter comparison
│   ├── pair_comparison.py     # Pair validation comparison
│   ├── hbond_comparison.py    # Hydrogen bond comparison
│   └── residue_indices_comparison.py
├── generators/                 # JSON generation orchestration
│   ├── __init__.py
│   ├── legacy_generator.py    # Legacy JSON generation
│   ├── modern_generator.py   # Modern JSON generation
│   ├── regeneration_manager.py # Regenerate specific types
│   └── file_finder.py         # Unified file finding
├── utils/                      # Utilities
│   ├── __init__.py
│   ├── json_loader.py         # Unified JSON loading/normalization
│   ├── pdb_utils.py           # PDB file utilities
│   └── cache.py               # Result caching
└── cli/                        # Command-line interface
    ├── __init__.py
    ├── compare.py             # Compare command
    ├── batch.py               # Batch processing command
    └── regenerate.py          # Regeneration command
```

## Detailed Design

### 1. Unified Data Models (`models/`)

#### JSON Data Models (`json_models.py`)
```python
@dataclass
class JsonRecord:
    """Base class for all JSON records."""
    type: str  # "pdb_atoms", "base_frame_calc", "ls_fitting", etc.
    data: Dict[str, Any]
    
@dataclass
class AtomRecord(JsonRecord):
    """Atom record with normalized fields."""
    atom_idx: int
    atom_name: str
    residue_name: str
    chain_id: str
    residue_seq: int
    insertion: str
    xyz: Tuple[float, float, float]
    # ... other fields
    
@dataclass
class FrameRecord(JsonRecord):
    """Frame calculation record."""
    residue_idx: int
    legacy_residue_idx: Optional[int]
    residue_name: str
    chain_id: str
    residue_seq: int
    insertion: str
    rms_fit: float
    # ... other fields
```

#### Comparison Result Models (`comparison_models.py`)
```python
@dataclass
class ComparisonResult:
    """Unified comparison result."""
    pdb_id: str
    status: str  # "match", "mismatch", "error", "missing"
    atom_comparison: Optional[AtomComparison] = None
    frame_comparison: Optional[FrameComparison] = None
    step_comparison: Optional[StepComparison] = None
    # ... other comparisons
    errors: List[str] = field(default_factory=list)
    timestamp: float = field(default_factory=time.time)
    
    def is_match(self) -> bool:
        """Check if all comparisons match."""
        return self.status == "match"
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization."""
        ...
    
    def to_json(self) -> str:
        """Convert to JSON string."""
        ...
```

#### Batch Models (`batch_models.py`)
```python
@dataclass
class BatchJob:
    """Single job in a batch."""
    pdb_id: str
    pdb_file: Path
    legacy_files: Dict[str, Path]  # type -> file path
    modern_files: Dict[str, Path]
    status: str = "pending"  # "pending", "running", "completed", "failed"
    result: Optional[ComparisonResult] = None
    error: Optional[str] = None
    
@dataclass
class BatchResult:
    """Result of batch processing."""
    batch_id: str
    total_jobs: int
    completed_jobs: int
    failed_jobs: int
    results: List[ComparisonResult] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    start_time: float = field(default_factory=time.time)
    end_time: Optional[float] = None
    
    def success_rate(self) -> float:
        """Calculate success rate."""
        if self.total_jobs == 0:
            return 0.0
        return (self.completed_jobs - self.failed_jobs) / self.total_jobs
```

### 2. Unified JSON Loading (`utils/json_loader.py`)

```python
class JsonLoader:
    """Unified JSON loading and normalization."""
    
    @staticmethod
    def load_file(file_path: Path) -> List[Dict]:
        """Load JSON file and normalize to list of records."""
        ...
    
    @staticmethod
    def normalize_records(data: Any, record_type: str) -> List[Dict]:
        """Normalize various JSON formats to list of records."""
        # Handle:
        # - List of records
        # - Dict with "calculations" key
        # - Dict with "type" key
        # - Single record dict
        ...
    
    @staticmethod
    def extract_by_type(records: List[Dict], record_type: str) -> List[Dict]:
        """Extract records of specific type."""
        ...
```

### 3. Unified File Finding (`generators/file_finder.py`)

```python
class JsonFileFinder:
    """Unified file finding for legacy and modern JSON."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.legacy_base = project_root / "data" / "json_legacy"
        self.modern_base = project_root / "data" / "json_modern"
    
    def find_legacy_file(self, pdb_id: str, record_type: str) -> Optional[Path]:
        """Find legacy JSON file for a record type."""
        ...
    
    def find_modern_file(self, pdb_id: str, record_type: str) -> Optional[Path]:
        """Find modern JSON file for a record type."""
        ...
    
    def find_all_files(self, pdb_id: str, record_types: List[str]) -> Dict[str, Dict[str, Optional[Path]]]:
        """Find all files for a PDB and record types.
        
        Returns:
            {
                "legacy": {"pdb_atoms": Path, "ls_fitting": Path, ...},
                "modern": {"pdb_atoms": Path, "ls_fitting": Path, ...}
            }
        """
        ...
    
    def check_missing(self, pdb_id: str, record_types: List[str]) -> Dict[str, List[str]]:
        """Check which files are missing.
        
        Returns:
            {
                "legacy": ["ls_fitting", "base_frame_calc"],
                "modern": ["pdb_atoms"]
            }
        """
        ...
```

### 4. Regeneration Manager (`generators/regeneration_manager.py`)

```python
class RegenerationManager:
    """Manage regeneration of specific JSON types."""
    
    def __init__(self, project_root: Path, legacy_exe: Optional[Path] = None,
                 modern_exe: Optional[Path] = None):
        self.project_root = project_root
        self.legacy_exe = legacy_exe
        self.modern_exe = modern_exe
        self.file_finder = JsonFileFinder(project_root)
    
    def check_what_needs_regeneration(self, pdb_id: str, 
                                     record_types: List[str]) -> Dict[str, List[str]]:
        """Check what needs to be regenerated."""
        missing = self.file_finder.check_missing(pdb_id, record_types)
        return missing
    
    def regenerate_legacy(self, pdb_id: str, record_types: List[str]) -> bool:
        """Regenerate legacy JSON for specific types."""
        # Legacy generates all types together, so regenerate all
        ...
    
    def regenerate_modern(self, pdb_id: str, record_types: List[str]) -> bool:
        """Regenerate modern JSON for specific types.
        
        Uses --stage flag to generate only requested types.
        """
        pdb_file = self.project_root / "data" / "pdb" / f"{pdb_id}.pdb"
        output_dir = self.project_root / "data" / "json_modern"
        
        # Map record types to stages
        stage_map = {
            "pdb_atoms": "atoms",
            "residue_indices": "residue_indices",
            "ls_fitting": "ls_fitting",
            "base_frame_calc": "frames",
            "frame_calc": "frames",
        }
        
        stages = [stage_map.get(rt, rt) for rt in record_types]
        unique_stages = list(set(stages))
        
        for stage in unique_stages:
            cmd = [str(self.modern_exe), str(pdb_file), str(output_dir), f"--stage={stage}"]
            subprocess.run(cmd, ...)
    
    def regenerate_all_missing(self, pdb_id: str, record_types: List[str]) -> bool:
        """Regenerate all missing files."""
        missing = self.check_what_needs_regeneration(pdb_id, record_types)
        
        if missing.get("legacy"):
            self.regenerate_legacy(pdb_id, record_types)
        if missing.get("modern"):
            self.regenerate_modern(pdb_id, missing["modern"])
        
        return True
```

### 5. Unified Parallel Executor (`core/parallel_executor.py`)

```python
class ParallelExecutor:
    """Unified parallel execution with progress tracking."""
    
    def __init__(self, max_workers: int = None, use_processes: bool = True):
        self.max_workers = max_workers or os.cpu_count()
        self.use_processes = use_processes
        self.executor_class = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    
    def execute_batch(self, jobs: List[Callable], 
                     progress_callback: Optional[Callable] = None) -> List[Any]:
        """Execute batch of jobs in parallel.
        
        Args:
            jobs: List of callables to execute
            progress_callback: Optional callback(completed, total, result)
        
        Returns:
            List of results in same order as jobs
        """
        results = [None] * len(jobs)
        
        with self.executor_class(max_workers=self.max_workers) as executor:
            future_to_index = {executor.submit(job): i for i, job in enumerate(jobs)}
            
            completed = 0
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                try:
                    results[index] = future.result()
                    completed += 1
                    if progress_callback:
                        progress_callback(completed, len(jobs), results[index])
                except Exception as e:
                    results[index] = {"error": str(e)}
                    completed += 1
                    if progress_callback:
                        progress_callback(completed, len(jobs), {"error": str(e)})
        
        return results
```

### 6. Batch Processor (`core/batch_processor.py`)

```python
class BatchProcessor:
    """Process batches with resume support."""
    
    def __init__(self, batch_file: Path, result_file: Path):
        self.batch_file = batch_file
        self.result_file = result_file
        self.jobs: List[BatchJob] = []
        self.results: BatchResult = None
    
    def load_batch(self) -> List[BatchJob]:
        """Load batch from file."""
        ...
    
    def save_batch(self, jobs: List[BatchJob]):
        """Save batch to file."""
        ...
    
    def load_results(self) -> BatchResult:
        """Load results from file."""
        ...
    
    def save_results(self, results: BatchResult):
        """Save results incrementally."""
        ...
    
    def process_batch(self, comparator: JsonComparator,
                    executor: ParallelExecutor,
                    regenerate_missing: bool = False) -> BatchResult:
        """Process entire batch with resume support."""
        # Load existing results
        if self.result_file.exists():
            self.results = self.load_results()
        else:
            self.results = BatchResult(batch_id=self.batch_file.stem, ...)
        
        # Filter to pending jobs
        pending_jobs = [j for j in self.jobs if j.status == "pending"]
        
        # Process in batches
        for job in pending_jobs:
            try:
                # Check/regenerate missing files
                if regenerate_missing:
                    regen_manager = RegenerationManager(...)
                    regen_manager.regenerate_all_missing(job.pdb_id, ...)
                
                # Compare
                result = comparator.compare_pdb(job.pdb_id, ...)
                job.result = result
                job.status = "completed" if result.is_match() else "mismatch"
                
                # Save incrementally
                self.results.results.append(result)
                self.save_results(self.results)
                self.save_batch(self.jobs)
                
            except Exception as e:
                job.status = "failed"
                job.error = str(e)
                self.results.errors.append(f"{job.pdb_id}: {e}")
                self.save_results(self.results)
                self.save_batch(self.jobs)
        
        return self.results
```

### 7. Refactored JsonComparator (`core/comparator.py`)

```python
class JsonComparator:
    """Unified comparison engine."""
    
    def __init__(self, project_root: Path, file_finder: Optional[JsonFileFinder] = None):
        self.project_root = project_root
        self.file_finder = file_finder or JsonFileFinder(project_root)
        self.json_loader = JsonLoader()
    
    def compare_pdb(self, pdb_id: str, 
                   record_types: List[str] = None,
                   auto_regenerate: bool = False) -> ComparisonResult:
        """Compare all record types for a PDB.
        
        Args:
            pdb_id: PDB identifier
            record_types: List of record types to compare (None = all)
            auto_regenerate: If True, regenerate missing files automatically
        
        Returns:
            ComparisonResult with all comparisons
        """
        if record_types is None:
            record_types = ["pdb_atoms", "residue_indices", "ls_fitting", 
                          "base_frame_calc", "frame_calc", "bpstep_params", ...]
        
        # Check for missing files
        missing = self.file_finder.check_missing(pdb_id, record_types)
        
        if auto_regenerate and (missing.get("legacy") or missing.get("modern")):
            regen_manager = RegenerationManager(self.project_root)
            regen_manager.regenerate_all_missing(pdb_id, record_types)
            # Re-check after regeneration
            missing = self.file_finder.check_missing(pdb_id, record_types)
        
        # Build result
        result = ComparisonResult(pdb_id=pdb_id, status="match")
        
        # Compare each type
        if "pdb_atoms" in record_types:
            result.atom_comparison = self._compare_atoms(pdb_id)
        
        if "ls_fitting" in record_types or "base_frame_calc" in record_types:
            result.frame_comparison = self._compare_frames(pdb_id, record_types)
        
        # ... other comparisons
        
        # Determine overall status
        if result.atom_comparison and not result.atom_comparison.is_match():
            result.status = "mismatch"
        if result.frame_comparison and not result.frame_comparison.is_match():
            result.status = "mismatch"
        
        return result
    
    def _compare_atoms(self, pdb_id: str) -> AtomComparison:
        """Compare atoms for a PDB."""
        legacy_file = self.file_finder.find_legacy_file(pdb_id, "pdb_atoms")
        modern_file = self.file_finder.find_modern_file(pdb_id, "pdb_atoms")
        
        if not legacy_file or not modern_file:
            return AtomComparison()  # Empty comparison
        
        legacy_records = self.json_loader.load_file(legacy_file)
        modern_records = self.json_loader.load_file(modern_file)
        
        return compare_atoms(legacy_records, modern_records, ...)
    
    def _compare_frames(self, pdb_id: str, record_types: List[str]) -> FrameComparison:
        """Compare frames for a PDB."""
        frame_types = [rt for rt in record_types if rt in ["ls_fitting", "base_frame_calc", "frame_calc"]]
        
        legacy_records = []
        modern_records = []
        
        for frame_type in frame_types:
            legacy_file = self.file_finder.find_legacy_file(pdb_id, frame_type)
            modern_file = self.file_finder.find_modern_file(pdb_id, frame_type)
            
            if legacy_file:
                records = self.json_loader.load_file(legacy_file)
                # Add type field
                for r in records:
                    r["type"] = frame_type
                legacy_records.extend(records)
            
            if modern_file:
                records = self.json_loader.load_file(modern_file)
                modern_records.extend(records)
        
        pdb_file = self.project_root / "data" / "pdb" / f"{pdb_id}.pdb"
        legacy_atoms_file = self.file_finder.find_legacy_file(pdb_id, "pdb_atoms")
        legacy_atoms = []
        if legacy_atoms_file:
            legacy_atoms = self.json_loader.load_file(legacy_atoms_file)
        
        return compare_frames(legacy_records, modern_records, pdb_file, 
                            legacy_atoms=legacy_atoms)
```

## Migration Plan

### Phase 1: Create New Structure (Non-Breaking)
1. Create new directory structure
2. Move existing code to new locations
3. Create unified models
4. Create unified loaders
5. **Keep old API working** (backward compatibility)

### Phase 2: Implement New Features
1. Implement `JsonFileFinder`
2. Implement `RegenerationManager`
3. Implement `BatchProcessor`
4. Refactor `JsonComparator` to use new components

### Phase 3: Update Test Scripts
1. Refactor `test_ls_fitting.py` to use new API
2. Refactor `test_frames_batch.py` to use new API
3. Refactor `test_atoms_batch.py` to use new API
4. Remove duplicate comparison logic

### Phase 4: CLI Tools
1. Create `compare` command
2. Create `batch` command
3. Create `regenerate` command

### Phase 5: Cleanup
1. Remove old duplicate code
2. Update documentation
3. Deprecate old API (with warnings)

## Benefits

1. **No More Duplication**: All comparison logic in one place
2. **Consistent API**: Same interface for all comparison types
3. **Easy to Extend**: Add new comparison types easily
4. **Regeneration Support**: Can regenerate specific types on demand
5. **Built-in Parallelization**: Threading/process support built-in
6. **Batch Processing**: Resume support, progress tracking
7. **Better Testing**: Can test comparison logic independently
8. **Maintainable**: Clear structure, single responsibility

## Example Usage

```python
from x3dna_json_compare import JsonComparator, BatchProcessor, RegenerationManager

# Single PDB comparison
comparator = JsonComparator(project_root=Path("."))
result = comparator.compare_pdb("100D", record_types=["ls_fitting", "base_frame_calc"])
print(f"Status: {result.status}")
print(f"Frames matched: {result.frame_comparison.matched}/{result.frame_comparison.total_legacy}")

# Batch processing
batch_processor = BatchProcessor(
    batch_file=Path("batch.json"),
    result_file=Path("results.json")
)
results = batch_processor.process_batch(
    comparator=comparator,
    executor=ParallelExecutor(max_workers=10),
    regenerate_missing=True
)
print(f"Success rate: {results.success_rate():.2%}")

# Regeneration
regen_manager = RegenerationManager(project_root=Path("."))
regen_manager.regenerate_modern("100D", ["ls_fitting", "base_frame_calc"])
```

## File Renaming (Optional)

The user mentioned we can rename the module. Options:
- Keep `x3dna_json_compare` (current, descriptive)
- Rename to `json_compare` (shorter, less specific)
- Rename to `comparison` (very short, might conflict)

**Recommendation**: Keep `x3dna_json_compare` - it's descriptive and already established.

