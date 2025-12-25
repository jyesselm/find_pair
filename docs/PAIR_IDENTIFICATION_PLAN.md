# Base Pair Identification Algorithm - Implementation Plan

## Goal

Design and prototype a new base pair identification algorithm that:
1. Uses the **same frame-based approach** as the modern C++ code (least-squares template fitting)
2. Extends classification to **Leontis-Westhof (LW) categories** using idealized templates
3. Combines geometric validation with hydrogen bond scoring

## Key Principle: Reuse Modern C++ Algorithm

The prototype must use the **identical algorithm** as the modern C++ code:
- **Frame Calculation**: Least-squares fitting of ring atoms to standard base templates
- **Geometric Validation**: dorg, d_v, plane_angle, dNN, overlap_area checks
- **Quality Scoring**: `quality_score = dorg + 1.5*d_v + plane_angle/180`
- **Pair Selection**: Greedy mutual best match

The extension is adding LW classification on top of validated pairs.

---

## 0. Modern C++ Algorithm Reference

Before implementing, understand the exact algorithm from the modern C++ code:

### Frame Calculation (BaseFrameCalculator)

Each nucleotide residue gets a **ReferenceFrame** via least-squares template fitting:

```
1. Load standard base template (Atomic_A.pdb, Atomic_C.pdb, etc.)
   - Ring atoms: C2, C4, C5, C6, N1, N3 (+ N7, C8, N9 for purines)

2. Match experimental atoms to template atoms
   - Minimum 3 atoms required

3. Least-squares fitting (Kabsch algorithm with quaternions)
   - Compute centroids of both point sets
   - Build covariance matrix
   - Find optimal rotation via eigenvalue decomposition
   - Compute translation

4. Store ReferenceFrame:
   - origin: base center (centroid after fitting)
   - rotation: 3x3 matrix with axes as columns
     - x_axis: column 0
     - y_axis: column 1
     - z_axis: column 2 (BASE PLANE NORMAL - critical!)
```

### Geometric Validation (BasePairValidator)

For two residues with frames (frame1, frame2):

```python
# 1. Origin distance
dorg = |frame1.origin - frame2.origin|
# Constraint: 0.0 <= dorg <= 15.0 Å

# 2. Direction vectors (dot products of axes)
dir_x = frame1.x_axis · frame2.x_axis
dir_y = frame1.y_axis · frame2.y_axis
dir_z = frame1.z_axis · frame2.z_axis
# Watson-Crick: dir_z < 0 (opposite z-axes)

# 3. Average origin and z-axis
oave = (frame1.origin + frame2.origin) / 2
if dir_z > 0:
    zave = normalize(frame1.z_axis + frame2.z_axis)
else:
    zave = normalize(frame2.z_axis - frame1.z_axis)

# 4. Vertical distance (along helix axis)
d_v = |dorg · zave|
# Constraint: 0.0 <= d_v <= 2.5 Å

# 5. Plane angle (between base planes)
plane_angle = angle_0_to_90(frame1.z_axis, frame2.z_axis)
# Constraint: 0.0 <= plane_angle <= 65.0°

# 6. Glycosidic N distance
dNN = |N1/N9_pos1 - N1/N9_pos2|
# Constraint: dNN >= 4.5 Å

# 7. Overlap area (ring clash detection)
# Constraint: overlap_area < 0.01

# 8. Quality score (lower = better)
quality_score = dorg + 1.5 * d_v + plane_angle / 180.0
```

### Pair Selection (Greedy Mutual Best Match)

```
Phase 1: Validate ALL possible pairs
  - For each (i, j) pair within 15 Å: run validation
  - Cache results

Phase 2: Greedy selection
  For each unpaired residue i:
    1. Find best partner j (lowest quality_score with valid checks)
    2. Find best partner of j
    3. If best_partner(j) == i (mutual): select pair (i, j)
    4. Mark both as matched
```

---

## 1. Directory Structure

```
prototypes/pair_identification/
    __init__.py

    # Phase 1: Load existing data + cache pairs
    frame_loader.py            # Load pre-computed frames from modern C++ JSON
    pair_cache.py              # Cache potential pairs within distance threshold
    cache_schema.py            # JSON schema definitions
    geometric_validator.py     # dorg, d_v, plane_angle, dNN checks (from frames)

    # Phase 2: Templates
    template_generator.py      # Generate idealized pair templates from DSSR data
    template_registry.py       # Load and query templates by LW class + sequence
    clustering.py              # Cluster similar pairs, compute representatives

    # Phase 3: Identification
    pair_identifier.py         # Main classification algorithm
    template_matcher.py        # Frame + geometry based template matching
    hbond_scorer.py            # H-bond quality scoring
    confidence.py              # Confidence score calculation

    # Utilities
    geometry.py                # RMSD, coordinate transforms
    atom_selection.py          # Select atoms for alignment (ring, hbond, etc.)
    
    # Testing
    tests/
        test_pair_cache.py
        test_template_generator.py
        test_template_matcher.py
        test_pair_identifier.py
        test_integration.py
    
    # Data (generated)
    data/
        potential_pairs/       # Cached potential pairs per PDB
        templates/             # Idealized pair templates
            cWW/
            tWW/
            cWH/
            ...
        templates.json         # Template registry index

scripts/
    generate_pair_cache.py     # CLI to generate potential pair cache
    generate_templates.py      # CLI to generate idealized templates
    validate_templates.py      # Validate templates against source pairs
    run_identification.py      # Run identification on test structures
```

---

## 2. Module Breakdown

### 2.0 Load Pre-computed Frames from Modern C++ JSON

**Responsibility**: Load reference frames that were already computed by the modern C++ code.

**Source Data**: `data/json/frame_calc/{PDB}.json` or `data/json_legacy/frame_calc/{PDB}.json`

**Key Classes**:

```python
@dataclass
class ReferenceFrame:
    """Reference frame for a nucleotide base (loaded from C++ JSON)."""
    origin: np.ndarray          # 3D position (base center)
    rotation: np.ndarray        # 3x3 rotation matrix (axes as columns)
    rmsd_fit: float             # RMSD from template fitting

    @property
    def x_axis(self) -> np.ndarray:
        """X-axis (column 0 of rotation matrix)."""
        return self.rotation[:, 0]

    @property
    def y_axis(self) -> np.ndarray:
        """Y-axis (column 1 of rotation matrix)."""
        return self.rotation[:, 1]

    @property
    def z_axis(self) -> np.ndarray:
        """Z-axis (column 2) - BASE PLANE NORMAL."""
        return self.rotation[:, 2]

    @classmethod
    def from_json(cls, data: Dict) -> 'ReferenceFrame':
        """Load from C++ JSON format."""
        return cls(
            origin=np.array(data["origin"]),
            rotation=np.array(data["rotation"]),
            rmsd_fit=data.get("rmsd_fit", 0.0)
        )


class FrameLoader:
    """Load pre-computed frames from modern C++ JSON output."""

    def __init__(self, json_dir: Path):
        self.json_dir = json_dir
        self.frame_calc_dir = json_dir / "frame_calc"

    def load_frames(self, pdb_id: str) -> Dict[str, ReferenceFrame]:
        """
        Load all frames for a PDB from the C++ JSON output.

        Returns: Dict mapping res_id -> ReferenceFrame
        """
        json_path = self.frame_calc_dir / f"{pdb_id}.json"
        if not json_path.exists():
            raise FileNotFoundError(f"No frame data for {pdb_id}")

        with open(json_path) as f:
            data = json.load(f)

        frames = {}
        for entry in data:
            res_id = self._build_res_id(entry)
            frames[res_id] = ReferenceFrame(
                origin=np.array(entry["origin"]),
                rotation=np.array(entry["orien"]),  # 3x3 matrix
                rmsd_fit=entry.get("rms_fit", 0.0)
            )
        return frames

    def _build_res_id(self, entry: Dict) -> str:
        """Build residue ID from JSON entry."""
        chain = entry["chain_id"]
        name = entry["residue_name"].strip()
        seq = entry["residue_seq"]
        return f"{chain}-{name}-{seq}"
```

**No frame calculation needed** - we use the frames already computed by the C++ code.

### 2.1 Phase 1: Pair Cache (`pair_cache.py`)

**Responsibility**: Find all residue pairs within distance threshold and cache their geometric/H-bond data.

**Key Classes**:

```python
@dataclass
class AtomCoords:
    """Coordinates for key atoms of a residue."""
    c1_prime: Optional[np.ndarray] = None      # C1' position
    n1_or_n9: Optional[np.ndarray] = None      # N1 (pyrimidines) or N9 (purines)
    ring_atoms: Dict[str, np.ndarray] = field(default_factory=dict)  # C2, C4, C5, C6, N1, N3, etc.
    hbond_atoms: Dict[str, np.ndarray] = field(default_factory=dict)  # N6, O6, N4, O2, O4, etc.

@dataclass
class CachedPair:
    """Cached data for a potential base pair (includes pre-computed frames)."""
    # Residue identifiers
    res1_id: str                  # e.g., "A-G-5"
    res2_id: str                  # e.g., "A-C-72"
    res1_name: str                # e.g., "G" (1-letter)
    res2_name: str                # e.g., "C"

    # PRE-COMPUTED FRAMES (from Phase 0)
    frame1: ReferenceFrame        # Frame for residue 1
    frame2: ReferenceFrame        # Frame for residue 2

    # Coordinates (for H-bond detection and template matching)
    res1_atoms: AtomCoords
    res2_atoms: AtomCoords

    # GEOMETRIC VALIDATION (computed from frames - matches C++)
    dorg: float                   # Origin distance: |frame1.origin - frame2.origin|
    d_v: float                    # Vertical distance: |dorg · zave|
    plane_angle: float            # Angle between base planes (0-90°)
    dNN: float                    # N1/N9 to N1/N9 distance
    dir_x: float                  # x-axis dot product
    dir_y: float                  # y-axis dot product
    dir_z: float                  # z-axis dot product (< 0 for WC)
    quality_score: float          # dorg + 1.5*d_v + plane_angle/180

    # VALIDATION FLAGS (matches C++ BasePairValidator)
    distance_check: bool          # 0 <= dorg <= 15
    d_v_check: bool               # 0 <= d_v <= 2.5
    plane_angle_check: bool       # 0 <= plane_angle <= 65
    dNN_check: bool               # dNN >= 4.5
    overlap_check: bool           # overlap_area < 0.01
    hbond_check: bool             # num_base_hb >= 1
    is_valid: bool                # All checks passed

    # H-bonds found between residues
    hbonds: List[Dict]            # [{donor, acceptor, distance, angle}, ...]
    num_base_hb: int              # Number of base H-bonds
    num_o2_hb: int                # Number of O2' H-bonds

class GeometricValidator:
    """Compute geometric validation metrics from frames (matches C++ BasePairValidator)."""

    # Constants from C++ validation_constants.hpp
    MAX_DORG = 15.0
    MAX_D_V = 2.5
    MAX_PLANE_ANGLE = 65.0
    MIN_DNN = 4.5
    OVERLAP_THRESHOLD = 0.01
    D_V_WEIGHT = 1.5
    PLANE_ANGLE_DIVISOR = 180.0

    def validate(
        self,
        frame1: ReferenceFrame,
        frame2: ReferenceFrame,
        n1n9_pos1: np.ndarray,
        n1n9_pos2: np.ndarray,
    ) -> Dict:
        """
        Compute all validation metrics from frames.

        Returns dict with: dorg, d_v, plane_angle, dNN, dir_x/y/z,
                          quality_score, and all check flags.
        """
        # 1. Origin distance
        dorg_vec = frame1.origin - frame2.origin
        dorg = np.linalg.norm(dorg_vec)

        # 2. Direction vectors
        dir_x = np.dot(frame1.x_axis, frame2.x_axis)
        dir_y = np.dot(frame1.y_axis, frame2.y_axis)
        dir_z = np.dot(frame1.z_axis, frame2.z_axis)

        # 3. Average z-axis
        if dir_z > 0:
            zave = frame1.z_axis + frame2.z_axis
        else:
            zave = frame2.z_axis - frame1.z_axis
        zave = zave / np.linalg.norm(zave)

        # 4. Vertical distance
        d_v = abs(np.dot(dorg_vec, zave))

        # 5. Plane angle (0-90°)
        dot = np.clip(np.dot(frame1.z_axis, frame2.z_axis), -1, 1)
        angle = np.degrees(np.arccos(abs(dot)))
        plane_angle = min(angle, 180 - angle)

        # 6. N1/N9 distance
        dNN = np.linalg.norm(n1n9_pos1 - n1n9_pos2)

        # 7. Quality score
        quality_score = dorg + self.D_V_WEIGHT * d_v + plane_angle / self.PLANE_ANGLE_DIVISOR

        # 8. Validation checks
        return {
            "dorg": dorg,
            "d_v": d_v,
            "plane_angle": plane_angle,
            "dNN": dNN,
            "dir_x": dir_x,
            "dir_y": dir_y,
            "dir_z": dir_z,
            "quality_score": quality_score,
            "distance_check": 0 <= dorg <= self.MAX_DORG,
            "d_v_check": 0 <= d_v <= self.MAX_D_V,
            "plane_angle_check": 0 <= plane_angle <= self.MAX_PLANE_ANGLE,
            "dNN_check": dNN >= self.MIN_DNN,
        }


class PairCache:
    """Cache of potential pairs for a structure (loads pre-computed frames from C++ JSON)."""

    def __init__(self, pdb_id: str, json_dir: Path):
        self.pdb_id = pdb_id
        self.json_dir = json_dir
        self.pairs: List[CachedPair] = []
        self.frames: Dict[str, ReferenceFrame] = {}  # res_id -> frame
        self.atoms: Dict[str, AtomCoords] = {}       # res_id -> coords

    def build_cache(self, max_distance: float = 15.0) -> None:
        """
        Build cache using pre-computed frames from C++ JSON.

        Algorithm:
        1. LOAD frames from existing JSON (no calculation!)
        2. LOAD atom coordinates from pdb_atoms JSON
        3. Find all pairs within max_distance using KDTree on frame origins
        4. For each pair: compute geometric validation from frames
        5. Compute H-bonds for each pair
        6. Store results
        """
        # Load pre-computed frames from C++ JSON
        frame_loader = FrameLoader(self.json_dir)
        self.frames = frame_loader.load_frames(self.pdb_id)

        # Load atom coordinates for H-bond detection
        self._load_atom_coords()

        # Build spatial index on frame origins
        res_ids = list(self.frames.keys())
        origins = np.array([self.frames[r].origin for r in res_ids])
        kdtree = KDTree(origins)

        # Find pairs within distance
        validator = GeometricValidator()
        for i, res1_id in enumerate(res_ids):
            neighbors = kdtree.query_ball_point(origins[i], r=max_distance)
            for j in neighbors:
                if j <= i:
                    continue

                res2_id = res_ids[j]
                frame1 = self.frames[res1_id]
                frame2 = self.frames[res2_id]

                # Get N1/N9 positions from loaded atoms
                n1n9_1 = self._get_n1n9_position(res1_id)
                n1n9_2 = self._get_n1n9_position(res2_id)

                # Validate using frames (same algorithm as C++)
                validation = validator.validate(frame1, frame2, n1n9_1, n1n9_2)

                # Compute H-bonds between residues
                hbonds = self._compute_hbonds(res1_id, res2_id)
                validation["hbond_check"] = len(hbonds) >= 1
                validation["is_valid"] = all([
                    validation["distance_check"],
                    validation["d_v_check"],
                    validation["plane_angle_check"],
                    validation["dNN_check"],
                    validation["hbond_check"],
                ])

                # Create cached pair
                pair = CachedPair(
                    res1_id=res1_id,
                    res2_id=res2_id,
                    frame1=frame1,
                    frame2=frame2,
                    hbonds=hbonds,
                    **validation,
                )
                self.pairs.append(pair)

    def _load_atom_coords(self):
        """Load atom coordinates from pdb_atoms JSON."""
        pdb_atoms_path = self.json_dir / "pdb_atoms" / f"{self.pdb_id}.json"
        with open(pdb_atoms_path) as f:
            data = json.load(f)
        # Parse atoms and group by residue
        ...

    def to_json(self) -> Dict:
        """Serialize to JSON (frames stored as origin + rotation)."""
        ...

    @classmethod
    def from_json(cls, data: Dict) -> 'PairCache':
        """Deserialize from JSON."""
        ...

    def save(self, path: Path):
        """Save to JSON file."""
        ...

    @classmethod
    def load(cls, path: Path) -> 'PairCache':
        """Load from JSON file."""
        ...
```

**Algorithm**:
1. Parse PDB with gemmi
2. Extract nucleotide residues with key atoms (C1', N1/N9, ring, H-bond atoms)
3. Build spatial index (KDTree) on C1' positions
4. For each residue, find all neighbors within 15 A
5. Filter: different residues, both nucleotides, min distance check
6. For each pair, compute all potential H-bonds using existing `HBondOptimizer`
7. Cache to JSON

### 2.2 Phase 2: Template Generator (`template_generator.py`)

**Responsibility**: Generate idealized pair templates from DSSR-classified pairs.

**Key Classes**:

```python
@dataclass
class PairInstance:
    """A single pair instance from DSSR data."""
    pdb_id: str
    nt1_id: str                   # e.g., "A.G1"
    nt2_id: str                   # e.g., "A.C72"
    sequence: str                 # e.g., "GC"
    lw_class: str                 # e.g., "cWW"
    n1n9_dist: float
    interbase_angle: float
    planarity: float
    hbonds_num: int
    hbonds_desc: str              # e.g., "O6-N4[2.83],N1-N3[2.88]"
    atom_coords: Dict[str, np.ndarray]  # All atom coords for this pair

@dataclass
class IdealizedTemplate:
    """An idealized base pair template."""
    sequence: str                 # e.g., "GC"
    lw_class: str                 # e.g., "cWW"
    
    # Idealized coordinates (residue 1 in standard frame, residue 2 positioned)
    res1_atoms: Dict[str, np.ndarray]
    res2_atoms: Dict[str, np.ndarray]
    
    # Expected properties
    expected_n1n9_dist: float
    expected_angle: float
    expected_planarity: float
    expected_hbonds: List[Tuple[str, str]]  # [(donor, acceptor), ...]
    
    # Statistics from source data
    sample_size: int
    rmsd_to_sources: float       # How well this template represents sources
    
    # Atoms to use for alignment
    alignment_atoms_1: List[str]  # e.g., ["C2", "C4", "C6", "N1", "N3"]
    alignment_atoms_2: List[str]  # Same for residue 2

class TemplateGenerator:
    """Generate idealized templates from DSSR data."""
    
    def __init__(self, dssr_json_dir: Path):
        self.dssr_dir = dssr_json_dir
        self.pair_instances: Dict[str, List[PairInstance]] = {}  # (seq, lw) -> instances
    
    def collect_pairs_from_dssr(self) -> None:
        """Scan all DSSR JSON files and collect classified pairs."""
        # For each DSSR JSON:
        #   - Extract pairs with LW classification
        #   - Get atom coords from original PDB
        #   - Group by (sequence, lw_class)
        ...
    
    def cluster_pairs(
        self,
        sequence: str,
        lw_class: str,
        max_clusters: int = 5,
    ) -> List[List[PairInstance]]:
        """Cluster similar pairs based on geometry."""
        # Use RMSD-based clustering
        ...
    
    def compute_idealized(
        self,
        instances: List[PairInstance],
    ) -> IdealizedTemplate:
        """Compute idealized template from multiple instances."""
        # 1. Superimpose all instances onto first
        # 2. Compute average atom positions
        # 3. Optionally: optimize H-bond geometry
        # 4. Flatten to z=0 plane (optional, for planarity)
        ...
    
    def generate_all_templates(self, output_dir: Path) -> Dict:
        """Generate templates for all sequence/LW combinations."""
        ...
```

**Algorithm**:
1. Scan 6288 DSSR JSON files in `data/json_dssr/`
2. Extract all pairs with LW classification
3. Group by (sequence, lw_class): "GC-cWW", "AU-cWW", "GA-tHS", etc.
4. For each group with >= 10 instances:
   - Cluster by geometric similarity (RMSD)
   - Select largest cluster or most representative
   - Superimpose all pairs in cluster onto canonical frame
   - Average atom positions
   - Optionally: project to z=0 plane for perfect planarity
   - Validate by checking RMSD back to source pairs
5. Save templates as PDB files (like existing `basepair-idealized/cWW/GC.pdb`)
6. Create registry JSON with metadata

### 2.3 Phase 3: Pair Identifier (`pair_identifier.py`)

**Responsibility**: Classify base pairs using template matching and H-bond scoring.

**Key Classes**:

```python
@dataclass
class MatchResult:
    """Result of matching a pair against a template."""
    template_seq: str
    template_lw: str
    rmsd: float                   # RMSD after superposition
    aligned_atoms: int            # Number of atoms used
    hbond_score: float            # H-bond quality score (0-1)
    hbonds_found: List[Dict]      # Actual H-bonds detected
    hbonds_expected: List[Tuple[str, str]]  # Expected from template
    hbonds_missing: List[Tuple[str, str]]   # Expected but not found
    hbonds_extra: List[Dict]      # Found but not expected

@dataclass
class PairClassification:
    """Final classification of a base pair."""
    res1_id: str
    res2_id: str
    sequence: str                 # e.g., "GC"
    
    # Best match
    lw_class: str                 # e.g., "cWW"
    confidence: float             # 0.0 - 1.0
    rmsd: float
    hbond_score: float
    
    # Alternative matches (if close)
    alternatives: List[MatchResult]
    
    # If no good match
    is_nonstandard: bool
    closest_template: Optional[str]  # e.g., "GC-cWW" even if poor match

class TemplateMatcher:
    """Match pairs against idealized templates."""
    
    def __init__(self, template_registry: TemplateRegistry):
        self.registry = template_registry
    
    def superimpose_pair(
        self,
        pair: CachedPair,
        template: IdealizedTemplate,
        atom_selection: str = "ring",  # "ring", "all", "hbond"
    ) -> Tuple[float, np.ndarray]:
        """
        Superimpose pair onto template and return RMSD.
        
        Returns:
            (rmsd, rotation_matrix)
        """
        ...
    
    def match(
        self,
        pair: CachedPair,
        template: IdealizedTemplate,
    ) -> MatchResult:
        """Compute full match result including H-bonds."""
        ...
    
    def find_best_match(
        self,
        pair: CachedPair,
        max_rmsd: float = 2.5,
        min_hbond_score: float = 0.3,
    ) -> Optional[MatchResult]:
        """Find best matching template for a pair."""
        ...

class HBondScorer:
    """Score H-bond quality for pair matching."""
    
    def __init__(self):
        # H-bond geometry parameters
        self.ideal_distance = 2.9  # Angstroms
        self.max_distance = 3.5
        self.ideal_angle = 160.0   # degrees (D-H...A)
        self.min_angle = 120.0
    
    def score_hbond(
        self,
        donor_pos: np.ndarray,
        acceptor_pos: np.ndarray,
        h_pos: Optional[np.ndarray] = None,
    ) -> float:
        """Score a single H-bond (0-1)."""
        ...
    
    def score_pair_hbonds(
        self,
        found_hbonds: List[Dict],
        expected_hbonds: List[Tuple[str, str]],
    ) -> Tuple[float, List, List, List]:
        """
        Score H-bonds for a pair match.
        
        Returns:
            (score, found, missing, extra)
        """
        ...

class PairIdentifier:
    """Main pair identification algorithm."""
    
    def __init__(
        self,
        template_registry: TemplateRegistry,
        hbond_scorer: HBondScorer,
    ):
        self.registry = template_registry
        self.matcher = TemplateMatcher(template_registry)
        self.hbond_scorer = hbond_scorer
    
    def classify_pair(
        self,
        pair: CachedPair,
    ) -> PairClassification:
        """
        Classify a single pair.
        
        Algorithm:
        1. Get sequence (e.g., "GC")
        2. Retrieve all templates for this sequence
        3. For each template:
           a. Superimpose and compute RMSD
           b. Score H-bonds
           c. Compute combined score
        4. Select best match above threshold
        5. If no good match: mark as non-standard
        """
        ...
    
    def classify_structure(
        self,
        pair_cache: PairCache,
        min_confidence: float = 0.5,
    ) -> List[PairClassification]:
        """Classify all pairs in a structure."""
        ...
```

**Classification Algorithm**:

```
For each potential pair (res1, res2):
    
    1. SEQUENCE LOOKUP
       sequence = get_sequence(res1, res2)  # e.g., "GC", "AU"
       templates = registry.get_templates_for_sequence(sequence)
    
    2. TEMPLATE MATCHING
       matches = []
       for template in templates:
           # Superimpose pair onto template
           rmsd = superimpose(pair, template, atoms="ring")
           
           if rmsd > MAX_RMSD_CUTOFF:  # 2.5 A
               continue
           
           # Score H-bonds
           found_hbonds = pair.hbonds
           expected_hbonds = template.expected_hbonds
           hbond_score = score_hbonds(found_hbonds, expected_hbonds)
           
           # Combined score
           score = compute_combined_score(rmsd, hbond_score)
           matches.append(MatchResult(template, rmsd, hbond_score, score))
       
       # Sort by combined score
       matches.sort(key=lambda m: m.score, reverse=True)
    
    3. CLASSIFICATION
       if matches and matches[0].score >= MIN_SCORE:
           best = matches[0]
           classification = PairClassification(
               lw_class=best.template.lw_class,
               confidence=compute_confidence(best, matches),
               rmsd=best.rmsd,
               hbond_score=best.hbond_score,
               alternatives=matches[1:3],  # Next 2 best
           )
       else:
           # No good match - non-standard pair
           classification = PairClassification(
               is_nonstandard=True,
               closest_template=matches[0].template if matches else None,
           )
    
    return classification
```

---

## 3. Data Structures

### 3.1 Potential Pairs Cache JSON Schema

```json
{
  "pdb_id": "1EHZ",
  "generated": "2024-12-24T00:00:00Z",
  "parameters": {
    "max_distance": 15.0,
    "validation_constants": {
      "max_dorg": 15.0,
      "max_d_v": 2.5,
      "max_plane_angle": 65.0,
      "min_dNN": 4.5,
      "d_v_weight": 1.5,
      "plane_angle_divisor": 180.0
    }
  },
  "frames": {
    "A-G-1": {
      "origin": [53.768, 42.0, 52.577],
      "rotation": [
        [0.123, -0.456, 0.789],
        [0.234, 0.567, -0.123],
        [-0.345, 0.678, 0.456]
      ],
      "rmsd_fit": 0.012
    },
    "A-C-72": {
      "origin": [44.234, 38.5, 51.2],
      "rotation": [...],
      "rmsd_fit": 0.015
    }
  },
  "residues": {
    "A-G-1": {
      "chain": "A",
      "name": "G",
      "seq": 1,
      "atoms": {
        "C1'": [53.768, 42.0, 52.577],
        "N9": [...],
        "C2": [...],
        "N1": [...],
        "O6": [...],
        "N2": [...]
      }
    }
  },
  "pairs": [
    {
      "res1_id": "A-G-1",
      "res2_id": "A-C-72",
      "res1_name": "G",
      "res2_name": "C",

      "geometric_validation": {
        "dorg": 10.581,
        "d_v": 0.42,
        "plane_angle": 11.674,
        "dNN": 8.851,
        "dir_x": -0.997,
        "dir_y": -0.985,
        "dir_z": -0.892,
        "quality_score": 11.275,
        "distance_check": true,
        "d_v_check": true,
        "plane_angle_check": true,
        "dNN_check": true,
        "overlap_check": true,
        "hbond_check": true,
        "is_valid": true
      },

      "hbonds": [
        {
          "donor_res": "res1",
          "donor_atom": "N1",
          "acceptor_res": "res2",
          "acceptor_atom": "N3",
          "distance": 2.88,
          "angle": 165.2
        },
        {
          "donor_res": "res1",
          "donor_atom": "N2",
          "acceptor_res": "res2",
          "acceptor_atom": "O2",
          "distance": 2.84,
          "angle": 158.7
        },
        {
          "donor_res": "res2",
          "donor_atom": "N4",
          "acceptor_res": "res1",
          "acceptor_atom": "O6",
          "distance": 2.83,
          "angle": 162.1
        }
      ],
      "num_base_hb": 3,
      "num_o2_hb": 0
    }
  ]
}
```

### 3.2 Template Registry JSON Schema

```json
{
  "version": "1.0",
  "generated": "2024-12-24T00:00:00Z",
  "templates": [
    {
      "sequence": "GC",
      "lw_class": "cWW",
      "pdb_file": "cWW/GC.pdb",
      "sample_size": 15234,

      "idealized_frames": {
        "res1": {
          "origin": [0.0, 4.5, 0.0],
          "rotation": [[1,0,0], [0,1,0], [0,0,1]]
        },
        "res2": {
          "origin": [0.0, -4.5, 0.0],
          "rotation": [[1,0,0], [0,-1,0], [0,0,-1]]
        }
      },

      "expected_geometry": {
        "dorg": 9.0,
        "d_v": 0.15,
        "plane_angle": 10.0,
        "dNN": 8.9,
        "dir_z": -1.0
      },

      "expected_hbonds": [
        {"donor": "G.N1", "acceptor": "C.N3", "distance": 2.90},
        {"donor": "G.N2", "acceptor": "C.O2", "distance": 2.85},
        {"donor": "C.N4", "acceptor": "G.O6", "distance": 2.83}
      ],

      "alignment_atoms": {
        "res1": ["C2", "C4", "C5", "C6", "N1", "N3", "N7", "C8", "N9"],
        "res2": ["C2", "C4", "C5", "C6", "N1", "N3"]
      },

      "statistics": {
        "rmsd_to_sources": 0.35,
        "dorg_std": 0.3,
        "plane_angle_std": 5.0
      }
    }
  ],
  "by_sequence": {
    "GC": ["cWW", "tWW", "cWH"],
    "AU": ["cWW", "tWW", "cHH"],
    "AG": ["tHS", "cWS", "tSS"]
  },
  "by_lw_class": {
    "cWW": ["GC", "CG", "AU", "UA", "GU", "UG"],
    "tWW": ["GC", "CG", "AU", "UA"],
    "tHS": ["AG", "GA", "AA"]
  }
}
```

### 3.3 Classification Output JSON Schema

```json
{
  "pdb_id": "1EHZ",
  "classified_pairs": [
    {
      "res1_id": "A-G-1",
      "res2_id": "A-C-72",
      "sequence": "GC",
      "lw_class": "cWW",
      "confidence": 0.95,
      "rmsd": 0.42,
      "hbond_score": 0.92,
      "hbonds_found": 3,
      "hbonds_expected": 3,
      "hbonds_missing": [],
      "is_nonstandard": false,
      "alternatives": [
        {"lw_class": "tWW", "rmsd": 1.8, "score": 0.45}
      ]
    },
    {
      "res1_id": "A-A-14",
      "res2_id": "A-G-21",
      "sequence": "AG",
      "lw_class": "tHS",
      "confidence": 0.72,
      ...
    },
    {
      "res1_id": "A-A-9",
      "res2_id": "A-A-23",
      "sequence": "AA",
      "is_nonstandard": true,
      "closest_template": "AA-cWW",
      "closest_rmsd": 3.5,
      ...
    }
  ],
  "summary": {
    "total_pairs": 76,
    "classified": 72,
    "nonstandard": 4,
    "by_lw_class": {
      "cWW": 45,
      "tWW": 8,
      "cWH": 5,
      ...
    }
  }
}
```

---

## 4. Algorithm Pseudocode

### 4.1 Cache Generation

```
function generate_pair_cache(pdb_path, output_path):
    structure = load_pdb(pdb_path)
    residues = extract_nucleotide_residues(structure)
    
    # Build spatial index
    c1_positions = [res.atoms["C1'"] for res in residues]
    kdtree = KDTree(c1_positions)
    
    pairs = []
    for i, res1 in enumerate(residues):
        # Find neighbors within 15 A
        neighbor_indices = kdtree.query_ball_point(c1_positions[i], r=15.0)
        
        for j in neighbor_indices:
            if j <= i:  # Avoid duplicates
                continue
            
            res2 = residues[j]
            
            # Compute distances
            c1_dist = distance(res1.c1_prime, res2.c1_prime)
            n1n9_dist = compute_n1n9_distance(res1, res2)
            
            if c1_dist < 4.0:  # Too close
                continue
            
            # Compute H-bonds
            hbonds = find_hbonds_between(res1, res2)
            
            # Compute interbase angle
            angle = compute_interbase_angle(res1, res2)
            
            pairs.append(CachedPair(
                res1_id=res1.res_id,
                res2_id=res2.res_id,
                c1_c1_dist=c1_dist,
                n1n9_dist=n1n9_dist,
                hbonds=hbonds,
                interbase_angle=angle,
            ))
    
    cache = PairCache(pdb_id, pairs)
    cache.save(output_path)
```

### 4.2 Template Generation

```
function generate_templates(dssr_dir, output_dir):
    # Collect all pairs from DSSR
    all_pairs = {}  # (sequence, lw_class) -> [PairInstance, ...]
    
    for dssr_json in glob(dssr_dir / "*.json"):
        data = load_json(dssr_json)
        pdb_id = get_pdb_id(dssr_json)
        
        for pair in data.get("pairs", []):
            if "LW" not in pair:
                continue
            
            lw = pair["LW"]  # e.g., "cWW"
            seq = get_sequence(pair)  # e.g., "GC"
            
            # Load atom coordinates from PDB
            coords = load_pair_coords(pdb_id, pair["nt1"], pair["nt2"])
            
            instance = PairInstance(
                pdb_id=pdb_id,
                nt1_id=pair["nt1"],
                nt2_id=pair["nt2"],
                sequence=seq,
                lw_class=lw,
                n1n9_dist=pair["N1N9_dist"],
                interbase_angle=pair["interBase_angle"],
                hbonds_desc=pair["hbonds_desc"],
                atom_coords=coords,
            )
            
            key = (seq, lw)
            if key not in all_pairs:
                all_pairs[key] = []
            all_pairs[key].append(instance)
    
    # Generate templates for each (sequence, lw) combination
    templates = []
    for (seq, lw), instances in all_pairs.items():
        if len(instances) < 10:  # Need enough samples
            continue
        
        # Cluster by geometry
        clusters = cluster_by_rmsd(instances, max_clusters=3)
        main_cluster = max(clusters, key=len)
        
        # Compute idealized template
        template = compute_average_structure(main_cluster)
        
        # Validate
        rmsd_to_sources = validate_template(template, main_cluster)
        
        template.sample_size = len(main_cluster)
        template.rmsd_to_sources = rmsd_to_sources
        templates.append(template)
        
        # Save PDB
        save_template_pdb(template, output_dir / lw / f"{seq}.pdb")
    
    # Save registry
    save_registry(templates, output_dir / "templates.json")
```

### 4.3 Pair Classification

```
function classify_pair(pair, registry):
    sequence = pair.sequence  # e.g., "GC"
    
    # Get all templates for this sequence
    templates = registry.get_templates(sequence)
    
    if not templates:
        # Try reverse sequence
        templates = registry.get_templates(reverse(sequence))
    
    if not templates:
        return NonStandardClassification(pair)
    
    # Match against each template
    matches = []
    for template in templates:
        # Superimpose pair onto template
        rmsd, transform = superimpose(
            pair_atoms=pair.get_ring_atoms(),
            template_atoms=template.get_ring_atoms(),
        )
        
        if rmsd > 2.5:  # Cutoff
            continue
        
        # Score H-bonds
        hbond_result = score_hbonds(
            found=pair.hbonds,
            expected=template.expected_hbonds,
        )
        
        # Combined score
        # RMSD contribution: 0-1 (lower RMSD = higher score)
        rmsd_score = max(0, 1 - rmsd / 2.5)
        
        # H-bond contribution: 0-1
        hbond_score = hbond_result.score
        
        # Weight: RMSD 40%, H-bonds 60%
        combined = 0.4 * rmsd_score + 0.6 * hbond_score
        
        matches.append(MatchResult(
            template=template,
            rmsd=rmsd,
            hbond_score=hbond_score,
            combined_score=combined,
        ))
    
    # Sort by combined score
    matches.sort(key=lambda m: m.combined_score, reverse=True)
    
    if not matches or matches[0].combined_score < 0.5:
        return NonStandardClassification(
            pair,
            closest=matches[0] if matches else None,
        )
    
    # Compute confidence
    best = matches[0]
    if len(matches) > 1:
        # Confidence based on gap to next best
        gap = best.combined_score - matches[1].combined_score
        confidence = min(1.0, 0.5 + gap)
    else:
        confidence = best.combined_score
    
    return PairClassification(
        lw_class=best.template.lw_class,
        confidence=confidence,
        rmsd=best.rmsd,
        hbond_score=best.hbond_score,
        alternatives=matches[1:3],
    )
```

---

## 5. File Formats

### 5.1 Idealized Template PDB Format

Based on existing files in `basepair-idealized/`:

```
REMARK   1 Idealized base pair: G-C-cWW
REMARK   2 LW classification: cWW
REMARK   3 Sequence: GC
REMARK   4 Sample size: 15234
REMARK   5 All atoms in z=0 plane (perfectly planar)
REMARK   6 H-bonds (expected):
REMARK   7   N1-N3: d=2.90 ang=165/160
REMARK   7   N2-O2: d=2.85 ang=162/158
REMARK   7   N4-O6: d=2.83 ang=160/165
ATOM      1  C1'   G A   1      -2.479   5.346   0.000  1.00  0.00            C
ATOM      2  N9    G A   1      -1.291   4.498   0.000  1.00  0.00            N
... (all ring and H-bond atoms for residue 1)
ATOM     12  C1'   C B   2      -2.399  -5.343   0.000  1.00  0.00            C
ATOM     13  N1    C B   2      -1.222  -4.462   0.000  1.00  0.00            N
... (all ring and H-bond atoms for residue 2)
END
```

**Key points**:
- Residue 1: chain A, position 1
- Residue 2: chain B, position 2
- All z-coordinates = 0 (perfectly planar)
- Includes only atoms needed for alignment and H-bond verification

### 5.2 Atoms to Include in Templates

**Ring atoms** (for RMSD alignment):
- Purines (A, G): C2, C4, C5, C6, C8, N1, N3, N7, N9
- Pyrimidines (C, U, T): C2, C4, C5, C6, N1, N3

**H-bond atoms** (for validation):
- Donors: N1 (G), N2 (G), N3 (U/T), N4 (C), N6 (A)
- Acceptors: O2 (C/U/T), O4 (U/T), O6 (G), N1 (A), N3 (A/C), N7 (A/G)

**Anchor atom**:
- C1' (for distance calculations)

---

## 6. Validation Strategy

### 6.1 Template Validation

**Goal**: Ensure idealized templates accurately represent source pairs.

**Metrics**:
1. **RMSD to sources**: Average RMSD of template to all source pairs in cluster
   - Target: < 0.5 A for main cluster members
2. **H-bond preservation**: Template should have expected H-bond geometry
   - Target: All expected H-bonds have distance < 3.5 A, angle > 120 degrees
3. **Cross-validation**: Train on 80% of pairs, test on 20%
   - Target: Test set RMSD similar to training set

**Validation script**:
```bash
python scripts/validate_templates.py \
    --templates prototypes/pair_identification/data/templates/ \
    --dssr-dir data/json_dssr/ \
    --output validation_report.json
```

### 6.2 Classification Validation

**Goal**: Validate classification accuracy against DSSR ground truth.

**Test sets**:
1. **100-PDB test set**: Existing test set with known DSSR classifications
2. **Cross-validation**: Split DSSR data 80/20, train templates on 80%, test on 20%

**Metrics**:
1. **Classification accuracy**: % of pairs classified correctly (same LW class as DSSR)
2. **Precision/Recall per LW class**: Especially for less common classes
3. **Non-standard detection**: Are truly unusual pairs flagged correctly?

**Validation script**:
```bash
python scripts/validate_classification.py \
    --test-set data/test_pdbs_100.txt \
    --output classification_report.json
```

---

## 7. Implementation Steps

### Phase 1: Load Frames + Build Pair Cache (Days 1-2)

**Step 1.1: Implement `frame_loader.py`**
- [ ] Define `ReferenceFrame` dataclass with origin, rotation, rmsd_fit
- [ ] Implement x_axis, y_axis, z_axis properties
- [ ] Implement `FrameLoader.load_frames()` to read from `data/json/frame_calc/{PDB}.json`
- [ ] Handle residue ID mapping (chain-name-seq format)
- Test: Load frames for 1EHZ, verify 76 frames loaded

**Step 1.2: Implement `geometric_validator.py`**
- [ ] Compute dorg, d_v, plane_angle, dNN from loaded frames
- [ ] Compute dir_x, dir_y, dir_z
- [ ] Compute quality_score = dorg + 1.5*d_v + plane_angle/180
- [ ] Apply validation thresholds (matches C++ constants)
- Test: Compare validation results to C++ base_pair JSON for 1EHZ

**Step 1.3: Implement `pair_cache.py`**
- [ ] Load frames using FrameLoader
- [ ] Load atom coordinates from pdb_atoms JSON
- [ ] Build KDTree on frame origins for spatial search
- [ ] Apply GeometricValidator to each pair within 15 Å
- [ ] Implement H-bond detection using existing `HBondOptimizer`
- [ ] Implement `to_json()` and `from_json()`
- Test: `test_pair_cache.py` with 1EHZ, verify geometric values match C++

**Step 1.4: Create cache generation script**
- [ ] Implement `scripts/generate_pair_cache.py` with parallel processing
- [ ] Run on 100-PDB test set
- Test: All 100 PDBs produce valid cache files

### Phase 2: Template Generation (Days 3-5)

**Step 2.1: Implement `clustering.py`**
- [ ] Implement RMSD-based clustering (superimpose on ring atoms)
- [ ] Implement cluster representative selection
- [ ] Use frames from DSSR JSON (already has frame data per pair)
- Test: Clustering produces reasonable groups for cWW pairs

**Step 2.2: Implement `template_generator.py`**
- [ ] Extract pairs from DSSR JSON with LW classification (`data/json_dssr/*.json`)
- [ ] Each DSSR pair already has frame data - use it directly!
- [ ] Group by (sequence, LW class)
- [ ] For each group with >= 10 instances:
  - Cluster by geometry similarity (dorg, plane_angle, dir_z)
  - Select main cluster
  - Compute average frame origin and rotation
  - Store idealized frame for the pair
- [ ] Compute expected geometry from averaged frames (dorg, d_v, plane_angle, dir_z)
- [ ] Extract expected H-bonds from DSSR `hbonds_desc` field
- Test: Generate templates for top 5 LW classes, verify frames are correct

**Step 2.3: Validate idealized templates against source pairs**
- [ ] For each template, compute geometry difference to all source pairs
- [ ] Verify expected H-bonds match actual H-bonds in sources
- [ ] Compare expected geometry (dorg, plane_angle) to source statistics
- Target: Geometry within 1 std of sources, H-bond preservation > 95%

**Step 2.4: Implement `template_registry.py`**
- [ ] Load templates with idealized frames and expected geometry
- [ ] Implement sequence/LW lookup
- [ ] Implement geometry-based matching
- Test: Registry correctly indexes all templates

**Step 2.5: Generate full template set**
- [ ] Run on all 6288 DSSR JSONs
- [ ] Validate all templates against sources
- [ ] Generate summary statistics
- Target: Templates for all 16 sequences x 18 LW classes = 288 possible combinations
- Report coverage: which (sequence, LW) combinations have templates

### Phase 3: Pair Identification (Days 6-8)

**Step 3.1: Implement `template_matcher.py`** (Frame-based matching)
- [ ] Implement frame superposition:
  - Transform pair frames to template canonical frame
  - Compare pair geometry (dorg, d_v, plane_angle) to template expected geometry
- [ ] Implement RMSD calculation on ring atoms after superposition
- [ ] Implement `match()` returning geometry match + RMSD
- Test: cWW pairs match cWW templates best (low RMSD, geometry match)

**Step 3.2: Implement `hbond_scorer.py`**
- [ ] Score H-bond presence (expected vs found)
- [ ] Score H-bond geometry (distance, angle)
- [ ] Penalize missing expected H-bonds
- [ ] Penalize unexpected H-bonds (optional)
- Test: Known good pairs score high, bad pairs score low

**Step 3.3: Implement `pair_identifier.py`**
- [ ] For each cached pair:
  1. Get sequence (e.g., "GC")
  2. Get all templates for sequence
  3. For each template:
     - Compute frame-based geometry match
     - Compute RMSD after superposition
     - Score H-bonds
     - Combined score = 0.3*geometry + 0.3*RMSD + 0.4*hbonds
  4. Select best match above threshold
  5. Compute confidence from gap to second-best
- [ ] Handle non-standard pairs (no good template match)
- Test: 1EHZ classifications match DSSR

**Step 3.4: Validate against DSSR ground truth**
- [ ] Run on 100-PDB test set
- [ ] Compare LW classifications to DSSR
- [ ] Compute accuracy per LW class
- Target: > 95% accuracy for cWW, > 80% overall

### Phase 4: Validation & Integration (Days 9-10)

**Step 4.1: Template validation**
- [ ] Run validation on all templates
- [ ] Fix any templates with poor RMSD
- Target: Mean RMSD < 0.5 A for all templates

**Step 4.2: Classification validation**
- [ ] Run on 100-PDB test set
- [ ] Compare to DSSR ground truth
- [ ] Analyze failures
- Target: > 90% accuracy for common LW classes

**Step 4.3: Performance optimization**
- [ ] Profile classification speed
- [ ] Add caching where needed
- Target: < 1 second per structure for typical RNA

---

## 8. Expected Challenges and Solutions

### Challenge 1: Template Coverage

**Problem**: Some sequence/LW combinations may have few examples in DSSR.

**Solution**:
- Set minimum sample size (10 instances) for template generation
- For rare combinations, use broader clustering or borrow from similar pairs
- Report confidence based on template quality

### Challenge 2: Modified Bases

**Problem**: Modified bases (5MC, PSU, etc.) may not match standard templates.

**Solution**:
- Map modified bases to parent base for template matching
- Use parent base ring atoms for superposition
- Note modification in output but classify based on geometry

### Challenge 3: Poor Geometry Structures

**Problem**: Low-resolution structures may have distorted base pairs.

**Solution**:
- Use robust RMSD threshold (2.5 A)
- Weight H-bonds more heavily than RMSD for poor structures
- Flag pairs with high RMSD but good H-bonds as "confident but distorted"

### Challenge 4: Symmetric Pairs

**Problem**: Some pairs (e.g., AA-cWW) may be classified with either orientation.

**Solution**:
- Try both orientations during matching
- Report the better match
- Store orientation information in output

### Challenge 5: Edge Cases

**Problem**: Some pairs may match multiple templates equally well.

**Solution**:
- Report all close matches as alternatives
- Lower confidence when ambiguous
- Consider using additional features (planarity, N-N distance) as tiebreakers

---

## 9. Dependencies

**Python packages**:
- `gemmi` - PDB parsing
- `numpy` - Numerical operations
- `scipy` - Spatial algorithms (KDTree, clustering)
- `scikit-learn` - Clustering algorithms (optional)
- `click` - CLI interface
- `pytest` - Testing

**Existing code to reuse**:
- `prototypes/hbond_optimizer/geometry.py` - H slot/LP slot prediction
- `prototypes/hbond_optimizer/optimizer.py` - H-bond detection
- `x3dna_json_compare/res_id_utils.py` - Residue ID utilities

---

## 10. Validation Against DSSR

Use the existing comparison script (`scripts/compare_legacy_dssr_pairs.py`) to validate:

### Frame Comparison
```bash
# Verify our frames match DSSR frames for matching pairs
python scripts/compare_legacy_dssr_pairs.py --frames --max-pdbs 500
```

Expected results (from existing comparison):
- 99.9% of matching pairs have frame origins within 0.01 Å
- Mean origin difference: 0.005 Å

### Pair Matching Baseline
```bash
# Overall pair matching between legacy and DSSR
python scripts/compare_legacy_dssr_pairs.py --summary
```

Baseline (4,120 PDBs):
- Legacy match rate: 77.8% of legacy pairs found in DSSR
- DSSR match rate: 85.3% of DSSR pairs found in legacy

### LW Classification Validation
After implementation, add to comparison script:
```bash
python scripts/compare_legacy_dssr_pairs.py --lw-classification --pdb 1EHZ
```

Compare our LW classification to DSSR's LW field.

---

## 11. Success Criteria

1. **Frame Loading + Validation (Phase 1)**:
   - Successfully load frames from C++ JSON for all test PDBs
   - Geometric validation (dorg, d_v, plane_angle) computed from loaded frames
   - Validation results match C++ base_pair JSON output

2. **Template Generation (Phase 2)**:
   - Templates for >= 90% of DSSR LW classes with >= 10 instances
   - Expected geometry (dorg, plane_angle) within 1 std of source statistics
   - Idealized templates preserve expected H-bonds (> 95%)

3. **Classification Accuracy (Phase 3)**:
   - >= 95% accuracy for cWW pairs (most common)
   - >= 85% accuracy for all LW classes with >= 100 instances
   - >= 80% overall accuracy
   - Geometry-based matching works well for classification

4. **Performance**:
   - Cache generation: < 5 seconds per PDB (loading frames + finding pairs)
   - Classification: < 1 second per PDB (using cache)

5. **Non-standard Detection**:
   - Correctly flag pairs that don't match any template
   - Report closest template with geometry metrics
   - Report which validation checks failed (dorg, d_v, plane_angle, hbonds)
