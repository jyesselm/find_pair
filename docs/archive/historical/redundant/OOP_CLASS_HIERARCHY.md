# OOP Class Hierarchy - Visual Overview

## Core Domain Model

```
┌─────────────────────────────────────────────────────────────┐
│                        Structure ✅                        │
│  - pdb_id: string                                           │
│  - chains: vector<Chain>                                    │
│                                                              │
│  + add_chain()                                              │
│  + num_residues()                                            │
│  + nucleotides()                                             │
│  + to_json_legacy()                                          │
│                                                              │
│  Note: Base pairs are NOT part of Structure (derived data)  │
└─────────────────────────────────────────────────────────────┘
                            │
                            │ contains
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                         Chain ✅                            │
│  - chain_id: char                                            │
│  - residues: vector<Residue>                                 │
│                                                              │
│  + add_residue()                                            │
│  + sequence()                                                │
│  + nucleotides()                                             │
│  + to_json()                                                  │
└─────────────────────────────────────────────────────────────┘
                            │
                            │ contains
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                        Residue ✅                           │
│  - name: string                                             │
│  - seq_num: int                                             │
│  - chain_id: char                                            │
│  - type: ResidueType                                         │
│  - atoms: vector<Atom>                                       │
│  - reference_frame: optional<ReferenceFrame>               │
│                                                              │
│  + add_atom()                                                │
│  + find_atom()                                               │
│  + ring_atoms()                                              │
│  + is_nucleotide()                                          │
│  + one_letter_code()                                         │
│  + to_json_legacy()                                          │
└─────────────────────────────────────────────────────────────┘
                            │
                            │ contains
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                          Atom ✅                            │
│  - name: string                                             │
│  - position: Vector3D                                       │
│  - residue_name: string                                      │
│  - chain_id: char                                            │
│  - residue_seq: int                                          │
│  - record_type: char                                         │
│                                                              │
│  + distance_to()                                             │
│  + is_ring_atom()                                           │
│  + is_hydrogen_bond_donor()                                 │
│  + is_hydrogen_bond_acceptor()                              │
│  + to_json_legacy()                                          │
└─────────────────────────────────────────────────────────────┘
```

## Base Pair & Reference Frame

```
┌─────────────────────────────────────────────────────────────┐
│                       BasePair ✅                           │
│  - residue_idx1: size_t                                      │
│  - residue_idx2: size_t                                     │
│  - type: BasePairType                                        │
│  - bp_type: string                                           │
│  - frame1: optional<ReferenceFrame>                        │
│  - frame2: optional<ReferenceFrame>                         │
│  - hbonds: vector<hydrogen_bond>                            │
│                                                              │
│  + origin_distance()                                         │
│  + plane_angle()                                             │
│  + n_n_distance()                                            │
│  + direction_dot_product()                                  │
│  + to_json_legacy()                                          │
│                                                              │
│  Note: Standalone class, not part of Structure              │
└─────────────────────────────────────────────────────────────┘
                            │
                            │ uses
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                    ReferenceFrame ✅                         │
│  - rotation: Matrix3D                                       │
│  - origin: Vector3D                                         │
│                                                              │
│  + x_axis()                                                  │
│  + y_axis()                                                  │
│  + z_axis()                                                  │
│  + direction_dot_product()                                  │
│  + to_json_legacy()                                          │
└─────────────────────────────────────────────────────────────┘
```

## Geometry Classes ✅

```
┌─────────────────────────────────────────────────────────────┐
│                        Vector3D ✅                          │
│  - x, y, z: double                                            │
│                                                              │
│  + length()                                                  │
│  + normalized()                                              │
│  + dot()                                                     │
│  + cross()                                                   │
│  + operator+(), operator-(), operator*()                  │
│  + distance_to()                                             │
│  + to_json()                                                 │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                        Matrix3D ✅                          │
│  - data: array<double, 9>  (row-major)                     │
│                                                              │
│  + operator*()                                               │
│  + transpose()                                               │
│  + inverse()                                                 │
│  + determinant()                                             │
│  + at(), column(), row()                                     │
│  + to_json_legacy()                                          │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                   LeastSquaresFitter ✅                       │
│                                                              │
│  + fit(points1, points2) -> FitResult                      │
│                                                              │
│  FitResult:                                                  │
│    - rotation: Matrix3D                                      │
│    - translation: Vector3D                                    │
│    - rms: double                                             │
│                                                              │
│  Algorithm: Quaternion-based with Jacobi eigenvalue         │
│  decomposition (matches original C code)                     │
└─────────────────────────────────────────────────────────────┘
```
<｜tool▁calls▁begin｜><｜tool▁call▁begin｜>
run_terminal_cmd

## Algorithm Classes

```
┌─────────────────────────────────────────────────────────────┐
│                   BaseFrameCalculator                        │
│  - template_path: path                                       │
│  - fitter: unique_ptr<LeastSquaresFitter>                    │
│                                                              │
│  + calculate_frame(residue) -> ReferenceFrame                │
│  + calculate_all_frames(structure)                           │
│  + calculate_with_metrics() -> FrameCalculationResult        │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                     BasePairFinder                          │
│  - strategy: PairFindingStrategy                             │
│  - hb_validator: unique_ptr<HydrogenBondValidator>          │
│                                                              │
│  + find_pairs(structure) -> vector<BasePair>                │
│  + validate_pair(res1, res2) -> bool                         │
│  + validate_with_details() -> ValidationResult               │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                 HydrogenBondValidator                        │
│  - max_distance: double                                      │
│  - min_angle: double                                         │
│                                                              │
│  + find_hydrogen_bonds(res1, res2) -> vector<HydrogenBond>  │
│  + has_sufficient_hbonds() -> bool                           │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                   ParameterCalculator                        │
│                                                              │
│  + calculate_step_parameters(pair1, pair2)                  │
│    -> BasePairStepParameters                                 │
│  + calculate_helical_parameters(pair1, pair2)               │
│    -> HelicalParameters                                     │
│  + calculate_midstep_frame(frame1, frame2)                 │
│    -> ReferenceFrame                                         │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                      HelixDetector                          │
│  - helix_break_distance: double                             │
│                                                              │
│  + detect_helices(structure) -> vector<Helix>                │
│  + reorder_base_pairs(structure)                             │
│  + analyze_context(pair, all_pairs) -> BasePairContext      │
└─────────────────────────────────────────────────────────────┘
```

## Protocol Classes

```
┌─────────────────────────────────────────────────────────────┐
│                      ProtocolBase                           │
│  - config: ConfigManager*                                    │
│                                                              │
│  + execute(structure)  (pure virtual)                       │
│  + set_config_manager()                                     │
└─────────────────────────────────────────────────────────────┘
                            ▲
                            │
            ┌───────────────┴───────────────┐
            │                               │
┌───────────────────────────┐  ┌───────────────────────────┐
│    FindPairProtocol       │  │    AnalyzeProtocol        │
│  - single_strand: bool    │  │  - calculate_torsions     │
│  - find_all_pairs: bool   │  │  - simple_parameters       │
│  - divide_helices: bool   │  │  - circular: bool          │
│                           │  │                           │
│  - frame_calculator       │  │  - frame_calculator        │
│  - pair_finder           │  │  - param_calculator        │
│  - helix_detector         │  │                           │
│                           │  │  - step_params            │
│  - base_pairs             │  │  - helical_params         │
│  - helices                │  │                           │
│                           │  │  + step_parameters()      │
│  + base_pairs()           │  │  + helical_parameters()    │
│  + helices()              │  │                           │
└───────────────────────────┘  └───────────────────────────┘
```

## I/O Classes

```
┌─────────────────────────────────────────────────────────────┐
│                        PdbParser                            │
│  - include_hetatm: bool                                     │
│  - include_waters: bool                                      │
│                                                              │
│  + parse_file(path) -> Structure                            │
│  + parse_stream(stream) -> Structure                        │
│  + parse_string(string) -> Structure                        │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                        PdbWriter                            │
│                                                              │
│  + write_file(structure, path)                              │
│  + write_stream(structure, stream)                          │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                       JsonWriter                            │
│                                                              │
│  + write_structure(structure, path)                          │
│  + write_base_pairs(pairs, path)                            │
│  + write_parameters(params, path)                            │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                       JsonReader                             │
│                                                              │
│  + read_structure(path) -> Structure                        │
│  + read_base_pairs(path) -> vector<BasePair>                │
└─────────────────────────────────────────────────────────────┘
```

## Configuration

```
┌─────────────────────────────────────────────────────────────┐
│                      ConfigManager                          │
│  (Singleton)                                                │
│                                                              │
│  - thresholds: ParameterThresholds                         │
│  - x3dna_home: path                                         │
│  - include_hetatm: bool                                      │
│  - include_waters: bool                                      │
│                                                              │
│  + instance() -> ConfigManager&                             │
│  + load_from_file(path)                                      │
│  + load_from_json(json)                                      │
│  + set_defaults()                                            │
│  + thresholds() -> ParameterThresholds&                     │
└─────────────────────────────────────────────────────────────┘
                            │
                            │ uses
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                  ParameterThresholds                        │
│  - min_dNN: double = 4.5                                     │
│  - max_dNN: double = 1e18                                    │
│  - max_dorg: double = 15.0                                   │
│  - max_plane_angle: double = 65.0                            │
│  - min_base_hb: int = 1                                      │
│  - hb_dist1: double = 4.0                                   │
│  - helix_break: double = 7.5                                 │
│  - max_dv: double = 2.5                                     │
└─────────────────────────────────────────────────────────────┘
```

## Parameter Structures ✅

```
┌─────────────────────────────────────────────────────────────┐
│              BasePairStepParameters ✅                      │
│  - shift: double                                             │
│  - slide: double                                             │
│  - rise: double                                              │
│  - tilt: double                                              │
│  - roll: double                                              │
│  - twist: double                                             │
│  - midstep_frame: optional<ReferenceFrame>                  │
│                                                              │
│  + as_array() -> array<double, 6>                           │
│  + from_array() -> BasePairStepParameters                    │
│  + approximately_equal(other, tolerance) -> bool            │
│  + to_json_legacy() -> json                                 │
│  + from_json_legacy() -> BasePairStepParameters              │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                   HelicalParameters ✅                       │
│  - x_displacement: double                                    │
│  - y_displacement: double                                   │
│  - rise: double                                              │
│  - inclination: double                                       │
│  - tip: double                                               │
│  - twist: double                                             │
│  - midstep_frame: optional<ReferenceFrame>                  │
│                                                              │
│  + as_array() -> array<double, 6>                           │
│  + from_array() -> HelicalParameters                        │
│  + approximately_equal(other, tolerance) -> bool            │
│  + to_json_legacy() -> json                                 │
│  + from_json_legacy() -> HelicalParameters                  │
└─────────────────────────────────────────────────────────────┘
```
<｜tool▁calls▁begin｜><｜tool▁call▁begin｜>
read_file

## Enumeration Types

```cpp
enum class ResidueType {
    ADENINE,
    CYTOSINE,
    GUANINE,
    THYMINE,
    URACIL,
    AMINO_ACID,
    OTHER
};

enum class PurinePyrimidine {
    PURINE = 1,
    PYRIMIDINE = 0,
    NOT_BASE = -1
};

enum class BasePairType {
    WATSON_CRICK,
    WOBBLE,
    HOOGSTEEN,
    REVERSE_HOOGSTEEN,
    UNKNOWN
};

enum class PairFindingStrategy {
    BEST_PAIR,      // Greedy mutual best match
    ALL_PAIRS,      // Exhaustive search
    DISTANCE_BASED  // Simple distance-based
};
```

## Data Flow Diagram

```
┌─────────────┐
│  PDB File   │
└──────┬──────┘
       │
       ▼
┌─────────────┐
│ PdbParser   │───► Structure
└─────────────┘      ├── Chains
                     ├── Residues
                     └── Atoms
                          │
                          ▼
              ┌───────────────────────┐
              │ FindPairProtocol       │
              │  ├── BaseFrameCalc     │───► ReferenceFrames
              │  ├── BasePairFinder    │───► BasePairs
              │  └── HelixDetector     │───► Helices
              └───────────────────────┘
                          │
                          ▼
              ┌───────────────────────┐
              │ AnalyzeProtocol        │
              │  ├── BaseFrameCalc     │───► Recalculate Frames
              │  └── ParamCalculator   │───► Step Parameters
              └───────────────────────┘
                          │
                          ▼
              ┌───────────────────────┐
              │ JsonWriter             │───► JSON Output
              └───────────────────────┘
```

## Key Relationships

1. **Composition**: Structure → Chain → Residue → Atom
2. **Aggregation**: Structure → BasePair (references Residues)
3. **Association**: BasePair → ReferenceFrame (uses)
4. **Dependency**: Algorithms depend on Core objects
5. **Inheritance**: Protocols inherit from ProtocolBase
6. **Singleton**: ConfigManager (single instance)

## Design Patterns Used

1. **Singleton**: ConfigManager
2. **Strategy**: PairFindingStrategy in BasePairFinder
3. **Factory**: PdbParser creates Structure objects
4. **Template Method**: ProtocolBase defines execute() interface
5. **RAII**: All resource management through constructors/destructors
6. **Builder**: Complex object construction (optional, for future)

---

*This hierarchy provides a clear, modular structure for the modernized X3DNA library.*

