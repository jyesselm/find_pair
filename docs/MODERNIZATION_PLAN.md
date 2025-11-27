# Modern C++ Modular Library Design Plan
## Converting X3DNA to Modern C++ with Strong OOP

## Executive Summary

This document outlines a comprehensive plan to convert the legacy X3DNA v2.4 codebase (~15,000 lines of C code) into a modern, modular C++ library with strong object-oriented design principles. The new architecture will be testable, maintainable, and extensible.

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Core Library Structure](#core-library-structure)
3. [OOP Class Hierarchy](#oop-class-hierarchy)
4. [Module Breakdown](#module-breakdown)
5. [Configuration Management](#configuration-management)
6. [PDB Parser & Structure Objects](#pdb-parser--structure-objects)
7. [Protocols & Algorithms](#protocols--algorithms)
8. [Applications Layer](#applications-layer)
9. [Testing Strategy](#testing-strategy)
10. [Migration Path](#migration-path)

---

## Architecture Overview

### Design Principles

1. **Separation of Concerns**: Clear boundaries between data structures, algorithms, I/O, and applications
2. **Dependency Inversion**: High-level modules depend on abstractions, not concrete implementations
3. **Single Responsibility**: Each class has one clear purpose
4. **Open/Closed**: Open for extension, closed for modification
5. **Interface Segregation**: Small, focused interfaces
6. **RAII**: Resource management through constructors/destructors
7. **Modern C++**: C++17/20 features, smart pointers, STL containers

### Library Structure

```
libx3dna/
├── include/x3dna/
│   ├── core/              # Core data structures
│   ├── io/                # I/O operations (PDB, JSON, etc.)
│   ├── geometry/          # Geometric calculations
│   ├── algorithms/        # Core algorithms
│   ├── protocols/         # Protocol implementations
│   └── config/            # Configuration management
├── src/                   # Implementation files
├── apps/                  # Application executables
├── tests/                 # Test suite
└── examples/              # Example usage
```

---

## Core Library Structure

### Directory Layout

```
find_pair_2/
├── include/
│   └── x3dna/
│       ├── core/              # Core domain objects
│       │   ├── Structure.hpp
│       │   ├── Chain.hpp
│       │   ├── Residue.hpp
│       │   ├── Atom.hpp
│       │   ├── BasePair.hpp
│       │   ├── ReferenceFrame.hpp
│       │   └── Parameters.hpp
│       ├── io/                 # I/O operations
│       │   ├── PdbParser.hpp
│       │   ├── PdbWriter.hpp
│       │   ├── JsonReader.hpp
│       │   └── JsonWriter.hpp
│       ├── geometry/           # Geometric utilities
│       │   ├── Vector3D.hpp
│       │   ├── Matrix3D.hpp
│       │   ├── Quaternion.hpp
│       │   ├── LeastSquares.hpp
│       │   └── GeometryUtils.hpp
│       ├── algorithms/         # Core algorithms
│       │   ├── BaseFrameCalculator.hpp
│       │   ├── BasePairFinder.hpp
│       │   ├── ParameterCalculator.hpp
│       │   ├── HelixDetector.hpp
│       │   └── HydrogenBondValidator.hpp
│       ├── protocols/          # Protocol implementations
│       │   ├── FindPairProtocol.hpp
│       │   ├── AnalyzeProtocol.hpp
│       │   └── ProtocolBase.hpp
│       └── config/             # Configuration
│           ├── ConfigManager.hpp
│           ├── ConfigReader.hpp
│           └── ParameterThresholds.hpp
├── src/                        # Implementation
│   └── [mirrors include structure]
├── apps/                       # Application executables
│   ├── find_pair_app.cpp
│   └── analyze_app.cpp
├── tests/                      # Test suite
│   ├── unit/                   # Unit tests
│   ├── integration/           # Integration tests
│   └── regression/            # Regression tests (vs JSON)
│       └── compare_json.py    # Compare with legacy JSON
└── examples/                   # Example code
    └── basic_usage.cpp
```

---

## OOP Class Hierarchy

### Core Domain Objects

#### 1. Atom
```cpp
namespace x3dna::core {

class Atom {
public:
    // Constructors
    Atom(const std::string& name, const Vector3D& position);
    Atom(const std::string& name, double x, double y, double z);
    
    // Getters
    const std::string& name() const noexcept { return name_; }
    const Vector3D& position() const noexcept { return position_; }
    double x() const noexcept { return position_.x(); }
    double y() const noexcept { return position_.y(); }
    double z() const noexcept { return position_.z(); }
    
    // Setters
    void set_position(const Vector3D& pos) { position_ = pos; }
    void set_name(const std::string& name) { name_ = name; }
    
    // Utilities
    double distance_to(const Atom& other) const;
    bool is_ring_atom() const;
    bool is_hydrogen_bond_donor() const;
    bool is_hydrogen_bond_acceptor() const;
    
    // JSON serialization
    nlohmann::json to_json() const;
    static Atom from_json(const nlohmann::json& j);
    
private:
    std::string name_;      // Atom name (e.g., " C1'", " N1 ")
    Vector3D position_;     // 3D coordinates
};

// JSON format matches: {"atom_name": " C1'", "xyz": [x, y, z]}

} // namespace x3dna::core
```

#### 2. Residue
```cpp
namespace x3dna::core {

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

class Residue {
public:
    // Constructors
    Residue(const std::string& name, int sequence_number, char chain_id);
    
    // Getters
    const std::string& name() const noexcept { return name_; }
    int sequence_number() const noexcept { return sequence_number_; }
    char chain_id() const noexcept { return chain_id_; }
    ResidueType type() const noexcept { return type_; }
    PurinePyrimidine ry_classification() const noexcept { return ry_; }
    
    // Atom management
    void add_atom(const Atom& atom);
    const std::vector<Atom>& atoms() const noexcept { return atoms_; }
    std::vector<Atom>& atoms() noexcept { return atoms_; }
    
    // Atom queries
    std::optional<Atom> find_atom(const std::string& name) const;
    std::vector<Atom> ring_atoms() const;
    std::vector<Atom> hydrogen_bond_atoms() const;
    
    // Base identification
    char one_letter_code() const;
    bool is_nucleotide() const;
    bool is_purine() const;
    bool is_pyrimidine() const;
    
    // Reference frame
    void set_reference_frame(const ReferenceFrame& frame);
    const std::optional<ReferenceFrame>& reference_frame() const { return frame_; }
    
    // JSON serialization
    nlohmann::json to_json() const;
    static Residue from_json(const nlohmann::json& j);
    
private:
    std::string name_;              // Residue name (e.g., "  A", "  C")
    int sequence_number_;            // PDB sequence number
    char chain_id_;                  // Chain identifier
    ResidueType type_;               // Residue classification
    PurinePyrimidine ry_;            // Purine/Pyrimidine classification
    std::vector<Atom> atoms_;       // Atoms in this residue
    std::optional<ReferenceFrame> frame_;  // Base reference frame
};

// JSON format matches legacy structure with atoms array

} // namespace x3dna::core
```

#### 3. Chain
```cpp
namespace x3dna::core {

class Chain {
public:
    // Constructors
    explicit Chain(char id);
    
    // Getters
    char id() const noexcept { return id_; }
    const std::vector<Residue>& residues() const noexcept { return residues_; }
    std::vector<Residue>& residues() noexcept { return residues_; }
    
    // Residue management
    void add_residue(const Residue& residue);
    void add_residue(Residue&& residue);
    
    // Queries
    size_t size() const noexcept { return residues_.size(); }
    bool empty() const noexcept { return residues_.empty(); }
    const Residue& operator[](size_t index) const { return residues_[index]; }
    Residue& operator[](size_t index) { return residues_[index]; }
    
    // Sequence operations
    std::string sequence() const;  // One-letter codes
    std::vector<Residue> nucleotides() const;  // Filter to nucleotides only
    
    // JSON serialization
    nlohmann::json to_json() const;
    static Chain from_json(const nlohmann::json& j);
    
private:
    char id_;
    std::vector<Residue> residues_;
};

} // namespace x3dna::core
```

#### 4. Structure
```cpp
namespace x3dna::core {

class Structure {
public:
    // Constructors
    Structure() = default;
    explicit Structure(const std::string& pdb_id);
    
    // Getters
    const std::string& pdb_id() const noexcept { return pdb_id_; }
    const std::vector<Chain>& chains() const noexcept { return chains_; }
    std::vector<Chain>& chains() noexcept { return chains_; }
    
    // Chain management
    void add_chain(const Chain& chain);
    void add_chain(Chain&& chain);
    Chain& get_chain(char chain_id);
    const Chain& get_chain(char chain_id) const;
    
    // Queries
    size_t num_chains() const noexcept { return chains_.size(); }
    size_t num_residues() const;
    size_t num_atoms() const;
    
    // Residue access
    std::vector<Residue*> all_residues();
    std::vector<const Residue*> all_residues() const;
    std::vector<Residue*> nucleotides();
    std::vector<const Residue*> nucleotides() const;
    
    // Base pair management (after detection)
    void set_base_pairs(const std::vector<BasePair>& pairs);
    const std::vector<BasePair>& base_pairs() const noexcept { return base_pairs_; }
    
    // JSON serialization
    nlohmann::json to_json() const;
    static Structure from_json(const nlohmann::json& j);
    
    // Legacy JSON format support (for compatibility with data/json_legacy/*.json)
    nlohmann::json to_json_legacy() const;  // Matches original JSON structure
    static Structure from_json_legacy(const nlohmann::json& j);
    
private:
    std::string pdb_id_;
    std::vector<Chain> chains_;
    std::vector<BasePair> base_pairs_;  // Detected base pairs
};

} // namespace x3dna::core
```

#### 5. ReferenceFrame
```cpp
namespace x3dna::core {

class ReferenceFrame {
public:
    // Constructors
    ReferenceFrame();
    ReferenceFrame(const Matrix3D& rotation, const Vector3D& origin);
    
    // Getters
    const Matrix3D& rotation() const noexcept { return rotation_; }
    const Vector3D& origin() const noexcept { return origin_; }
    Matrix3D& rotation() noexcept { return rotation_; }
    Vector3D& origin() noexcept { return origin_; }
    
    // Axis access
    Vector3D x_axis() const;
    Vector3D y_axis() const;
    Vector3D z_axis() const;  // Normal to base plane
    
    // Transformations
    Vector3D transform_point(const Vector3D& point) const;
    ReferenceFrame transform_frame(const ReferenceFrame& frame) const;
    
    // Utilities
    ReferenceFrame inverse() const;
    double direction_dot_product(const ReferenceFrame& other) const;  // z-axis dot product
    
    // Serialization
    std::array<double, 9> rotation_as_array() const;  // Flattened 3x3 matrix
    std::array<double, 3> origin_as_array() const;
    
    // JSON serialization (matches legacy format)
    nlohmann::json to_json() const;
    static ReferenceFrame from_json(const nlohmann::json& j);
    
    // Legacy format support (for compatibility with JSON files)
    // Legacy format: {"orien": [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]], "org": [x, y, z]}
    nlohmann::json to_json_legacy() const;  // Uses "orien" and "org" keys
    static ReferenceFrame from_json_legacy(const nlohmann::json& j);
    
private:
    Matrix3D rotation_;  // 3x3 rotation matrix
    Vector3D origin_;   // Origin point
};

} // namespace x3dna::core
```

#### 6. BasePair
```cpp
namespace x3dna::core {

enum class BasePairType {
    WATSON_CRICK,
    WOBBLE,
    HOOGSTEEN,
    REVERSE_HOOGSTEEN,
    UNKNOWN
};

class BasePair {
public:
    // Constructors
    BasePair(Residue* residue1, Residue* residue2, BasePairType type);
    BasePair(size_t residue_idx1, size_t residue_idx2, BasePairType type);
    
    // Getters
    Residue* residue1() const noexcept { return residue1_; }
    Residue* residue2() const noexcept { return residue2_; }
    size_t residue_index1() const noexcept { return residue_idx1_; }
    size_t residue_index2() const noexcept { return residue_idx2_; }
    BasePairType type() const noexcept { return type_; }
    
    // Reference frames
    void set_frame1(const ReferenceFrame& frame);
    void set_frame2(const ReferenceFrame& frame);
    const ReferenceFrame& frame1() const { return frame1_; }
    const ReferenceFrame& frame2() const { return frame2_; }
    
    // Validation
    bool is_valid() const noexcept { return is_valid_; }
    void set_valid(bool valid) { is_valid_ = valid; }
    
    // Hydrogen bonds
    void set_hydrogen_bonds(const std::vector<HydrogenBond>& hbonds);
    const std::vector<HydrogenBond>& hydrogen_bonds() const { return hbonds_; }
    size_t num_hydrogen_bonds() const { return hbonds_.size(); }
    
    // Distance/angle metrics
    double origin_distance() const;
    double plane_angle() const;
    double n_n_distance() const;
    
    // JSON serialization
    nlohmann::json to_json() const;
    static BasePair from_json(const nlohmann::json& j);
    
    // Legacy format support
    // Legacy format: {"type": "base_pair", "base_i": 1, "base_j": 72, "bp_type": "AT",
    //                 "dir_xyz": [x, y, z], "orien_i": [[...]], "orien_j": [[...]], 
    //                 "org_i": [x, y, z], "org_j": [x, y, z]}
    nlohmann::json to_json_legacy() const;
    static BasePair from_json_legacy(const nlohmann::json& j);
    
private:
    Residue* residue1_;
    Residue* residue2_;
    size_t residue_idx1_;
    size_t residue_idx2_;
    BasePairType type_;
    ReferenceFrame frame1_;
    ReferenceFrame frame2_;
    bool is_valid_;
    std::vector<HydrogenBond> hbonds_;
};

struct HydrogenBond {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance;
    double angle;
    bool is_standard;
    char type;  // '-' for standard, ' ' for non-standard
    int linkage_type;
    
    // JSON serialization
    nlohmann::json to_json() const;
    static HydrogenBond from_json(const nlohmann::json& j);
    
    // Legacy format: {"hbond_idx": 1, "donor_atom": " N1 ", "acceptor_atom": " N3 ",
    //                 "distance": 2.95, "type": "-", "linkage_type": 18}
    nlohmann::json to_json_legacy() const;
    static HydrogenBond from_json_legacy(const nlohmann::json& j);
};

} // namespace x3dna::core
```

#### 7. Parameters
```cpp
namespace x3dna::core {

struct BasePairStepParameters {
    double shift;   // x-displacement
    double slide;   // y-displacement
    double rise;    // z-displacement
    double tilt;    // rotation about x (degrees)
    double roll;    // rotation about y (degrees)
    double twist;   // rotation about z (degrees)
    
    // Comparison
    bool operator==(const BasePairStepParameters& other) const;
    bool approximately_equal(const BasePairStepParameters& other, double tolerance = 0.001) const;
    
    // Serialization
    std::array<double, 6> as_array() const;
    static BasePairStepParameters from_array(const std::array<double, 6>& arr);
    
    // JSON serialization (matches legacy format)
    nlohmann::json to_json() const;
    static BasePairStepParameters from_json(const nlohmann::json& j);
    
    // Legacy format: {"params": {"Shift": ..., "Slide": ..., ...}}
    nlohmann::json to_json_legacy() const;
    static BasePairStepParameters from_json_legacy(const nlohmann::json& j);
};

struct HelicalParameters {
    double x_displacement;
    double y_displacement;
    double rise;
    double inclination;
    double tip;
    double twist;
    
    std::array<double, 6> as_array() const;
    
    // JSON serialization (matches legacy format)
    nlohmann::json to_json() const;
    static HelicalParameters from_json(const nlohmann::json& j);
    
    // Legacy format: {"params": [x_displacement, y_displacement, ...]}
    nlohmann::json to_json_legacy() const;
    static HelicalParameters from_json_legacy(const nlohmann::json& j);
};

} // namespace x3dna::core
```

---

## Module Breakdown

### 1. Configuration Management

#### ConfigManager
```cpp
namespace x3dna::config {

class ConfigManager {
public:
    // Singleton access
    static ConfigManager& instance();
    
    // Configuration loading
    void load_from_file(const std::filesystem::path& config_path);
    void load_from_json(const nlohmann::json& json);
    void set_defaults();
    
    // Parameter thresholds
    const ParameterThresholds& thresholds() const { return thresholds_; }
    ParameterThresholds& thresholds() { return thresholds_; }
    
    // Paths
    void set_x3dna_home(const std::filesystem::path& path);
    std::filesystem::path x3dna_home() const { return x3dna_home_; }
    std::filesystem::path standard_base_path() const;
    
    // Options
    bool include_hetatm() const { return include_hetatm_; }
    void set_include_hetatm(bool value) { include_hetatm_ = value; }
    bool include_waters() const { return include_waters_; }
    void set_include_waters(bool value) { include_waters_ = value; }
    
    // Legacy mode (for exact compatibility with legacy code)
    // When enabled, breaks some OOP principles for exact matching
    bool legacy_mode() const { return legacy_mode_; }
    void set_legacy_mode(bool value) { legacy_mode_ = value; }
    
private:
    ConfigManager() = default;
    ParameterThresholds thresholds_;
    std::filesystem::path x3dna_home_;
    bool include_hetatm_ = false;
    bool include_waters_ = false;
    bool legacy_mode_ = false;  // Enable legacy compatibility mode
};

struct ParameterThresholds {
    double min_dNN = 4.5;              // Minimum N-N distance
    double max_dNN = 1e18;             // Maximum N-N distance
    double max_dorg = 15.0;            // Maximum origin distance
    double max_plane_angle = 65.0;     // Maximum plane angle (degrees)
    int min_base_hb = 1;               // Minimum hydrogen bonds
    double hb_dist1 = 4.0;              // H-bond distance threshold
    double helix_break = 7.5;          // Helix break distance
    double max_dv = 2.5;               // Maximum vertical distance
};

} // namespace x3dna::config
```

### 2. PDB Parser

#### PdbParser
```cpp
namespace x3dna::io {

class PdbParser {
public:
    // Parsing
    Structure parse_file(const std::filesystem::path& pdb_path);
    Structure parse_stream(std::istream& stream);
    Structure parse_string(const std::string& pdb_content);
    
    // Options
    void set_include_hetatm(bool value) { include_hetatm_ = value; }
    void set_include_waters(bool value) { include_waters_ = value; }
    
    // Error handling
    class ParseError : public std::runtime_error {
    public:
        ParseError(const std::string& message, size_t line_number);
        size_t line_number() const { return line_number_; }
    private:
        size_t line_number_;
    };
    
private:
    void parse_atom_line(const std::string& line, Structure& structure);
    void parse_hetatm_line(const std::string& line, Structure& structure);
    ResidueType identify_residue_type(const std::string& resname);
    bool include_hetatm_ = false;
    bool include_waters_ = false;
};

} // namespace x3dna::io
```

### 3. Algorithms

#### BaseFrameCalculator
```cpp
namespace x3dna::algorithms {

class BaseFrameCalculator {
public:
    // Frame calculation
    ReferenceFrame calculate_frame(Residue& residue);
    ReferenceFrame calculate_frame(const Residue& residue) const;
    
    // Batch calculation
    void calculate_all_frames(Structure& structure);
    
    // Quality metrics
    struct FrameCalculationResult {
        ReferenceFrame frame;
        double rms_fit;
        std::vector<std::string> matched_atoms;
        std::string template_file;
    };
    
    FrameCalculationResult calculate_with_metrics(Residue& residue);
    
    // Configuration
    void set_template_path(const std::filesystem::path& path);
    
private:
    std::filesystem::path template_path_;
    std::unique_ptr<LeastSquaresFitter> fitter_;
    
    Structure load_standard_base(ResidueType type);
    std::vector<Atom> match_ring_atoms(const Residue& residue, const Structure& standard);
    ReferenceFrame fit_frame(const std::vector<Atom>& experimental, 
                           const std::vector<Atom>& standard);
};

} // namespace x3dna::algorithms
```

#### BasePairFinder
```cpp
namespace x3dna::algorithms {

enum class PairFindingStrategy {
    BEST_PAIR,      // Greedy mutual best match
    ALL_PAIRS,      // Exhaustive search
    DISTANCE_BASED  // Simple distance-based
};

class BasePairFinder {
public:
    // Finding pairs
    std::vector<BasePair> find_pairs(Structure& structure);
    std::vector<BasePair> find_pairs(const Structure& structure) const;
    
    // Strategy
    void set_strategy(PairFindingStrategy strategy) { strategy_ = strategy; }
    
    // Validation
    bool validate_pair(const Residue& res1, const Residue& res2) const;
    
    // Results
    struct ValidationResult {
        bool is_valid;
        BasePairType type;
        double dir_z;
        double dorg;
        double plane_angle;
        double dNN;
        std::vector<HydrogenBond> hbonds;
    };
    
    ValidationResult validate_with_details(const Residue& res1, const Residue& res2) const;
    
private:
    PairFindingStrategy strategy_ = PairFindingStrategy::BEST_PAIR;
    std::unique_ptr<HydrogenBondValidator> hb_validator_;
    
    std::vector<BasePair> find_best_pairs(Structure& structure);
    std::vector<BasePair> find_all_pairs(const Structure& structure);
    BasePair find_best_partner(const Residue& residue, const Structure& structure) const;
};

} // namespace x3dna::algorithms
```

#### ParameterCalculator
```cpp
namespace x3dna::algorithms {

class ParameterCalculator {
public:
    // Step parameters
    BasePairStepParameters calculate_step_parameters(
        const BasePair& pair1, 
        const BasePair& pair2
    );
    
    BasePairStepParameters calculate_step_parameters(
        const ReferenceFrame& frame1,
        const ReferenceFrame& frame2
    );
    
    // Helical parameters
    HelicalParameters calculate_helical_parameters(
        const BasePair& pair1,
        const BasePair& pair2
    );
    
    // Batch calculation
    std::vector<BasePairStepParameters> calculate_all_step_parameters(
        const std::vector<BasePair>& pairs
    );
    
    // Midstep frame
    ReferenceFrame calculate_midstep_frame(
        const ReferenceFrame& frame1,
        const ReferenceFrame& frame2
    );
    
private:
    // Core algorithm (matches original bpstep_par)
    void bpstep_par_impl(
        const Matrix3D& r1, const Vector3D& o1,
        const Matrix3D& r2, const Vector3D& o2,
        BasePairStepParameters& params,
        ReferenceFrame& midstep_frame
    );
};

} // namespace x3dna::algorithms
```

#### HydrogenBondValidator
```cpp
namespace x3dna::algorithms {

class HydrogenBondValidator {
public:
    // Validation
    std::vector<HydrogenBond> find_hydrogen_bonds(
        const Residue& res1,
        const Residue& res2
    ) const;
    
    bool has_sufficient_hbonds(
        const Residue& res1,
        const Residue& res2,
        int min_count = 1
    ) const;
    
    // Configuration
    void set_max_distance(double distance) { max_distance_ = distance; }
    void set_min_angle(double angle) { min_angle_ = angle; }
    
private:
    double max_distance_ = 4.0;
    double min_angle_ = 120.0;  // degrees
    
    bool is_donor_acceptor_pair(const Atom& atom1, const Atom& atom2) const;
    double calculate_hbond_angle(const Atom& donor, const Atom& acceptor) const;
};

} // namespace x3dna::algorithms
```

#### HelixDetector
```cpp
namespace x3dna::algorithms {

struct Helix {
    std::vector<size_t> base_pair_indices;  // Indices into base_pairs vector
    size_t start_index;
    size_t end_index;
    bool is_circular;
};

class HelixDetector {
public:
    // Detection
    std::vector<Helix> detect_helices(const Structure& structure);
    std::vector<Helix> detect_helices(const std::vector<BasePair>& pairs);
    
    // Ordering
    void reorder_base_pairs(Structure& structure);
    void ensure_five_to_three_ordering(std::vector<BasePair>& pairs);
    
    // Analysis
    struct BasePairContext {
        std::vector<size_t> neighbors;
        bool is_helix_start;
        bool is_helix_end;
        double coplanarity;
    };
    
    BasePairContext analyze_context(const BasePair& pair, const std::vector<BasePair>& all_pairs);
    
private:
    double helix_break_distance_ = 7.5;
    
    bool are_coplanar(const BasePair& p1, const BasePair& p2) const;
    bool is_circular(const std::vector<BasePair>& pairs) const;
};

} // namespace x3dna::algorithms
```

### 4. Protocols

#### ProtocolBase
```cpp
namespace x3dna::protocols {

class ProtocolBase {
public:
    virtual ~ProtocolBase() = default;
    
    // Execution
    virtual void execute(Structure& structure) = 0;
    virtual void execute(const Structure& structure) const = 0;
    
    // Configuration
    void set_config_manager(ConfigManager& config) { config_ = &config; }
    ConfigManager& config() { return *config_; }
    const ConfigManager& config() const { return *config_; }
    
    // Error handling
    class ProtocolError : public std::runtime_error {
    public:
        ProtocolError(const std::string& message);
    };
    
protected:
    ConfigManager* config_ = nullptr;
};

} // namespace x3dna::protocols
```

#### FindPairProtocol
```cpp
namespace x3dna::protocols {

class FindPairProtocol : public ProtocolBase {
public:
    // Execution
    void execute(Structure& structure) override;
    
    // Options
    void set_single_strand_mode(bool value) { single_strand_ = value; }
    void set_find_all_pairs(bool value) { find_all_pairs_ = value; }
    void set_divide_helices(bool value) { divide_helices_ = value; }
    
    // Results
    const std::vector<BasePair>& base_pairs() const { return base_pairs_; }
    const std::vector<Helix>& helices() const { return helices_; }
    
private:
    bool single_strand_ = false;
    bool find_all_pairs_ = false;
    bool divide_helices_ = false;
    
    std::unique_ptr<BaseFrameCalculator> frame_calculator_;
    std::unique_ptr<BasePairFinder> pair_finder_;
    std::unique_ptr<HelixDetector> helix_detector_;
    
    std::vector<BasePair> base_pairs_;
    std::vector<Helix> helices_;
    
    void calculate_frames(Structure& structure);
    void find_pairs(Structure& structure);
    void detect_helices(Structure& structure);
    void reorder_pairs(Structure& structure);
};

} // namespace x3dna::protocols
```

#### AnalyzeProtocol
```cpp
namespace x3dna::protocols {

class AnalyzeProtocol : public ProtocolBase {
public:
    // Execution
    void execute(Structure& structure) override;
    
    // Options
    void set_calculate_torsions(bool value) { calculate_torsions_ = value; }
    void set_simple_parameters(bool value) { simple_parameters_ = value; }
    void set_circular_structure(bool value) { circular_ = value; }
    
    // Results
    const std::vector<BasePairStepParameters>& step_parameters() const { return step_params_; }
    const std::vector<HelicalParameters>& helical_parameters() const { return helical_params_; }
    
private:
    bool calculate_torsions_ = false;
    bool simple_parameters_ = false;
    bool circular_ = false;
    
    std::unique_ptr<BaseFrameCalculator> frame_calculator_;
    std::unique_ptr<ParameterCalculator> param_calculator_;
    
    std::vector<BasePairStepParameters> step_params_;
    std::vector<HelicalParameters> helical_params_;
    
    void recalculate_frames(Structure& structure);
    void calculate_parameters(Structure& structure);
};

} // namespace x3dna::protocols
```

### 5. Geometry Utilities

#### Vector3D, Matrix3D, etc.
```cpp
namespace x3dna::geometry {

class Vector3D {
public:
    Vector3D() : x_(0), y_(0), z_(0) {}
    Vector3D(double x, double y, double z) : x_(x), y_(y), z_(z) {}
    
    double x() const { return x_; }
    double y() const { return y_; }
    double z() const { return z_; }
    
    double length() const;
    Vector3D normalized() const;
    double dot(const Vector3D& other) const;
    Vector3D cross(const Vector3D& other) const;
    
    Vector3D operator+(const Vector3D& other) const;
    Vector3D operator-(const Vector3D& other) const;
    Vector3D operator*(double scalar) const;
    
private:
    double x_, y_, z_;
};

class Matrix3D {
public:
    Matrix3D();  // Identity
    Matrix3D(const std::array<double, 9>& values);  // Row-major
    
    Vector3D operator*(const Vector3D& vec) const;
    Matrix3D operator*(const Matrix3D& other) const;
    Matrix3D transpose() const;
    Matrix3D inverse() const;
    
    Vector3D column(size_t col) const;
    Vector3D row(size_t row) const;
    
    std::array<double, 9> as_array() const;  // Row-major
    
private:
    std::array<double, 9> data_;  // Row-major: [r11, r12, r13, r21, r22, r23, r31, r32, r33]
};

class LeastSquaresFitter {
public:
    struct FitResult {
        Matrix3D rotation;
        Vector3D translation;
        double rms;
    };
    
    FitResult fit(
        const std::vector<Vector3D>& points1,
        const std::vector<Vector3D>& points2
    );
};

} // namespace x3dna::geometry
```

---

## Applications Layer

### find_pair_app
```cpp
#include <x3dna/core/Structure.hpp>
#include <x3dna/io/PdbParser.hpp>
#include <x3dna/protocols/FindPairProtocol.hpp>
#include <x3dna/config/ConfigManager.hpp>

int main(int argc, char* argv[]) {
    // Parse command line
    // ...
    
    // Load configuration
    auto& config = ConfigManager::instance();
    config.load_from_file("config.json");
    
    // Parse PDB
    PdbParser parser;
    parser.set_include_hetatm(include_hetatm);
    Structure structure = parser.parse_file(pdb_path);
    
    // Execute protocol
    FindPairProtocol protocol;
    protocol.set_config_manager(config);
    protocol.execute(structure);
    
    // Output results
    // ...
    
    return 0;
}
```

### analyze_app
```cpp
#include <x3dna/core/Structure.hpp>
#include <x3dna/io/PdbParser.hpp>
#include <x3dna/protocols/AnalyzeProtocol.hpp>

int main(int argc, char* argv[]) {
    // Parse command line and input file
    // ...
    
    // Load structure and base pairs
    Structure structure = load_from_input_file(input_path);
    
    // Execute protocol
    AnalyzeProtocol protocol;
    protocol.execute(structure);
    
    // Output parameters
    // ...
    
    return 0;
}
```

---

## Testing Strategy

### Test Structure

```
tests/
├── unit/                       # Unit tests
│   ├── core/
│   │   ├── test_atom.cpp
│   │   ├── test_residue.cpp
│   │   ├── test_chain.cpp
│   │   ├── test_structure.cpp
│   │   ├── test_basepair.cpp
│   │   └── test_referenceframe.cpp
│   ├── geometry/
│   │   ├── test_vector3d.cpp
│   │   ├── test_matrix3d.cpp
│   │   └── test_leastsquares.cpp
│   ├── algorithms/
│   │   ├── test_baseframecalculator.cpp
│   │   ├── test_basepairfinder.cpp
│   │   ├── test_parametercalculator.cpp
│   │   └── test_hydrogenbondvalidator.cpp
│   └── io/
│       └── test_pdbparser.cpp
├── integration/                # Integration tests
│   ├── test_findpair_protocol.cpp
│   └── test_analyze_protocol.cpp
└── regression/                 # Regression tests
    ├── test_vs_json_legacy.cpp
    └── compare_json.py         # Python script to compare outputs
```

### Regression Testing Against JSON

```cpp
// tests/regression/test_vs_json_legacy.cpp

#include <gtest/gtest.h>
#include <x3dna/core/Structure.hpp>
#include <x3dna/io/PdbParser.hpp>
#include <x3dna/protocols/FindPairProtocol.hpp>
#include <x3dna/protocols/AnalyzeProtocol.hpp>
#include <nlohmann/json.hpp>

class JsonRegressionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Load JSON reference files
        std::ifstream json_file("data/json_legacy/1EHZ.json");
        json_file >> reference_json_;
    }
    
    nlohmann::json reference_json_;
};

TEST_F(JsonRegressionTest, ReferenceFramesMatch) {
    // Parse PDB
    PdbParser parser;
    Structure structure = parser.parse_file("data/pdb/1EHZ.pdb");
    
    // Execute find_pair
    FindPairProtocol find_pair;
    find_pair.execute(structure);
    
    // Compare reference frames
    auto ref_frames = reference_json_["ref_frames"];
    for (size_t i = 0; i < structure.base_pairs().size(); ++i) {
        const auto& pair = structure.base_pairs()[i];
        const auto& ref_frame = ref_frames[i];
        
        // Compare rotation matrices
        auto rotation = pair.frame1().rotation().as_array();
        auto ref_rotation = ref_frame["orien1"];
        for (size_t j = 0; j < 9; ++j) {
            EXPECT_NEAR(rotation[j], ref_rotation[j], 0.001);
        }
        
        // Compare origins
        auto origin = pair.frame1().origin().as_array();
        auto ref_origin = ref_frame["org1"];
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(origin[j], ref_origin[j], 0.001);
        }
    }
}

TEST_F(JsonRegressionTest, StepParametersMatch) {
    // Parse and analyze
    PdbParser parser;
    Structure structure = parser.parse_file("data/pdb/1EHZ.pdb");
    
    FindPairProtocol find_pair;
    find_pair.execute(structure);
    
    AnalyzeProtocol analyze;
    analyze.execute(structure);
    
    // Compare step parameters
    auto ref_params = reference_json_["step_parameters"];
    const auto& params = analyze.step_parameters();
    
    for (size_t i = 0; i < params.size(); ++i) {
        const auto& param = params[i];
        const auto& ref_param = ref_params[i];
        
        EXPECT_NEAR(param.shift, ref_param["shift"], 0.001);
        EXPECT_NEAR(param.slide, ref_param["slide"], 0.001);
        EXPECT_NEAR(param.rise, ref_param["rise"], 0.001);
        EXPECT_NEAR(param.tilt, ref_param["tilt"], 0.001);
        EXPECT_NEAR(param.roll, ref_param["roll"], 0.001);
        EXPECT_NEAR(param.twist, ref_param["twist"], 0.001);
    }
}
```

### Python Comparison Script

```python
# tests/regression/compare_json.py

import json
import sys
from pathlib import Path

def compare_json_files(legacy_path, new_path, tolerance=0.001):
    """Compare legacy JSON with new JSON output."""
    with open(legacy_path) as f:
        legacy = json.load(f)
    with open(new_path) as f:
        new = json.load(f)
    
    errors = []
    
    # Compare reference frames
    if "ref_frames" in legacy and "ref_frames" in new:
        for i, (leg_frame, new_frame) in enumerate(zip(legacy["ref_frames"], new["ref_frames"])):
            # Compare rotation matrices
            for j in range(9):
                diff = abs(leg_frame["orien1"][j] - new_frame["orien1"][j])
                if diff > tolerance:
                    errors.append(f"Frame {i}, orien1[{j}]: diff={diff}")
    
    # Compare step parameters
    if "step_parameters" in legacy and "step_parameters" in new:
        for i, (leg_param, new_param) in enumerate(zip(legacy["step_parameters"], new["step_parameters"])):
            for param_name in ["shift", "slide", "rise", "tilt", "roll", "twist"]:
                diff = abs(leg_param[param_name] - new_param[param_name])
                if diff > tolerance:
                    errors.append(f"Step {i}, {param_name}: diff={diff}")
    
    return errors

if __name__ == "__main__":
    legacy_path = sys.argv[1]
    new_path = sys.argv[2]
    
    errors = compare_json_files(legacy_path, new_path)
    
    if errors:
        print(f"Found {len(errors)} differences:")
        for error in errors[:10]:  # Print first 10
            print(f"  {error}")
        sys.exit(1)
    else:
        print("All values match within tolerance!")
        sys.exit(0)
```

---

## JSON Serialization & Legacy Format Compatibility

### Overview

All core structure classes support JSON serialization using the `nlohmann/json` library. The implementation provides two formats:

1. **Modern format**: Clean, structured JSON optimized for the new C++ API
2. **Legacy format**: Exact match with `data/json_legacy/*.json` files for regression testing

### JSON Library Dependency

```cpp
#include <nlohmann/json.hpp>
using json = nlohmann::json;
```

### Legacy JSON Structure Analysis

Based on analysis of `data/json_legacy/*.json` files, the legacy format includes:

#### Top-Level Structure
```json
{
  "pdb_file": "data/pdb/100D.pdb",
  "pdb_name": "100D",
  "calculations": [
    // Array of calculation records
  ]
}
```

#### Calculation Record Types

1. **base_frame_calc**: Base frame calculation metadata
   ```json
   {
     "type": "base_frame_calc",
     "residue_idx": 1,
     "base_type": "C",
     "standard_template": "/path/to/Atomic_C.pdb",
     "rms_fit": 0.012090,
     "num_matched_atoms": 6,
     "matched_atoms": [" N3 ", " C2 ", ...]
   }
   ```

2. **ls_fitting**: Least-squares fitting results
   ```json
   {
     "type": "ls_fitting",
     "residue_idx": 1,
     "num_points": 6,
     "rms_fit": 0.012090,
     "rotation_matrix": [[...], [...], [...]],
     "translation": [x, y, z]
   }
   ```

3. **frame_calc**: Complete frame calculation with matched coordinates
   ```json
   {
     "type": "frame_calc",
     "residue_idx": 1,
     "base_type": "C",
     "template_file": "/path/to/Atomic_C.pdb",
     "rms_fit": 0.012090,
     "num_matched_atoms": 6,
     "matched_coordinates": [
       {
         "atom_idx": 1,
         "std_xyz": [x, y, z],
         "exp_xyz": [x, y, z]
       },
       ...
     ]
   }
   ```

4. **base_pair**: Base pair information
   ```json
   {
     "type": "base_pair",
     "base_i": 1,
     "base_j": 72,
     "bp_type": "AT",
     "dir_xyz": [x, y, z],
     "orien_i": [[...], [...], [...]],
     "orien_j": [[...], [...], [...]],
     "org_i": [x, y, z],
     "org_j": [x, y, z]
   }
   ```

5. **pair_validation**: Base pair validation results
   ```json
   {
     "type": "pair_validation",
     "base_i": 1,
     "base_j": 72,
     "is_valid": 1,
     "bp_type_id": 2,
     "direction_vectors": {
       "dir_x": 0.798801,
       "dir_y": 0.487792,
       "dir_z": -0.979314
     },
     "calculated_values": {
       "dorg": 8.234567,
       "d_v": 0.123456,
       "plane_angle": 12.345678,
       "dNN": 10.234567,
       "quality_score": 8.456789
     },
     "validation_checks": {
       "distance_check": true,
       "d_v_check": true,
       "plane_angle_check": true,
       "dNN_check": true
     },
     "thresholds": {...}
   }
   ```

6. **hbond_list**: Hydrogen bond information
   ```json
   {
     "type": "hbond_list",
     "base_i": 1,
     "base_j": 72,
     "num_hbonds": 2,
     "hb_info_string": "[2] N1-N3 2.95 N6-O4 2.87",
     "hbonds": [
       {
         "hbond_idx": 1,
         "donor_atom": " N1 ",
         "acceptor_atom": " N3 ",
         "distance": 2.950000,
         "type": "-",
         "linkage_type": 18
       },
       ...
     ]
   }
   ```

7. **bpstep_params**: Base pair step parameters
   ```json
   {
     "type": "bpstep_params",
     "bp_idx1": 1,
     "bp_idx2": 2,
     "params": {
       "Shift": 0.303694,
       "Slide": -1.587609,
       "Rise": 3.311063,
       "Tilt": 2.300071,
       "Roll": 7.464288,
       "Twist": 33.206652
     },
     "mst_org": [x, y, z],
     "mst_orien": [[...], [...], [...]]
   }
   ```

8. **helical_params**: Helical parameters
   ```json
   {
     "type": "helical_params",
     "bp_idx1": 1,
     "bp_idx2": 2,
     "params": [x_displacement, y_displacement, rise, inclination, tip, twist],
     "mst_orgH": [x, y, z],
     "mst_orienH": [[...], [...], [...]]
   }
   ```

9. **pdb_atoms**: PDB atom data
   ```json
   {
     "type": "pdb_atoms",
     "num_atoms": 100,
     "atoms": [
       {
         "atom_idx": 1,
         "atom_name": " C1'",
         "residue_name": "  A",
         "chain_id": "A",
         "residue_seq": 1,
         "xyz": [x, y, z]
       },
       ...
     ]
   }
   ```

10. **ref_frame**: Reference frame for a residue
    ```json
    {
      "type": "ref_frame",
      "residue_idx": 1,
      "orien": [[...], [...], [...]],
      "org": [x, y, z]
    }
    ```

11. **base_pairs**: Base pair indices
    ```json
    {
      "type": "base_pairs",
      "ds": 1,
      "num_bp": 10,
      "pair_num": [[1, 72], [2, 71], ...]
    }
    ```

12. **bp_sequence**: Base pair sequence
    ```json
    {
      "type": "bp_sequence",
      "ds": 1,
      "num_bp": 10,
      "bp_seq": [["A", "T"], ["C", "G"], ...]
    }
    ```

13. **ry_classification**: Purine/Pyrimidine classification
    ```json
    {
      "type": "ry_classification",
      "num_residue": 20,
      "ry_values": [1, 0, 1, 0, ...]
    }
    ```

14. **ring_atoms**: Ring atom indices
    ```json
    {
      "type": "ring_atoms",
      "residue_idx": 1,
      "num_ring_atoms": 6,
      "ring_atom_indices": [14, 13, 11, 10, 17, 16]
    }
    ```

15. **distance_checks**: Distance/angle measurements
    ```json
    {
      "type": "distance_checks",
      "base_i": 1,
      "base_j": 72,
      "values": {
        "dorg": 8.234567,
        "dNN": 10.234567,
        "plane_angle": 12.345678,
        "d_v": 0.123456,
        "overlap_area": 0.456789
      }
    }
    ```

### Implementation Strategy

#### 1. Core Classes JSON Support

All core classes implement:
- `to_json()`: Modern format
- `to_json_legacy()`: Legacy format (matches `data/json_legacy/*.json`)
- `from_json()`: Parse modern format
- `from_json_legacy()`: Parse legacy format

#### 2. Example Implementation

```cpp
// Atom JSON serialization
nlohmann::json Atom::to_json() const {
    return {
        {"name", name_},
        {"position", position_.to_json()}
    };
}

nlohmann::json Atom::to_json_legacy() const {
    return {
        {"atom_name", name_},
        {"xyz", {position_.x(), position_.y(), position_.z()}}
    };
}

// ReferenceFrame JSON serialization
nlohmann::json ReferenceFrame::to_json_legacy() const {
    auto rot = rotation_.as_array();
    auto org_arr = origin_.as_array();
    
    // Convert to 3x3 nested array format
    std::vector<std::vector<double>> orien_matrix = {
        {rot[0], rot[1], rot[2]},
        {rot[3], rot[4], rot[5]},
        {rot[6], rot[7], rot[8]}
    };
    
    return {
        {"orien", orien_matrix},
        {"org", {org_arr[0], org_arr[1], org_arr[2]}}
    };
}

// BasePairStepParameters JSON serialization
nlohmann::json BasePairStepParameters::to_json_legacy() const {
    return {
        {"params", {
            {"Shift", shift},
            {"Slide", slide},
            {"Rise", rise},
            {"Tilt", tilt},
            {"Roll", roll},
            {"Twist", twist}
        }}
    };
}
```

#### 3. Structure-Level JSON Export

```cpp
nlohmann::json Structure::to_json_legacy() const {
    nlohmann::json result = {
        {"pdb_file", pdb_file_path_},
        {"pdb_name", pdb_id_},
        {"calculations", nlohmann::json::array()}
    };
    
    auto& calculations = result["calculations"];
    
    // Add all calculation records in order
    for (const auto& residue : all_residues()) {
        if (residue->reference_frame()) {
            calculations.push_back({
                {"type", "ref_frame"},
                {"residue_idx", residue->index()},
                {"orien", residue->reference_frame()->to_json_legacy()["orien"]},
                {"org", residue->reference_frame()->to_json_legacy()["org"]}
            });
        }
    }
    
    for (const auto& pair : base_pairs_) {
        calculations.push_back(pair.to_json_legacy());
    }
    
    // Add step parameters, etc.
    
    return result;
}
```

#### 4. Regression Testing Support

The JSON serialization enables:

1. **Reading legacy JSON files**:
   ```cpp
   Structure structure = Structure::from_json_legacy(legacy_json);
   ```

2. **Writing in legacy format**:
   ```cpp
   nlohmann::json output = structure.to_json_legacy();
   std::ofstream file("output.json");
   file << output.dump(2);  // Pretty print with 2-space indent
   ```

3. **Comparing outputs**:
   ```cpp
   auto our_json = structure.to_json_legacy();
   auto legacy_json = load_json_file("data/json_legacy/100D.json");
   
   // Compare specific fields
   assert(our_json["calculations"][0] == legacy_json["calculations"][0]);
   ```

### JSON Writer Utility Class

For convenience, provide a `JsonWriter` class that matches the original `json_writer.c` functionality:

```cpp
namespace x3dna::io {

class JsonWriter {
public:
    void init(const std::filesystem::path& pdb_file);
    void close();
    
    // Record functions matching original API
    void record_base_frame_calc(size_t residue_idx, char base_type, 
                                 const std::string& template_file,
                                 double rms_fit, 
                                 const std::vector<std::string>& matched_atoms);
    
    void record_ls_fitting(size_t residue_idx, size_t num_points, double rms_fit,
                          const Matrix3D& rotation, const Vector3D& translation);
    
    void record_base_pair(const BasePair& pair);
    void record_bpstep_params(size_t bp_idx1, size_t bp_idx2,
                             const BasePairStepParameters& params,
                             const ReferenceFrame& midstep_frame);
    
    void record_helical_params(size_t bp_idx1, size_t bp_idx2,
                               const HelicalParameters& params,
                               const ReferenceFrame& midstep_frame);
    
    void record_pair_validation(const BasePairFinder::ValidationResult& result);
    void record_hbond_list(const BasePair& pair);
    void record_pdb_atoms(const Structure& structure);
    void record_all_ref_frames(const Structure& structure);
    // ... etc
    
    nlohmann::json get_json() const { return json_data_; }
    void write_to_file(const std::filesystem::path& path) const;
    
private:
    nlohmann::json json_data_;
    std::filesystem::path pdb_file_;
    std::string pdb_name_;
};

} // namespace x3dna::io
```

### Benefits

1. **Regression Testing**: Can read/write legacy JSON format for comparison
2. **Data Interchange**: Easy serialization for storage/transmission
3. **Debugging**: Export structures at any point in processing
4. **Compatibility**: Maintain compatibility with existing tools that read legacy JSON
5. **Testing**: Generate test fixtures from real data

---

## Migration Path

### Phase 1: Foundation (Weeks 1-2)
1. Set up CMake build system
2. Implement core geometry classes (Vector3D, Matrix3D)
3. Implement core domain objects (Atom, Residue, Chain, Structure)
4. Basic unit tests for core classes

### Phase 2: I/O Layer (Weeks 3-4)
1. Implement PdbParser
2. Implement JSON reader/writer
3. Test with sample PDB files
4. Compare parsed structures with original

### Phase 3: Algorithms - Part 1 (Weeks 5-7)
1. Implement LeastSquaresFitter
2. Implement BaseFrameCalculator
3. Test frame calculation against JSON reference
4. Implement HydrogenBondValidator

### Phase 4: Algorithms - Part 2 (Weeks 8-10)
1. Implement BasePairFinder
2. Test pair finding against JSON reference
3. Implement HelixDetector
4. Test helix detection

### Phase 5: Parameter Calculation (Weeks 11-12)
1. Implement ParameterCalculator
2. Test against JSON reference (critical!)
3. Verify all 6 parameters match within tolerance

### Phase 6: Protocols (Weeks 13-14)
1. Implement ProtocolBase
2. Implement FindPairProtocol
3. Implement AnalyzeProtocol
4. Integration tests

### Phase 7: Applications (Week 15)
1. Implement find_pair_app
2. Implement analyze_app
3. Command-line argument parsing
4. Output formatting

### Phase 8: Testing & Validation (Weeks 16-17)
1. Comprehensive regression testing
2. Compare all outputs with JSON legacy files
3. Fix any discrepancies
4. Performance testing

### Phase 9: Documentation & Polish (Week 18)
1. API documentation (Doxygen)
2. User guide
3. Example code
4. Final review

---

## Key Design Decisions

### 1. Memory Management
- Use smart pointers (`std::unique_ptr`, `std::shared_ptr`) instead of raw pointers
- RAII for all resources
- STL containers instead of manual arrays

### 2. Indexing
- Use 0-based indexing throughout (modern C++ standard)
- Provide conversion utilities if needed for legacy compatibility

### 3. Error Handling
- Use exceptions for error conditions
- Provide detailed error messages with context
- No global error state

### 4. Configuration
- Singleton ConfigManager for global settings
- Immutable ParameterThresholds structure
- JSON-based configuration files

### 5. Testing
- Comprehensive unit tests for all classes
- Integration tests for protocols
- Regression tests comparing with JSON legacy files
- Use Google Test framework

### 6. API Design
- Const-correctness throughout
- Clear separation between mutable and immutable operations
- Factory functions for complex object creation
- Builder pattern for complex configurations

---

## Success Criteria

1. ✅ All core classes implemented with proper OOP design
2. ✅ PDB parser generates Structure objects correctly
3. ✅ Base pair detection matches original algorithm
4. ✅ Parameter calculation matches JSON reference within 0.001 tolerance
5. ✅ Comprehensive test coverage (>80%)
6. ✅ All regression tests pass
7. ✅ Clean, modern C++ code following best practices
8. ✅ Full API documentation
9. ✅ Performance comparable to or better than original

---

## Next Steps

1. Review and approve this design plan
2. Set up project structure and build system
3. Begin Phase 1 implementation
4. Establish continuous integration for testing
5. Create development branch and start implementation

---

## Appendix: Class Dependency Graph

```
ConfigManager
    ↑
    ├── PdbParser
    │       ↓
    │   Structure
    │       ├── Chain
    │       │   └── Residue
    │       │       └── Atom
    │       └── BasePair
    │           └── ReferenceFrame
    │
    ├── BaseFrameCalculator
    │       ├── LeastSquaresFitter
    │       └── GeometryUtils
    │
    ├── BasePairFinder
    │       └── HydrogenBondValidator
    │
    ├── ParameterCalculator
    │       └── GeometryUtils
    │
    ├── HelixDetector
    │
    └── Protocols
            ├── FindPairProtocol
            └── AnalyzeProtocol
```

---

*This plan provides a comprehensive roadmap for modernizing the X3DNA codebase into a maintainable, testable, and extensible C++ library.*

