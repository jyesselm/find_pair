# Stage 12: Advanced Features & Extensions

## Objectives

Add advanced analysis features, complete remaining functionality, implement performance optimizations, and create extension points for future development.

## Duration

**2-3 weeks**

## Dependencies

- ✅ Stage 0-11: Core functionality complete and refactored
- ✅ Clean OOP architecture from Stage 11

---

## Tasks

### Task 12.1: Complete Helix Detection System

**Current State**: `HelixDetector` exists but is minimal. Full helix analysis is needed.

**Solution**: Comprehensive helix detection with context analysis.

```cpp
/**
 * @class HelixDetector
 * @brief Detects helical regions from base pairs
 * 
 * Implements legacy helix detection:
 * - Consecutive base pair grouping
 * - 5'→3' strand orientation
 * - Helix numbering
 * - Circular structure handling
 */
class HelixDetector {
public:
    /**
     * @brief Detect helices from found base pairs
     * @param pairs Vector of base pairs
     * @param structure Structure for context
     * @return Detected helices
     */
    std::vector<Helix> detect(
        const std::vector<core::BasePair>& pairs,
        const core::Structure& structure
    );

    /**
     * @brief Reorder pairs within helices for 5'→3' consistency
     */
    void reorder_for_five_to_three(
        std::vector<core::BasePair>& pairs,
        const core::Structure& structure
    );

    /**
     * @brief Analyze base pair context (stacking, neighboring)
     */
    BasePairContext analyze_context(
        const core::BasePair& pair,
        const std::vector<core::BasePair>& all_pairs,
        const core::Structure& structure
    ) const;

    /**
     * @brief Group consecutive pairs into stems
     */
    std::vector<Stem> find_stems(
        const std::vector<core::BasePair>& pairs,
        const core::Structure& structure
    ) const;

    // Options
    void set_min_helix_length(size_t min_length) { min_helix_length_ = min_length; }
    void set_allow_circular(bool allow) { allow_circular_ = allow; }
    void set_max_bulge_size(size_t max_size) { max_bulge_size_ = max_size; }

private:
    size_t min_helix_length_ = 2;
    bool allow_circular_ = true;
    size_t max_bulge_size_ = 1;

    bool are_consecutive(const core::BasePair& p1, const core::BasePair& p2,
                         const core::Structure& structure) const;
    bool is_stacked(const core::BasePair& p1, const core::BasePair& p2) const;
};

/**
 * @struct Helix
 * @brief Represents a helical region
 */
struct Helix {
    size_t helix_id;
    std::vector<size_t> pair_indices;
    size_t strand1_start_residue;
    size_t strand1_end_residue;
    size_t strand2_start_residue;
    size_t strand2_end_residue;
    bool is_circular;
    double average_twist;
    double average_rise;
};

/**
 * @struct Stem
 * @brief Represents a double-stranded stem region
 */
struct Stem {
    std::vector<size_t> pair_indices;
    bool contains_bulge;
    std::vector<size_t> bulge_residues;
};

/**
 * @struct BasePairContext
 * @brief Context information for a base pair
 */
struct BasePairContext {
    std::optional<size_t> stacked_above;  // Index of pair stacked above
    std::optional<size_t> stacked_below;  // Index of pair stacked below
    bool is_terminal;
    bool is_in_helix;
    size_t helix_position;  // Position within helix (0 = first)
};
```

**Files**:
- `include/x3dna/algorithms/helix_detector.hpp` (update)
- `src/x3dna/algorithms/helix_detector.cpp` (complete implementation)
- `tests/unit/algorithms/test_helix_detector.cpp`

**Deliverable**: Complete helix detection system

---

### Task 12.2: Base Pair Step Analysis Pipeline

**Current State**: `ParameterCalculator` exists but not integrated into full pipeline.

**Solution**: Complete step parameter analysis with helical parameters.

```cpp
/**
 * @class StepAnalyzer
 * @brief Analyzes consecutive base pair steps
 * 
 * Calculates for each step:
 * - 6 step parameters (Shift, Slide, Rise, Tilt, Roll, Twist)
 * - 6 helical parameters (X-disp, Y-disp, h-Rise, Incl, Tip, h-Twist)
 * - Midstep reference frames
 */
class StepAnalyzer {
public:
    /**
     * @brief Analyze all steps in ordered pairs
     */
    std::vector<StepAnalysisResult> analyze_all_steps(
        const std::vector<core::BasePair>& ordered_pairs
    );

    /**
     * @brief Analyze single step between two pairs
     */
    StepAnalysisResult analyze_step(
        const core::BasePair& pair1,
        const core::BasePair& pair2
    );

    /**
     * @brief Calculate groove widths
     */
    GrooveWidths calculate_groove_widths(
        const std::vector<core::BasePair>& pairs,
        const core::Structure& structure
    ) const;

private:
    ParameterCalculator param_calc_;
};

/**
 * @struct StepAnalysisResult
 */
struct StepAnalysisResult {
    // Step parameters
    double shift, slide, rise;
    double tilt, roll, twist;
    
    // Helical parameters
    double x_disp, y_disp, h_rise;
    double inclination, tip, h_twist;
    
    // Reference frames
    core::ReferenceFrame midstep_frame;
    core::ReferenceFrame helical_midstep_frame;
    
    // Pair indices
    size_t pair1_idx;
    size_t pair2_idx;
};

/**
 * @struct GrooveWidths
 */
struct GrooveWidths {
    std::vector<double> minor_groove;
    std::vector<double> major_groove;
};
```

**Files**:
- `include/x3dna/algorithms/step_analyzer.hpp`
- `src/x3dna/algorithms/step_analyzer.cpp`
- `tests/unit/algorithms/test_step_analyzer.cpp`

**Deliverable**: Complete step analysis pipeline

---

### Task 12.3: Parallel Processing Support

**Current State**: All processing is single-threaded.

**Solution**: Add parallel processing for independent operations.

```cpp
/**
 * @class ParallelProcessor
 * @brief Parallel execution for independent operations
 * 
 * Uses std::async or thread pool for:
 * - Frame calculation (per residue)
 * - Pair validation (Phase 1)
 * - Parameter calculation (per step)
 */
class ParallelProcessor {
public:
    /**
     * @brief Configure thread count (0 = hardware concurrency)
     */
    void set_thread_count(size_t count);
    size_t thread_count() const;

    /**
     * @brief Calculate frames in parallel
     */
    void calculate_frames_parallel(
        core::Structure& structure,
        BaseFrameCalculator& calculator
    );

    /**
     * @brief Build pair candidate cache in parallel
     */
    void build_cache_parallel(
        PairCandidateCache& cache,
        const core::Structure& structure,
        const BasePairValidator& validator
    );

    /**
     * @brief Process multiple PDB files in parallel
     */
    template<typename ResultType>
    std::vector<ResultType> process_files_parallel(
        const std::vector<std::filesystem::path>& files,
        std::function<ResultType(const std::filesystem::path&)> processor
    );

private:
    size_t thread_count_ = 0;  // 0 = auto-detect
    
    // Thread-safe accumulator for results
    template<typename T>
    class ThreadSafeAccumulator { /* ... */ };
};

/**
 * @class BatchProcessor
 * @brief Process multiple structures efficiently
 */
class BatchProcessor {
public:
    struct BatchResult {
        std::string pdb_id;
        bool success;
        std::string error_message;
        std::vector<core::BasePair> pairs;
        std::chrono::milliseconds duration;
    };

    /**
     * @brief Process directory of PDB files
     */
    std::vector<BatchResult> process_directory(
        const std::filesystem::path& pdb_dir,
        const std::filesystem::path& output_dir
    );

    void set_parallel(bool enable) { parallel_ = enable; }
    void set_progress_callback(std::function<void(size_t, size_t)> cb);

private:
    bool parallel_ = true;
    std::function<void(size_t, size_t)> progress_callback_;
};
```

**Files**:
- `include/x3dna/algorithms/parallel_processor.hpp`
- `include/x3dna/algorithms/batch_processor.hpp`
- `src/x3dna/algorithms/parallel_processor.cpp`
- `src/x3dna/algorithms/batch_processor.cpp`

**Deliverable**: Parallel processing with 4-8x speedup on multi-core

---

### Task 12.4: Plugin Architecture for Custom Strategies

**Current State**: Strategies must be compiled into library.

**Solution**: Plugin system for runtime strategy loading.

```cpp
/**
 * @class StrategyRegistry
 * @brief Registry for pair selection strategies
 * 
 * Allows registration of custom strategies at runtime.
 */
class StrategyRegistry {
public:
    static StrategyRegistry& instance();

    /**
     * @brief Register a strategy factory
     */
    void register_strategy(
        const std::string& name,
        std::function<std::unique_ptr<IPairSelectionStrategy>()> factory
    );

    /**
     * @brief Create strategy by name
     */
    std::unique_ptr<IPairSelectionStrategy> create(
        const std::string& name
    ) const;

    /**
     * @brief List available strategies
     */
    std::vector<std::string> available_strategies() const;

private:
    std::map<std::string, 
             std::function<std::unique_ptr<IPairSelectionStrategy>()>> factories_;
};

/**
 * @brief Macro for easy strategy registration
 */
#define REGISTER_STRATEGY(StrategyClass, name) \
    static bool _##StrategyClass##_registered = []() { \
        StrategyRegistry::instance().register_strategy( \
            name, []() { return std::make_unique<StrategyClass>(); } \
        ); \
        return true; \
    }();

// Usage in .cpp file:
REGISTER_STRATEGY(MutualBestPairStrategy, "mutual_best")
REGISTER_STRATEGY(AllPairsStrategy, "all_pairs")
```

**Files**:
- `include/x3dna/algorithms/strategy_registry.hpp`
- `src/x3dna/algorithms/strategy_registry.cpp`

**Deliverable**: Plugin architecture for custom strategies

---

### Task 12.5: Enhanced Validation Options

**Current State**: Fixed validation thresholds.

**Solution**: Configurable validation profiles.

```cpp
/**
 * @class ValidationProfile
 * @brief Named validation parameter sets
 */
class ValidationProfile {
public:
    static ValidationParameters strict();      // Tight thresholds
    static ValidationParameters relaxed();     // Looser thresholds
    static ValidationParameters legacy();      // Exact legacy values
    static ValidationParameters rna_only();    // RNA-optimized
    static ValidationParameters dna_only();    // DNA-optimized

    /**
     * @brief Load profile from JSON file
     */
    static ValidationParameters from_file(
        const std::filesystem::path& path
    );

    /**
     * @brief Save profile to JSON file
     */
    static void to_file(
        const ValidationParameters& params,
        const std::filesystem::path& path
    );
};

/**
 * @class ValidationBuilder
 * @brief Fluent builder for validation parameters
 */
class ValidationBuilder {
public:
    ValidationBuilder& max_origin_distance(double d);
    ValidationBuilder& max_vertical_distance(double d);
    ValidationBuilder& max_plane_angle(double a);
    ValidationBuilder& min_hbonds(int n);
    ValidationBuilder& hbond_distance_range(double min, double max);
    
    ValidationParameters build() const;

private:
    ValidationParameters params_;
};

// Usage:
auto params = ValidationBuilder()
    .max_origin_distance(12.0)
    .max_plane_angle(45.0)
    .min_hbonds(2)
    .build();
```

**Files**:
- `include/x3dna/algorithms/validation_profile.hpp`
- `src/x3dna/algorithms/validation_profile.cpp`

**Deliverable**: Configurable validation with preset profiles

---

### Task 12.6: Statistics and Analytics

**Current State**: No statistical analysis built in.

**Solution**: Analytics module for structure analysis.

```cpp
/**
 * @class StructureAnalytics
 * @brief Statistical analysis of structure properties
 */
class StructureAnalytics {
public:
    /**
     * @brief Analyze base pair statistics
     */
    struct BasePairStats {
        size_t total_pairs;
        std::map<std::string, size_t> pair_type_counts;  // AT, GC, etc.
        double mean_quality_score;
        double std_quality_score;
        size_t watson_crick_count;
        size_t wobble_count;
        size_t non_canonical_count;
    };
    BasePairStats analyze_pairs(const std::vector<core::BasePair>& pairs) const;

    /**
     * @brief Analyze helix statistics
     */
    struct HelixStats {
        size_t total_helices;
        double mean_helix_length;
        double max_helix_length;
        double mean_twist;
        double mean_rise;
    };
    HelixStats analyze_helices(const std::vector<Helix>& helices) const;

    /**
     * @brief Analyze step parameter distributions
     */
    struct StepStats {
        struct ParamStats {
            double mean, stddev, min, max;
        };
        ParamStats shift, slide, rise;
        ParamStats tilt, roll, twist;
    };
    StepStats analyze_steps(const std::vector<StepAnalysisResult>& steps) const;

    /**
     * @brief Generate summary report
     */
    std::string generate_report(
        const core::Structure& structure,
        const std::vector<core::BasePair>& pairs,
        const std::vector<Helix>& helices
    ) const;
};

/**
 * @class ComparisonAnalytics
 * @brief Compare two structures or analysis results
 */
class ComparisonAnalytics {
public:
    struct ComparisonResult {
        size_t pairs_in_common;
        size_t pairs_only_in_first;
        size_t pairs_only_in_second;
        double jaccard_index;
        std::vector<std::pair<size_t, size_t>> pair_differences;
    };

    ComparisonResult compare_pairs(
        const std::vector<core::BasePair>& pairs1,
        const std::vector<core::BasePair>& pairs2
    ) const;
};
```

**Files**:
- `include/x3dna/analytics/structure_analytics.hpp`
- `include/x3dna/analytics/comparison_analytics.hpp`
- `src/x3dna/analytics/structure_analytics.cpp`
- `src/x3dna/analytics/comparison_analytics.cpp`

**Deliverable**: Analytics module for structure statistics

---

### Task 12.7: Output Format Extensions

**Current State**: JSON output only.

**Solution**: Multiple output format support.

```cpp
/**
 * @class IOutputFormatter
 * @brief Interface for output formatting
 */
class IOutputFormatter {
public:
    virtual ~IOutputFormatter() = default;

    virtual void write_pairs(
        const std::vector<core::BasePair>& pairs,
        std::ostream& out
    ) const = 0;

    virtual void write_helices(
        const std::vector<Helix>& helices,
        std::ostream& out
    ) const = 0;

    virtual void write_parameters(
        const std::vector<StepAnalysisResult>& steps,
        std::ostream& out
    ) const = 0;

    virtual std::string extension() const = 0;
};

/**
 * @brief JSON output (default)
 */
class JsonFormatter : public IOutputFormatter {
    // Current JSON format
};

/**
 * @brief CSV output for spreadsheet analysis
 */
class CsvFormatter : public IOutputFormatter {
    // CSV with headers
};

/**
 * @brief Legacy X3DNA .out format
 */
class LegacyFormatter : public IOutputFormatter {
    // Original X3DNA output format
};

/**
 * @brief XML format
 */
class XmlFormatter : public IOutputFormatter {
    // Structured XML output
};

/**
 * @class OutputManager
 * @brief Manages output to multiple formats
 */
class OutputManager {
public:
    void add_formatter(std::unique_ptr<IOutputFormatter> formatter);
    
    void write_all(
        const std::filesystem::path& base_path,
        const std::vector<core::BasePair>& pairs,
        const std::vector<Helix>& helices,
        const std::vector<StepAnalysisResult>& steps
    ) const;

private:
    std::vector<std::unique_ptr<IOutputFormatter>> formatters_;
};
```

**Files**:
- `include/x3dna/io/output_formatter.hpp`
- `include/x3dna/io/json_formatter.hpp`
- `include/x3dna/io/csv_formatter.hpp`
- `src/x3dna/io/csv_formatter.cpp`

**Deliverable**: Multi-format output support

---

### Task 12.8: Performance Benchmarking Suite

**Current State**: No systematic performance testing.

**Solution**: Benchmarking framework.

```cpp
/**
 * @class Benchmark
 * @brief Performance benchmarking utility
 */
class Benchmark {
public:
    struct Result {
        std::string name;
        std::chrono::nanoseconds mean_time;
        std::chrono::nanoseconds std_dev;
        size_t iterations;
        double ops_per_second;
    };

    /**
     * @brief Run benchmark
     */
    template<typename Func>
    Result run(const std::string& name, Func&& func, size_t iterations = 100);

    /**
     * @brief Run benchmark suite
     */
    std::vector<Result> run_suite(
        const std::vector<std::pair<std::string, std::function<void()>>>& benchmarks
    );

    /**
     * @brief Generate report
     */
    std::string generate_report(const std::vector<Result>& results) const;
};

// Predefined benchmarks
namespace benchmarks {
    void benchmark_frame_calculation(const std::filesystem::path& pdb_file);
    void benchmark_pair_finding(const std::filesystem::path& pdb_file);
    void benchmark_validation(const std::filesystem::path& pdb_file);
    void benchmark_json_writing(const std::filesystem::path& json_file);
}
```

**Files**:
- `include/x3dna/benchmark/benchmark.hpp`
- `src/x3dna/benchmark/benchmark.cpp`
- `tools/run_benchmarks.cpp`

**Deliverable**: Benchmark suite with baseline measurements

---

## Testing Plan

### Unit Tests

| Component | Test File | Key Tests |
|-----------|-----------|-----------|
| HelixDetector | `test_helix_detector.cpp` | Consecutive detection, reordering |
| StepAnalyzer | `test_step_analyzer.cpp` | Parameter calculation |
| ParallelProcessor | `test_parallel_processor.cpp` | Thread safety, correctness |
| StrategyRegistry | `test_strategy_registry.cpp` | Registration, creation |
| ValidationProfile | `test_validation_profile.cpp` | Profiles, serialization |
| StructureAnalytics | `test_analytics.cpp` | Statistics accuracy |
| CsvFormatter | `test_csv_formatter.cpp` | Format correctness |

### Integration Tests

- [ ] Parallel processing produces same results as sequential
- [ ] All output formats contain equivalent information
- [ ] Custom strategies can be registered and used
- [ ] Batch processing handles errors gracefully

### Performance Tests

- [ ] Parallel frame calculation 4x faster on 4 cores
- [ ] Batch processing 8x faster with parallel on 8 cores
- [ ] Memory usage scales linearly with structure size

---

## Success Criteria

- [ ] HelixDetector produces correct helix assignments
- [ ] Step analysis matches legacy parameter values
- [ ] Parallel processing gives 4x+ speedup on 4 cores
- [ ] All output formats pass validation
- [ ] Analytics produce correct statistics
- [ ] Benchmark suite runs without errors
- [ ] All tests pass
- [ ] Documentation complete

---

## Files Created

```
include/x3dna/algorithms/
├── helix_detector.hpp (complete)
├── step_analyzer.hpp
├── parallel_processor.hpp
├── batch_processor.hpp
├── strategy_registry.hpp
├── validation_profile.hpp

include/x3dna/analytics/
├── structure_analytics.hpp
├── comparison_analytics.hpp

include/x3dna/io/
├── output_formatter.hpp
├── json_formatter.hpp
├── csv_formatter.hpp
├── xml_formatter.hpp

include/x3dna/benchmark/
├── benchmark.hpp

src/x3dna/algorithms/
├── helix_detector.cpp (complete)
├── step_analyzer.cpp
├── parallel_processor.cpp
├── batch_processor.cpp
├── strategy_registry.cpp
├── validation_profile.cpp

src/x3dna/analytics/
├── structure_analytics.cpp
├── comparison_analytics.cpp

src/x3dna/io/
├── csv_formatter.cpp
├── xml_formatter.cpp

src/x3dna/benchmark/
├── benchmark.cpp

tools/
├── run_benchmarks.cpp
├── batch_process.cpp
```

---

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Thread safety bugs | High | Extensive testing, thread sanitizers |
| Performance regression | Medium | Benchmark suite, profiling |
| Plugin stability | Medium | Clear API contracts, versioning |
| Output format incompatibility | Low | Format validation tests |

---

## Future Extensions (Post Stage 12)

These features are out of scope but would build on Stage 12:

1. **RNA Tertiary Structure**: Base triples, pseudoknots, A-minor motifs
2. **Protein-DNA/RNA Contacts**: Interface analysis
3. **Molecular Dynamics Analysis**: Trajectory processing
4. **Web Service API**: REST API for remote analysis
5. **GUI Application**: Visual structure explorer
6. **Machine Learning Integration**: ML-based pair prediction

---

*Estimated Completion: Week 22-24*

