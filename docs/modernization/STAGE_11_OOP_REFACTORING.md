# Stage 11: OOP Refactoring & Clean Architecture

## Objectives

Refactor the existing implementation to follow strong OOP principles, improve testability, and create an extensible architecture that separates concerns properly.

## Duration

**2-3 weeks**

## Dependencies

- ✅ Stage 0-10: Core functionality complete
- ✅ All algorithms working and validated against legacy

## Motivation

The current `BasePairFinder` is 1000+ lines and violates several OOP principles:
1. **Single Responsibility Violation**: Handles validation, caching, selection, recording, and quality scoring all in one class
2. **Tight Coupling**: `JsonWriter` passed through multiple layers; `ParameterCalculator` embedded
3. **Hidden State**: Phase 1 validation results stored in maps passed as parameters
4. **Legacy Mixing**: Modern and legacy indexing logic intertwined throughout

---

## Tasks

### Task 11.1: Extract Pair Candidate Cache (Single Responsibility)

**Problem**: Phase 1 validation results and bp_type_ids are computed in `find_best_pairs()` and stored in maps passed around as function parameters.

**Solution**: Create `PairCandidateCache` class to manage validated candidates.

```cpp
/**
 * @class PairCandidateCache
 * @brief Caches validation results for all residue pairs
 * 
 * Pre-computes and caches validation for all candidate pairs,
 * ensuring consistency between validation and selection phases.
 */
class PairCandidateCache {
public:
    struct CandidateInfo {
        ValidationResult validation;
        int bp_type_id;
        double adjusted_quality_score;
        bool is_valid() const { return validation.is_valid; }
    };

    /**
     * @brief Build cache for all valid pairs in structure
     * @param structure Structure with frames calculated
     * @param validator Validator to use for pair checking
     */
    void build(const core::Structure& structure, 
               const BasePairValidator& validator);

    /**
     * @brief Get cached result for a pair (order-independent)
     * @return CandidateInfo or nullopt if pair wasn't cached
     */
    std::optional<CandidateInfo> get(int legacy_idx1, int legacy_idx2) const;

    /**
     * @brief Check if pair exists in cache
     */
    bool contains(int legacy_idx1, int legacy_idx2) const;

    /**
     * @brief Get all valid candidates for a residue
     */
    std::vector<std::pair<int, CandidateInfo>> get_candidates_for(int legacy_idx) const;

    /**
     * @brief Iterator over all valid pairs
     */
    auto valid_pairs() const;

private:
    // Normalize pair key (smaller index first)
    static std::pair<int, int> normalize(int i, int j);
    
    std::map<std::pair<int, int>, CandidateInfo> cache_;
    std::map<int, std::vector<int>> candidates_by_residue_;
};
```

**Files**:
- `include/x3dna/algorithms/pair_candidate_cache.hpp`
- `src/x3dna/algorithms/pair_candidate_cache.cpp`

**Deliverable**: PairCandidateCache class with unit tests

---

### Task 11.2: Implement Strategy Pattern for Pair Finding

**Problem**: `find_best_pairs()` and `find_all_pairs()` have duplicated logic; adding new strategies requires modifying `BasePairFinder`.

**Solution**: Strategy pattern with `IPairSelectionStrategy` interface.

```cpp
/**
 * @class IPairSelectionStrategy
 * @brief Interface for pair selection strategies
 */
class IPairSelectionStrategy {
public:
    virtual ~IPairSelectionStrategy() = default;

    /**
     * @brief Select pairs from validated candidates
     * @param cache Validated pair candidates
     * @param residue_lookup Maps legacy indices to residues
     * @return Selected pairs in order
     */
    virtual std::vector<SelectedPair> select(
        const PairCandidateCache& cache,
        const ResidueIndexMap& residue_lookup
    ) const = 0;

    /**
     * @brief Get strategy name for debugging/logging
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Output from selection process
 */
struct SelectedPair {
    int legacy_idx1;
    int legacy_idx2;
    const PairCandidateCache::CandidateInfo* info;
};

/**
 * @class MutualBestPairStrategy
 * @brief Greedy mutual best match (legacy find_bestpair)
 */
class MutualBestPairStrategy : public IPairSelectionStrategy {
public:
    std::vector<SelectedPair> select(
        const PairCandidateCache& cache,
        const ResidueIndexMap& residue_lookup
    ) const override;

    std::string name() const override { return "mutual_best_pair"; }

private:
    std::optional<int> find_best_partner(
        int legacy_idx,
        const PairCandidateCache& cache,
        const std::unordered_set<int>& matched
    ) const;
};

/**
 * @class AllPairsStrategy
 * @brief Return all valid pairs (exhaustive)
 */
class AllPairsStrategy : public IPairSelectionStrategy {
public:
    std::vector<SelectedPair> select(
        const PairCandidateCache& cache,
        const ResidueIndexMap& residue_lookup
    ) const override;

    std::string name() const override { return "all_pairs"; }
};

/**
 * @class DistanceBasedStrategy
 * @brief Simple distance-based selection (placeholder for future)
 */
class DistanceBasedStrategy : public IPairSelectionStrategy {
    // ...
};
```

**Files**:
- `include/x3dna/algorithms/strategies/pair_selection_strategy.hpp`
- `include/x3dna/algorithms/strategies/mutual_best_pair_strategy.hpp`
- `include/x3dna/algorithms/strategies/all_pairs_strategy.hpp`
- `src/x3dna/algorithms/strategies/mutual_best_pair_strategy.cpp`
- `src/x3dna/algorithms/strategies/all_pairs_strategy.cpp`

**Deliverable**: Strategy interface + 2-3 concrete strategies with tests

---

### Task 11.3: Implement Observer Pattern for Recording

**Problem**: `JsonWriter*` passed through 5+ method calls; recording logic scattered throughout finder.

**Solution**: Observer pattern with `IPairFindingObserver` interface.

```cpp
/**
 * @class IPairFindingObserver
 * @brief Observer interface for pair finding events
 * 
 * Allows decoupled recording of validation and selection events.
 */
class IPairFindingObserver {
public:
    virtual ~IPairFindingObserver() = default;

    // Validation phase events
    virtual void on_pair_validated(
        int idx1, int idx2,
        const ValidationResult& result,
        int bp_type_id
    ) = 0;

    virtual void on_hbonds_found(
        int idx1, int idx2,
        const std::vector<core::hydrogen_bond>& hbonds
    ) = 0;

    // Selection phase events
    virtual void on_partner_search(
        int idx,
        const std::vector<std::tuple<int, bool, double, int>>& candidates,
        int best_partner,
        double best_score
    ) = 0;

    virtual void on_mutual_best_check(
        int idx1, int idx2,
        int best_for_1, int best_for_2,
        bool is_mutual, bool was_selected
    ) = 0;

    virtual void on_pair_selected(const core::BasePair& pair) = 0;

    virtual void on_iteration_complete(
        int iteration,
        size_t matched_count,
        const std::vector<std::pair<int, int>>& new_pairs
    ) = 0;
};

/**
 * @class JsonRecordingObserver
 * @brief Observer that records events to JsonWriter
 */
class JsonRecordingObserver : public IPairFindingObserver {
public:
    explicit JsonRecordingObserver(io::JsonWriter& writer);
    
    // Implement all observer methods...
private:
    io::JsonWriter& writer_;
};

/**
 * @class NullObserver
 * @brief No-op observer for when recording isn't needed
 */
class NullObserver : public IPairFindingObserver {
    // All methods are no-ops
};
```

**Files**:
- `include/x3dna/algorithms/observers/pair_finding_observer.hpp`
- `include/x3dna/algorithms/observers/json_recording_observer.hpp`
- `src/x3dna/algorithms/observers/json_recording_observer.cpp`

**Deliverable**: Observer interface + JSON recording implementation

---

### Task 11.4: Create Residue Index Mapper

**Problem**: Legacy 1-based indices handled ad-hoc throughout; conversion logic duplicated.

**Solution**: Dedicated `ResidueIndexMap` class.

```cpp
/**
 * @class ResidueIndexMap
 * @brief Maps between legacy 1-based and modern 0-based indices
 * 
 * Provides consistent, tested index conversion and residue lookup.
 * Handles structures where legacy indices may have gaps.
 */
class ResidueIndexMap {
public:
    /**
     * @brief Build mapping from structure
     * @param structure Structure to index
     */
    void build(const core::Structure& structure);

    // Lookups
    const core::Residue* get_by_legacy_idx(int legacy_idx) const;
    const core::Residue* get_by_modern_idx(size_t modern_idx) const;

    // Conversions
    int to_legacy(size_t modern_idx) const;
    size_t to_modern(int legacy_idx) const;

    // Iteration (legacy order - matches original algorithm)
    int max_legacy_idx() const { return max_legacy_idx_; }
    
    /**
     * @brief Iterate valid legacy indices in order
     */
    class LegacyIterator {
        // ...
    };
    LegacyIterator begin() const;
    LegacyIterator end() const;

    // Filtering
    std::vector<int> nucleotide_legacy_indices() const;

private:
    std::map<int, const core::Residue*> by_legacy_;
    std::map<size_t, const core::Residue*> by_modern_;
    std::map<int, size_t> legacy_to_modern_;
    std::map<size_t, int> modern_to_legacy_;
    int max_legacy_idx_ = 0;
};
```

**Files**:
- `include/x3dna/algorithms/residue_index_map.hpp`
- `src/x3dna/algorithms/residue_index_map.cpp`

**Deliverable**: ResidueIndexMap with comprehensive tests

---

### Task 11.5: Refactor BasePairFinder as Facade

**Problem**: `BasePairFinder` does everything; should coordinate components.

**Solution**: Refactor to Facade that orchestrates the pipeline.

```cpp
/**
 * @class BasePairFinder
 * @brief Facade that orchestrates pair finding pipeline
 * 
 * Simplified interface that coordinates:
 * - ResidueIndexMap (index management)
 * - PairCandidateCache (validation caching)
 * - IPairSelectionStrategy (selection algorithm)
 * - IPairFindingObserver (event recording)
 */
class BasePairFinder {
public:
    /**
     * @brief Constructor with default strategy
     */
    explicit BasePairFinder(
        const ValidationParameters& params = ValidationParameters::defaults()
    );

    /**
     * @brief Find pairs with observer for recording
     */
    std::vector<core::BasePair> find_pairs(
        core::Structure& structure,
        IPairFindingObserver* observer = nullptr
    );

    // Configuration
    void set_strategy(std::unique_ptr<IPairSelectionStrategy> strategy);
    void set_parameters(const ValidationParameters& params);
    
    // For backward compatibility
    void set_strategy(PairFindingStrategy strategy_enum);
    std::vector<core::BasePair> find_pairs_with_recording(
        core::Structure& structure,
        io::JsonWriter* writer
    );

private:
    BasePairValidator validator_;
    std::unique_ptr<IPairSelectionStrategy> strategy_;
    QualityScoreCalculator quality_calculator_;
};
```

**Implementation Flow**:
```cpp
std::vector<BasePair> BasePairFinder::find_pairs(
    Structure& structure,
    IPairFindingObserver* observer
) {
    // 1. Build index mapping
    ResidueIndexMap index_map;
    index_map.build(structure);

    // 2. Build validation cache (Phase 1)
    PairCandidateCache cache;
    cache.build(structure, validator_);

    // 3. Notify observer of all validations
    if (observer) {
        for (const auto& [pair, info] : cache.valid_pairs()) {
            observer->on_pair_validated(
                pair.first, pair.second,
                info.validation, info.bp_type_id
            );
        }
    }

    // 4. Run selection strategy
    auto selected = strategy_->select(cache, index_map, observer);

    // 5. Convert to BasePair objects
    return create_base_pairs(selected, index_map, observer);
}
```

**Deliverable**: Refactored BasePairFinder under 200 lines

---

### Task 11.6: Extract Quality Score Calculator

**Problem**: `adjust_pair_quality()` and `calculate_bp_type_id()` embedded in BasePairFinder.

**Solution**: Dedicated `QualityScoreCalculator` class.

```cpp
/**
 * @class QualityScoreCalculator
 * @brief Calculates adjusted quality scores for pair selection
 * 
 * Encapsulates:
 * - adjust_pairQuality logic (H-bond based adjustment)
 * - bp_type_id calculation (Watson-Crick/Wobble detection)
 */
class QualityScoreCalculator {
public:
    /**
     * @brief Calculate adjusted quality score for selection
     * @param result Validation result with raw quality score
     * @param hbonds Hydrogen bonds for this pair
     * @param res1 First residue (for bp_type calculation)
     * @param res2 Second residue (for bp_type calculation)
     * @return Adjusted score (lower is better)
     */
    double calculate_selection_score(
        const ValidationResult& result,
        const std::vector<core::hydrogen_bond>& hbonds,
        const core::Residue& res1,
        const core::Residue& res2
    ) const;

    /**
     * @brief Calculate bp_type_id (matches legacy check_wc_wobble_pair)
     * @return -1 (unknown), 0 (invalid), 1 (wobble), 2 (Watson-Crick)
     */
    int calculate_bp_type_id(
        const core::Residue& res1,
        const core::Residue& res2,
        const ValidationResult& result
    ) const;

    /**
     * @brief Calculate H-bond quality adjustment
     */
    double adjust_pair_quality(
        const std::vector<core::hydrogen_bond>& hbonds
    ) const;

private:
    ParameterCalculator param_calculator_;
    
    // Watson-Crick pair list
    static const std::vector<std::string> WC_LIST;
    
    char get_base_letter(core::ResidueType type) const;
};
```

**Files**:
- `include/x3dna/algorithms/quality_score_calculator.hpp`
- `src/x3dna/algorithms/quality_score_calculator.cpp`

**Deliverable**: QualityScoreCalculator with unit tests

---

### Task 11.7: Update Protocol to Use New Architecture

**Problem**: `FindPairProtocol` uses old `BasePairFinder` interface.

**Solution**: Update to use observer pattern for recording.

```cpp
class FindPairProtocol : public ProtocolBase {
public:
    void execute(core::Structure& structure) override {
        // Create observer if recording enabled
        std::unique_ptr<IPairFindingObserver> observer;
        if (json_writer_) {
            observer = std::make_unique<JsonRecordingObserver>(*json_writer_);
        }

        // Calculate frames
        calculate_frames(structure);

        // Find pairs with observer
        base_pairs_ = pair_finder_.find_pairs(structure, observer.get());

        // Detect helices
        if (!single_strand_mode_) {
            detect_helices(structure);
        }
    }

private:
    BasePairFinder pair_finder_;
    // No more json_writer_ passed to find_pairs_with_recording
};
```

**Deliverable**: Updated protocol using observer pattern

---

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                      FindPairProtocol                           │
│  (Orchestrates: frames → validation → selection → helices)      │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│                       BasePairFinder                             │
│  (Facade: coordinates components)                               │
└───────┬───────────┬───────────────┬────────────────┬────────────┘
        │           │               │                │
        ▼           ▼               ▼                ▼
┌───────────┐ ┌───────────────┐ ┌───────────────┐ ┌────────────────┐
│ Residue   │ │ PairCandidate │ │ ISelection    │ │ IFinding       │
│ IndexMap  │ │ Cache         │ │ Strategy      │ │ Observer       │
└───────────┘ └───────────────┘ └───────┬───────┘ └───────┬────────┘
                      │                 │                 │
                      ▼                 ▼                 ▼
              ┌───────────────┐ ┌───────────────┐ ┌────────────────┐
              │ BasePair      │ │ MutualBest    │ │ JsonRecording  │
              │ Validator     │ │ Strategy      │ │ Observer       │
              └───────────────┘ ├───────────────┤ └────────────────┘
                                │ AllPairs      │
                                │ Strategy      │
                                └───────────────┘
```

---

## Testing Plan

### Unit Tests (Each Component)

| Component | Test File | Key Tests |
|-----------|-----------|-----------|
| PairCandidateCache | `test_pair_candidate_cache.cpp` | Build, lookup, iteration |
| ResidueIndexMap | `test_residue_index_map.cpp` | Conversion, gaps, iteration |
| MutualBestPairStrategy | `test_mutual_best_strategy.cpp` | Selection logic, tie-breaking |
| AllPairsStrategy | `test_all_pairs_strategy.cpp` | Exhaustive enumeration |
| QualityScoreCalculator | `test_quality_calculator.cpp` | H-bond adjustment, bp_type_id |
| JsonRecordingObserver | `test_json_observer.cpp` | Event recording |

### Integration Tests

- [ ] Refactored `BasePairFinder` produces identical output to original
- [ ] Strategy switching works correctly
- [ ] Observer recording matches legacy JSON format
- [ ] Performance is comparable (within 5%)

### Regression Tests

- [ ] All 3602 PDB files produce same base pairs
- [ ] JSON output byte-for-byte identical in legacy mode

---

## Success Criteria

- [ ] `BasePairFinder` reduced to < 200 lines
- [ ] Each class has single responsibility
- [ ] No `JsonWriter*` parameters passed through call chains
- [ ] New strategies can be added without modifying existing classes
- [ ] All unit tests pass
- [ ] All regression tests pass
- [ ] No performance regression > 5%

---

## Files Created

```
include/x3dna/algorithms/
├── pair_candidate_cache.hpp
├── residue_index_map.hpp
├── quality_score_calculator.hpp
├── strategies/
│   ├── pair_selection_strategy.hpp
│   ├── mutual_best_pair_strategy.hpp
│   └── all_pairs_strategy.hpp
└── observers/
    ├── pair_finding_observer.hpp
    └── json_recording_observer.hpp

src/x3dna/algorithms/
├── pair_candidate_cache.cpp
├── residue_index_map.cpp
├── quality_score_calculator.cpp
├── strategies/
│   ├── mutual_best_pair_strategy.cpp
│   └── all_pairs_strategy.cpp
└── observers/
    └── json_recording_observer.cpp

tests/unit/algorithms/
├── test_pair_candidate_cache.cpp
├── test_residue_index_map.cpp
├── test_quality_calculator.cpp
├── test_mutual_best_strategy.cpp
└── test_json_observer.cpp
```

---

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Breaking legacy compatibility | **CRITICAL** | Keep original as fallback, extensive regression tests |
| Performance regression | High | Profile before/after, cache aggressively |
| Over-engineering | Medium | Only extract what has clear benefit; YAGNI |
| Integration complexity | Medium | Incremental refactoring; one component at a time |

---

## Migration Path

1. **Phase 1**: Add new classes alongside existing code
2. **Phase 2**: Update `BasePairFinder` to use new components internally
3. **Phase 3**: Expose new interfaces in public API
4. **Phase 4**: Mark old interfaces as deprecated
5. **Phase 5**: Remove old code after verification

---

*Estimated Completion: Week 19-21*

