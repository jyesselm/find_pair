# X3DNA Find_Pair Cleanup and Restructure Plan

**Created**: December 2, 2025  
**Goal**: Massive reduction in files, centralized comparison framework, fix legacy indices matching

---

## Overview

This plan addresses the core issues preventing progress:

1. **Documentation bloat** - Too many outdated docs making it hard to understand the current state
2. **Tool proliferation** - Too many scripts and cpp files doing similar things
3. **Index matching problem** - The fundamental blocker that keeps recurring
4. **Comparison fragmentation** - No unified way to track what's working and what's not

---

## Phase 1: Documentation Cleanup

### Files to KEEP

**Core Documentation (6 files)**:
- `README.md` - Main entry point (will be streamlined)
- `docs/JSON_DATA_TYPES_AND_COMPARISONS.md` - JSON format and comparison methodology
- `docs/legacy/10_JSON_STRUCTURE.md` - Legacy JSON format
- `docs/legacy/04_ALGORITHMS.md` - How legacy code works
- `docs/CODE_FLOW.md` - How modern code works (to be updated)
- `docs/TESTING_GUIDE.md` - Testing workflow (to be simplified)

### Files to ARCHIVE/DELETE

**Remove Build Plans** (already built, not needed):
- `docs/BUILD_INSTRUCTIONS.md`
- `docs/QUICK_START.md`
- `docs/modernization/` (entire directory - 17 files)

**Remove Old Debugging Docs** (starting over):
- `docs/DEBUG_1TTT_RESIDUE_16.md`
- `docs/DEBUG_9CF3_RESIDUE_27.md`
- `docs/DEBUGGING_TOOLS.md`
- `docs/DEBUGGING_WORKFLOW.md`
- `docs/REF_FRAMES_*` (5 files)
- `docs/STEP_BY_STEP_*` (3 files)
- `docs/RMSD_RESIDUE_RECOGNITION.md`
- `docs/MATCHING_PLAN*.md` (2 files)
- `docs/LEGACY_ORDER_FIXES.md`
- `docs/COMPARISON_FIXES_SUMMARY.md`
- `docs/STEP_PARAMETERS_IMPLEMENTATION.md`

**Remove Redundant Documentation**:
- `docs/ALGORITHM_CRITICAL_GUIDE.md` (merge critical info into CODE_FLOW.md)
- `docs/COMPARISON_INDEX_RULES.md` (outdated)
- `docs/DATA_STRUCTURE.md` (simplify into README)
- `docs/FIX_INDICES_OPTION.md` (will be obsolete after fix)
- `docs/LEGACY_INDICES_GUIDE.md` (will be obsolete after fix)
- `docs/PROJECT_SUMMARY.md` (redundant with README)

**Archive Directory** (already archived):
- `docs/archive/` - Keep as is (already archived)

**Total Reduction**: ~35 docs → ~6 core docs

---

## Phase 2: Tool Consolidation

### Current Problem

**36 C++ tools in `tools/`** doing overlapping tasks:
- 15 comparison tools (`compare_*.cpp`)
- 8 debug tools (`debug_*.cpp`)
- 5 test tools (`test_*.cpp`)
- 8 other tools (generate, analyze, investigate, etc.)

**57 Python scripts in `scripts/`** doing overlapping tasks:
- 20+ archived scripts
- 15+ comparison scripts
- 10+ investigation scripts
- Cluster tools

### New Unified Framework

**Keep ONLY these tools**:

#### Core Production Tools (4 tools):
1. **`tools/generate_modern_json.cpp`** - Generate JSON from modern code
2. **`tools/generate_legacy_json.cpp`** (rename from org executable) - Generate JSON from legacy code
3. **`scripts/rebuild_json.py`** - Manage JSON regeneration
4. **`scripts/unified_compare.py`** (NEW) - Single comparison framework

#### Archive Everything Else:
- Move all other tools to `tools/archive/`
- Move archived scripts to `scripts/archive/` (already exists)
- Delete redundant scripts from `scripts/`

### New Unified Comparison Framework

**`scripts/unified_compare.py`** - Single point of entry for all comparisons:

```python
# Compare JSON data
python3 scripts/unified_compare.py json 1H4S --verbose

# Compare program outputs (.inp files, ref_frames.dat)
python3 scripts/unified_compare.py output 1H4S --verbose

# Compare specific record types
python3 scripts/unified_compare.py json 1H4S --record-type base_pair

# Batch comparison with CSV tracking
python3 scripts/unified_compare.py batch --test-set 100 --output-csv data/comparison_status.csv

# Resume from last run (uses CSV)
python3 scripts/unified_compare.py batch --test-set 100 --resume
```

**Features**:
- Single entry point for all comparisons
- CSV tracking of PDB status at each stage
- No caching by default (explicit `--cache` flag if needed)
- Clear progress tracking
- Unified reporting format

---

## Phase 3: Fix Legacy Index Matching (THE CRITICAL FIX)

### The Core Problem

**Current Issue**: Modern and legacy code assign different indices to residues, making it impossible to compare results accurately.

**Root Cause**: We don't track the order residues are read in and which ones get filtered out.

### The Solution: Residue Tracking Vector

**Implementation Plan**:

#### 1. Create `ResidueTracker` class in modern code:

```cpp
class ResidueTracker {
public:
    struct ResidueRecord {
        int read_index;           // Order read from PDB (0-based)
        int legacy_index;         // Index from legacy JSON (1-based)
        int modern_index;         // Index after filtering (0-based)
        bool filtered;            // Was it filtered out?
        std::string filter_reason;  // Why was it filtered?
        
        // PDB properties for matching
        std::string chain_id;
        int residue_seq;
        std::string insertion;
        std::string residue_name;
    };
    
    std::vector<ResidueRecord> residues;
    
    void add_residue(/* PDB properties */);
    void mark_filtered(int read_index, const std::string& reason);
    void assign_modern_index(int read_index, int modern_index);
    void load_legacy_indices(const std::string& legacy_json_path);
    void validate_match() const;  // Verify we match legacy exactly
    void export_mapping(const std::string& output_path) const;
};
```

#### 2. Integrate into PDB parsing:

```cpp
// In PDBParser::parse()
ResidueTracker tracker;

// As we read each residue from PDB
for (each residue in PDB) {
    tracker.add_residue(chain, seq, ins, name);
    // ... existing parsing logic ...
}

// After filtering
for (each residue that gets filtered) {
    tracker.mark_filtered(read_index, "no ring atoms");
    // or "failed frame calculation", etc.
}

// After final residue list is created
for (int i = 0; i < final_residues.size(); i++) {
    tracker.assign_modern_index(residue.read_index, i);
}

// Load legacy indices from JSON
tracker.load_legacy_indices("data/json_legacy/base_frame_calc/1H4S.json");

// VALIDATE - this must pass for all PDBs
tracker.validate_match();

// Export mapping for debugging
tracker.export_mapping("data/index_mapping/1H4S.json");
```

#### 3. Validation Requirements:

**Before ANY comparison work continues**, we must:

1. Generate index mappings for ALL test PDBs
2. Validate 100% match for all test PDBs
3. Document any PDBs that don't match and why
4. Fix any discrepancies before proceeding

**CSV Tracking File**: `data/index_validation_status.csv`

```csv
pdb_id,num_residues_read,num_legacy,num_modern,num_filtered,match_status,notes
1H4S,42,40,40,2,PASS,2 residues filtered (no ring atoms)
2BNA,24,24,24,0,PASS,All residues match
3G8T,1200,1180,1185,20,FAIL,Modern has 5 extra residues - investigate
```

---

## Phase 4: Unified Comparison Framework Implementation

### CSV Status Tracking

**File**: `data/comparison_status.csv`

Tracks progress through comparison stages:

```csv
pdb_id,index_match,atoms,frames,hbonds,pairs,step_params,last_updated,notes
1H4S,PASS,PASS,PASS,PASS,PASS,PASS,2025-12-02,100% match
2BNA,PASS,PASS,PASS,FAIL,PENDING,PENDING,2025-12-02,Hbond mismatch - 3 bonds different
3G8T,FAIL,PENDING,PENDING,PENDING,PENDING,PENDING,2025-12-02,Index mismatch - 5 extra residues
6CAQ,PASS,PASS,PASS,PASS,FAIL,PENDING,2025-12-01,4 pair mismatches
```

**Stages**:
1. `index_match` - Residue indices match legacy
2. `atoms` - PDB atom parsing matches
3. `frames` - Reference frame calculation matches
4. `hbonds` - Hydrogen bond detection matches
5. `pairs` - Base pair selection matches
6. `step_params` - Step parameters match

**Status Values**: `PASS`, `FAIL`, `PENDING`, `SKIP`

### Comparison Workflow

**New unified workflow**:

```bash
# Step 1: Validate indices for all PDBs
python3 scripts/unified_compare.py validate-indices --test-set 100

# Step 2: Compare stage by stage
python3 scripts/unified_compare.py batch --test-set 100 --stage atoms
python3 scripts/unified_compare.py batch --test-set 100 --stage frames
python3 scripts/unified_compare.py batch --test-set 100 --stage hbonds
python3 scripts/unified_compare.py batch --test-set 100 --stage pairs

# Step 3: Get summary
python3 scripts/unified_compare.py status

# Step 4: Investigate failures
python3 scripts/unified_compare.py investigate 3G8T --stage index_match
```

### No More Tool Proliferation

**RULE**: Before creating any new script or C++ tool:

1. Ask: "Can this be done with `unified_compare.py`?"
2. If yes → Add a flag/command to unified_compare.py
3. If no → Document why in a comment before creating

---

## Phase 5: Caching Strategy

### Current Problem
Caching has been confusing and may have caused false positives/negatives.

### New Strategy

**Default: NO CACHING**
- Every comparison regenerates JSON unless explicitly told otherwise
- This is slower but guarantees correctness

**Explicit Caching**:
```bash
# Enable caching only when explicitly requested
python3 scripts/unified_compare.py batch --test-set 100 --cache

# Clear cache
python3 scripts/unified_compare.py clear-cache
```

**Cache Location**: `data/.cache/`
- JSON cache: `data/.cache/json/`
- Result cache: `data/.cache/results/`

**Cache Invalidation**: Cache is invalidated when:
- Modern code is recompiled (check build timestamp)
- Legacy code is recompiled (check build timestamp)
- PDB file is modified

---

## Phase 6: File Structure After Cleanup

```
find_pair_2/
├── README.md                          # Streamlined main README
├── Makefile                           # Build rules
├── CMakeLists.txt                     # CMake config
│
├── docs/                              # MINIMAL documentation
│   ├── README.md                      # Docs navigation
│   ├── JSON_DATA_TYPES_AND_COMPARISONS.md  # JSON format reference
│   ├── CODE_FLOW.md                   # How modern code works
│   ├── TESTING_GUIDE.md               # Simple testing workflow
│   ├── legacy/                        # Legacy code docs
│   │   ├── 04_ALGORITHMS.md           # Legacy algorithms
│   │   └── 10_JSON_STRUCTURE.md       # Legacy JSON format
│   └── archive/                       # Archived docs (kept for history)
│
├── src/x3dna/                         # Modern C++ implementation
├── include/x3dna/                     # Modern C++ headers
├── org/                               # Legacy code (reference)
│
├── tools/                             # MINIMAL tools
│   ├── generate_modern_json.cpp       # Generate modern JSON
│   └── archive/                       # All other tools (archived)
│
├── scripts/                           # MINIMAL scripts
│   ├── unified_compare.py             # SINGLE comparison framework
│   ├── rebuild_json.py                # JSON regeneration
│   ├── download_pdbs.py               # PDB download utility
│   └── archive/                       # All other scripts (archived)
│
├── x3dna_json_compare/                # Python comparison library (keep as is)
│
├── data/                              # Data files
│   ├── pdb/                           # PDB files
│   ├── json/                          # Modern JSON output
│   ├── json_legacy/                   # Legacy JSON output (reference)
│   ├── index_mapping/                 # Residue index mappings (NEW)
│   ├── comparison_status.csv          # Unified status tracking (NEW)
│   ├── index_validation_status.csv    # Index validation tracking (NEW)
│   └── .cache/                        # Cache directory (NEW)
│
├── tests/                             # Test suite
└── build/                             # Build artifacts
```

**File Count Reduction**:
- Docs: ~50 files → ~7 files (86% reduction)
- Tools: 36 cpp files → 1 cpp file + archive (97% reduction)  
- Scripts: ~60 py files → 3 py files + archive (95% reduction)

---

## Phase 7: Implementation Order

### Week 1: Index Matching Foundation
1. ✅ Create `ResidueTracker` class
2. ✅ Integrate into PDB parsing
3. ✅ Add legacy index loading
4. ✅ Validate on test_set_10
5. ✅ Fix any discrepancies
6. ✅ Generate index mappings for test_set_100
7. ✅ Document any PDBs that don't match

**Success Criteria**: 100% index match validation for all test PDBs

### Week 2: Unified Comparison Framework
1. ✅ Create `scripts/unified_compare.py` skeleton
2. ✅ Implement CSV status tracking
3. ✅ Implement stage-by-stage comparison
4. ✅ Implement batch processing
5. ✅ Test on test_set_10
6. ✅ Run on test_set_100

**Success Criteria**: Can compare all test PDBs with clear status tracking

### Week 3: Tool Consolidation
1. ✅ Archive old C++ tools
2. ✅ Archive old Python scripts
3. ✅ Update Makefile to only build essential tools
4. ✅ Test that build still works
5. ✅ Update CMakeLists.txt

**Success Criteria**: Clean tool directory, everything still works

### Week 4: Documentation Cleanup
1. ✅ Archive old documentation
2. ✅ Streamline README.md
3. ✅ Update CODE_FLOW.md with current architecture
4. ✅ Simplify TESTING_GUIDE.md
5. ✅ Create docs/README.md navigation

**Success Criteria**: Clear, minimal documentation

### Week 5: Full System Validation
1. ✅ Run full comparison on test_set_100
2. ✅ Investigate ALL failures
3. ✅ Fix root causes
4. ✅ Achieve 100% match on test_set_100
5. ✅ Document remaining edge cases

**Success Criteria**: Clear path to 100% match on all test PDBs

---

## Critical Rules Going Forward

### 1. Index Matching is Sacred
- **NEVER** compare results without validated index matching
- **ALWAYS** check `index_validation_status.csv` first
- If indices don't match, **STOP** and fix before proceeding

### 2. No Tool Proliferation
- **ONE** comparison framework: `unified_compare.py`
- **ONE** JSON generator: `generate_modern_json.cpp`
- **ONE** status tracker: `comparison_status.csv`
- Before creating new tool: Add to unified framework

### 3. CSV Tracking is Truth
- `comparison_status.csv` is the single source of truth
- Update it after every comparison run
- Use it to resume work
- Use it to track progress

### 4. Caching is Explicit
- Default: NO caching
- Enable only with `--cache` flag
- Clear cache when in doubt

### 5. Documentation is Minimal
- Only keep docs that are actively used
- Archive everything else
- Update docs when code changes
- Delete redundant docs immediately

---

## Success Metrics

### Immediate (Week 1)
- ✅ 100% index match validation for test_set_10
- ✅ Index mapping JSON files for all test PDBs
- ✅ Clear documentation of any mismatches

### Short Term (Week 3)
- ✅ Unified comparison framework working
- ✅ CSV status tracking implemented
- ✅ Tool count reduced by 90%
- ✅ Doc count reduced by 80%

### Medium Term (Week 5)
- ✅ 100% match on atoms for test_set_100
- ✅ 100% match on frames for test_set_100
- ✅ 100% match on hbonds for test_set_100
- ✅ Clear understanding of remaining pair mismatches

### Long Term (Week 8)
- ✅ 100% match on ALL stages for test_set_100
- ✅ Validated on test_set_1000
- ✅ Production-ready modern code

---

## Questions to Answer Before Starting

1. ✅ Should we delete archived files or just move them?
   - **Answer**: Move to archive/ directories for history

2. ✅ Should we keep the cluster comparison tools?
   - **Answer**: Yes, but integrate into unified framework

3. ✅ Should we keep the x3dna_json_compare Python module?
   - **Answer**: Yes, it's the comparison library used by unified_compare.py

4. ✅ What's the minimum viable documentation set?
   - **Answer**: 6-7 core docs covering JSON format, code flow, testing

5. ✅ Should index validation block ALL other work?
   - **Answer**: YES - this is the foundation, must be 100% correct

---

## Next Immediate Steps

**Stop everything and focus on index matching:**

1. Implement `ResidueTracker` class
2. Integrate into modern code PDB parsing
3. Validate on one PDB (1H4S)
4. Validate on test_set_10
5. Only proceed when 100% validated

**No more:**
- Creating new comparison scripts
- Creating new debug tools
- Creating new documentation
- Trying to fix pair mismatches without index validation

**Start fresh with:**
- Solid index matching foundation
- Unified comparison framework
- Clear CSV status tracking
- Minimal, focused toolset

---

## End Goal

A clean, maintainable codebase where:
- ✅ We know EXACTLY which residues match between legacy and modern
- ✅ We have ONE tool for all comparisons
- ✅ We track progress clearly with CSV files
- ✅ We have minimal, essential documentation
- ✅ We can confidently iterate on fixes
- ✅ We achieve 100% match with legacy code

**The foundation must be solid before we can build on it.**

