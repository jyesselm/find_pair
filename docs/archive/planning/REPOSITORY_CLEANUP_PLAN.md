# Repository Cleanup and Reorganization Plan

**Created**: 2025-01-XX  
**Status**: Planning Phase  
**Goal**: Massive cleanup and reorganization to improve repository structure

---

## Executive Summary

This repository contains significant clutter:
- **374+ `.inp` files** in root directory (temporary test inputs)
- **5 `.out` files** in root directory (temporary outputs)
- **60+ markdown status/session files** in root (many duplicates/redundant)
- **50+ analysis/investigation scripts** in `scripts/` (many one-off)
- **Multiple redundant documentation files**
- **Temporary test directories** that should be cleaned
- **Duplicate JSON data** in various locations

This plan will systematically clean and reorganize everything.

---

## Phase 1: Remove Temporary Files

### 1.1 Remove `.inp` and `.out` Files from Root

**Action**: Delete all `.inp` and `.out` files from repository root.

**Rationale**: These are temporary test input/output files that can be regenerated. They should not be version controlled.

**Files to remove**:
- All `*_legacy.inp` files (374+ files)
- All `*_modern.inp` files  
- All `*.out` files (5 files: `1A34.out`, `6V9Q.out`, `7EH2.out`, etc.)

**Commands**:
```bash
# Remove all .inp files from root
find . -maxdepth 1 -name "*.inp" -type f -delete

# Remove all .out files from root
find . -maxdepth 1 -name "*.out" -type f -delete
```

**Impact**: Removes ~380 files from root, significantly cleaning up the directory.

### 1.2 Remove Temporary Directories

**Action**: Delete temporary test directories.

**Directories to remove**:
- `temp_compare/` - Contains temporary comparison files
- `base_frame_calc/` - Temporary test directory (already in .gitignore)

**Commands**:
```bash
rm -rf temp_compare/
rm -rf base_frame_calc/
```

### 1.3 Remove Temporary Data Files

**Action**: Remove temporary PDB, DAT, SCR, and obsolete JSON files from root.

**Files to remove**:
- `bestpairs.pdb` - Temporary PDB output
- `bp_order.dat` - Temporary data file
- `hel_regions.pdb` - Temporary PDB output
- `col_chains.scr` - Temporary script output
- `col_helices.scr` - Temporary script output
- `comparison_results_summary.json` - Can be regenerated
- `*.par` files in root (parameter files that should be in data/ if needed)

**Commands**:
```bash
cd /Users/jyesselman2/Dropbox/2_code/cpp/find_pair_2
rm -f bestpairs.pdb bp_order.dat hel_regions.pdb
rm -f col_chains.scr col_helices.scr
rm -f comparison_results_summary.json
# Review .par files - move to data/ if needed, else delete
```

---

## Phase 2: Consolidate Root Markdown Files

### 2.1 Identify Essential Root Files (Keep)

**Essential files to KEEP in root**:
1. `README.md` - Main project README
2. `Makefile` - Build system
3. `CMakeLists.txt` - Build configuration
4. `pyproject.toml` - Python package configuration
5. `.gitignore` - Git ignore rules
6. `.cursorrules` - Cursor configuration (if exists)

### 2.2 Archive Status/Session Files

**Action**: Move all status, session, and implementation summary files to `docs/archive/status_reports/`.

**Files to archive** (60+ files):
- `PROJECT_COMPLETE.md`
- `PROJECT_STATUS.md`
- `MODERNIZATION_STATUS.md`
- `IMPLEMENTATION_STATUS.md`
- `PROTOCOLS_STATUS.md`
- `PROTOCOLS_COMPLETE.md`
- `PROTOCOLS_FINAL_STATUS.md`
- `PROTOCOLS_IMPLEMENTATION_STATUS.md`
- `PROTOCOLS_IMPLEMENTATION_COMPLETE.md`
- `PROTOCOLS_READY_TO_COMMIT.md`
- `SESSION_COMPLETE_SUMMARY.md`
- `SESSION_SUMMARY.md`
- `COMPLETE_SESSION_SUMMARY.md`
- `IMPLEMENTATION_SESSION_SUMMARY.md`
- `IMPLEMENTATION_SUMMARY.md`
- `IMPLEMENTATION_ROADMAP.md`
- `NEXT_IMPLEMENTATION_STEPS.md`
- `WHAT_NEXT.md`
- `NEXT_ACTIONS.md`
- `RESTART_GUIDE.md`
- `COMMIT_GUIDE.md`
- `DOCUMENTATION_INDEX.md`
- `CLEANUP_COMPLETE.md`
- `CLEANUP_PLAN_FILES.md`
- All other `*_STATUS.md`, `*_SUMMARY.md`, `*_SESSION.md` files

**Commands**:
```bash
mkdir -p docs/archive/status_reports

# Move all status/session files
mv PROJECT_*.md docs/archive/status_reports/ 2>/dev/null
mv MODERNIZATION_STATUS.md docs/archive/status_reports/ 2>/dev/null
mv IMPLEMENTATION*.md docs/archive/status_reports/ 2>/dev/null
mv PROTOCOLS*.md docs/archive/status_reports/ 2>/dev/null
mv SESSION*.md docs/archive/status_reports/ 2>/dev/null
mv COMPLETE_SESSION*.md docs/archive/status_reports/ 2>/dev/null
mv NEXT*.md docs/archive/status_reports/ 2>/dev/null
mv WHAT_NEXT.md docs/archive/status_reports/ 2>/dev/null
mv RESTART_GUIDE.md docs/archive/status_reports/ 2>/dev/null
mv COMMIT_GUIDE.md docs/archive/status_reports/ 2>/dev/null
mv DOCUMENTATION_INDEX.md docs/archive/status_reports/ 2>/dev/null
mv CLEANUP*.md docs/archive/status_reports/ 2>/dev/null
```

### 2.3 Consolidate Comparison/Investigation Files

**Action**: Move comparison and investigation files to appropriate archive locations.

**Files to archive**:
- Comparison results → `docs/archive/comparison_reports/`:
  - `BATCH_COMPARISON_RESULTS.md`
  - `COMPARISON_SUMMARY.md`
  - `FINAL_100_PDB_COMPARISON_RESULTS.md`
  - `OUTPUT_COMPARISON.md`
  - `OUTPUT_COMPARISON_NO_JSON.md`
  - `LARGE_SWEEP_RESULTS.md`

- Investigation files → `docs/archive/investigation/`:
  - `INVESTIGATION_COMPLETE_SUMMARY.md`
  - `INVESTIGATION_ROOT_CAUSE.md`
  - `DIFFERENCE_INVESTIGATION_PLAN.md`
  - `PARAMETER_COMPARISON_INVESTIGATION.md`
  - `PARAMETER_DIFFERENCES_ANALYSIS.md`
  - `PARAMETER_MATCHING_SUMMARY.md`
  - `PARAMETER_VERIFICATION_SUMMARY.md`
  - `REFERENCE_FRAME_COMPARISON.md`
  - `HELICAL_PARAMS_VERIFICATION.md`
  - `MODERN_MISSES_ANALYSIS.md`
  - `3CME_EDGE_CASE_ANALYSIS.md`
  - `TIE_BREAKING_FIX_SUMMARY.md`
  - `TIE_BREAKING_FIX_FINAL_RESULTS.md`

**Commands**:
```bash
# Ensure archive directories exist
mkdir -p docs/archive/comparison_reports
mkdir -p docs/archive/investigation

# Move comparison files
mv BATCH_COMPARISON_RESULTS.md docs/archive/comparison_reports/
mv COMPARISON_SUMMARY.md docs/archive/comparison_reports/
mv FINAL_100_PDB_COMPARISON_RESULTS.md docs/archive/comparison_reports/
mv OUTPUT_COMPARISON*.md docs/archive/comparison_reports/
mv LARGE_SWEEP_RESULTS.md docs/archive/comparison_reports/

# Move investigation files
mv INVESTIGATION*.md docs/archive/investigation/
mv DIFFERENCE_INVESTIGATION_PLAN.md docs/archive/investigation/
mv PARAMETER*.md docs/archive/investigation/
mv REFERENCE_FRAME_COMPARISON.md docs/archive/investigation/
mv HELICAL_PARAMS_VERIFICATION.md docs/archive/investigation/
mv MODERN_MISSES_ANALYSIS.md docs/archive/investigation/
mv 3CME_EDGE_CASE_ANALYSIS.md docs/archive/investigation/
mv TIE_BREAKING*.md docs/archive/investigation/
```

---

## Phase 3: Reorganize Resources Directory

### 3.1 Create Resources Directory Structure

**Action**: Create a `resources/` directory at repository root to house all required runtime resources.

**Rationale**: 
- Separates required runtime resources from generated/test data
- Makes it clear what files are needed for the application to run
- Follows standard C++ project organization practices
- Makes it easier to package/distribute the application

**New structure**:
```
resources/
├── templates/          # Template files (Atomic_*.pdb)
├── config/             # Configuration files
└── test_sets/          # Test set definitions
```

### 3.2 Move Required Runtime Resources

**Action**: Move required runtime resources from `data/` to `resources/`.

**Resources to move** (these are REQUIRED at runtime):

1. **Templates** → `resources/templates/`:
   - `data/templates/Atomic_*.pdb` (all template files)
   - These are loaded by `StandardBaseTemplates` class at runtime
   - Required for base frame calculation

2. **Configuration** → `resources/config/`:
   - `data/config/special_residues.json`
   - Application configuration file

3. **Test Sets** → `resources/test_sets/`:
   - `data/test_sets/*.json` (test_set_10.json, test_set_100.json, etc.)
   - Test set definitions (needed for testing but not runtime)

**Commands**:
```bash
# Create resources directory structure
mkdir -p resources/templates
mkdir -p resources/config
mkdir -p resources/test_sets

# Move templates
mv data/templates/*.pdb resources/templates/

# Move config
mv data/config/*.json resources/config/

# Move test sets
mv data/test_sets/*.json resources/test_sets/

# Remove now-empty directories
rmdir data/templates
rmdir data/config
rmdir data/test_sets
```

### 3.3 Update Code to Use Resources Directory

**Action**: Update code to look for resources in new location.

**Files to update**:

1. **`src/x3dna/algorithms/standard_base_templates.cpp`**:
   - Update `find_template_path()` to check `resources/templates/` first
   - Keep backward compatibility with `data/templates/` if exists
   - Update default path logic

2. **Configuration loading** (if any code loads `special_residues.json`):
   - Update paths to use `resources/config/`
   - Keep backward compatibility with `data/config/` if exists

3. **Test set loading** (in scripts):
   - Update script paths to use `resources/test_sets/`
   - Keep backward compatibility with `data/test_sets/` if exists

**Approach**: 
- Prefer `resources/` directory
- Fall back to `data/` directory for backward compatibility
- Use environment variable or build-time path if needed

**Code changes needed**:
```cpp
// In standard_base_templates.cpp, update find_template_path():
std::filesystem::path StandardBaseTemplates::find_template_path() {
    // Check X3DNA_HOMEDIR environment variable first
    const char* x3dna_home = std::getenv("X3DNA_HOMEDIR");
    if (x3dna_home) {
        std::filesystem::path config_path = std::filesystem::path(x3dna_home) / "config";
        if (std::filesystem::exists(config_path)) {
            return config_path;
        }
    }

    // Check resources/templates directory (new location)
    std::filesystem::path resources_path = "resources/templates";
    if (std::filesystem::exists(resources_path)) {
        return resources_path;
    }

    // Fall back to data/templates directory (backward compatibility)
    std::filesystem::path local_path = "data/templates";
    if (std::filesystem::exists(local_path)) {
        return local_path;
    }

    // Return empty path if not found
    return std::filesystem::path();
}
```

### 3.4 Organize Parameter Files

**Action**: Review `.par` files in root and remove (these are generated output).

**Files to remove**:
- `auxiliary.par` - Generated output file
- `bp_helical.par` - Generated output file  
- `bp_step.par` - Generated output file
- `cf_7methods.par` - Generated output file

**Rationale**: These are generated output files from the legacy `analyze` program, not input resources. They should not be version controlled.

**Commands**:
```bash
rm -f auxiliary.par bp_helical.par bp_step.par cf_7methods.par
```

---

## Phase 4: Clean Up Data Directory

### 4.1 Remove Redundant JSON Directories

**Action**: Remove test/duplicate JSON directories.

**Directories to remove**:
- `data/json_test/` - Test JSON files (redundant with main json/ directory)
- `data/residue_ordering_legacy/` - Legacy test data (can be regenerated)
- `data/json/6CAQ.json/` - Duplicate directory structure (if exists)
- Any other duplicate/temporary JSON directories

**Commands**:
```bash
rm -rf data/json_test/
rm -rf data/residue_ordering_legacy/
# Check for and remove any duplicate structures in data/json/
```

### 4.2 Update Data Directory Structure

**Action**: Clarify that `data/` is now only for generated/test data.

**New `data/` structure**:
```
data/
├── pdb/               # Test PDB files (test inputs, not required resources)
├── json/              # Generated modern JSON output (regenerable)
└── json_legacy/       # Generated legacy JSON output (regenerable, for comparison)
```

**Rationale**: 
- `data/` contains test inputs and generated outputs
- All required runtime resources are now in `resources/`
- Clear separation between resources and test/generated data

---

## Phase 5: Organize Scripts Directory

### 5.1 Identify Core Scripts (Keep)

**Core scripts to KEEP in `scripts/`**:
1. `compare_json.py` - Main JSON comparison tool
2. `rebuild_json.py` - JSON regeneration tool
3. `README.md` - Scripts documentation

### 5.2 Archive One-Off Analysis Scripts

**Action**: Move one-off analysis/investigation scripts to `scripts/archive/`.

**Scripts to archive** (40+ files):
- All `analyze_*.py` scripts
- All `investigate_*.py` scripts
- All `compare_*.py` scripts (except `compare_json.py`)
- All `debug_*.py` scripts
- All `create_*.py` scripts
- All `generate_*.py` scripts (except if they're core tools)
- All `*.sh` scripts (shell scripts)

**Commands**:
```bash
mkdir -p scripts/archive

# Move one-off scripts
mv scripts/analyze_*.py scripts/archive/ 2>/dev/null
mv scripts/investigate_*.py scripts/archive/ 2>/dev/null
mv scripts/compare_*.py scripts/archive/ 2>/dev/null
# Keep compare_json.py in main directory
mv scripts/compare_json.py scripts/temp_compare_json.py
mv scripts/archive/compare_json.py scripts/ 2>/dev/null || true
mv scripts/temp_compare_json.py scripts/compare_json.py

mv scripts/debug_*.py scripts/archive/ 2>/dev/null
mv scripts/create_*.py scripts/archive/ 2>/dev/null
mv scripts/generate_*.py scripts/archive/ 2>/dev/null
mv scripts/*.sh scripts/archive/ 2>/dev/null

# Keep core tools
# (compare_json.py, rebuild_json.py, README.md already handled)
```

### 5.3 Update Script Paths for Resources

**Action**: Update scripts that reference `data/test_sets/` to use `resources/test_sets/`.

**Scripts to update**:
- Any scripts that load test set JSON files
- Comparison scripts that reference test sets
- Scripts in `x3dna_json_compare/` package

**Approach**: 
- Update paths to check `resources/test_sets/` first
- Keep backward compatibility with `data/test_sets/` if exists

### 5.4 Create Scripts Index

**Action**: Create `scripts/README.md` that documents:
- Which scripts are current/active
- Which scripts are archived
- How to use the main scripts

---

## Phase 6: Consolidate Documentation

### 6.1 Review docs/ Directory Structure

**Current structure**:
```
docs/
├── archive/              # Already has some archives
│   ├── comparison_reports/
│   ├── investigation/
│   └── pdb_analysis/
├── legacy/              # Keep - Legacy code docs
├── modernization/       # Keep - Stage-by-stage plan
└── [many .md files]     # Need to organize
```

### 6.2 Organize Active Documentation

**Action**: Consolidate active docs and archive outdated ones.

**Active docs to keep in `docs/` root**:
- `README.md` - Documentation index
- `TESTING_GUIDE.md` - Testing guide (primary)
- `BUILD_INSTRUCTIONS.md` - Build instructions
- `QUICK_START.md` - Quick start guide
- `SIMPLIFIED_DEVELOPMENT.md` - Development guide
- `JSON_DATA_TYPES_AND_COMPARISONS.md` - JSON reference
- `ALGORITHM_CRITICAL_GUIDE.md` - Algorithm reference
- `DEBUGGING_STRATEGIES.md` - Debugging guide
- `MODERNIZATION_PLAN.md` - Modernization plan
- `DATA_STRUCTURE.md` - Data directory structure

**Docs to archive to `docs/archive/`**:
- All `*_STATUS.md` files (already handled in Phase 2)
- All `*_SUMMARY.md` files (already handled)
- Investigation-specific docs (already handled)
- Older version docs that are superseded

---

## Phase 7: Update .gitignore

### 7.1 Add Patterns to .gitignore

**Action**: Ensure all temporary/test files are ignored.

**Additions needed**:
```gitignore
# Temporary test files
*.inp
*.out

# Temporary directories
temp_compare/
base_frame_calc/

# Temporary data files
bestpairs.pdb
bp_order.dat
hel_regions.pdb
col_*.scr
*.par

# Generated comparison results
comparison_results_summary.json
```

---

## Phase 8: Create Final Structure

### 8.1 Target Root Directory Structure

```
find_pair_2/
├── README.md                    # Main README (keep)
├── Makefile                     # Build system (keep)
├── CMakeLists.txt               # Build config (keep)
├── pyproject.toml               # Python config (keep)
├── .gitignore                   # Git ignore (keep, updated)
├── .cursorrules                 # Cursor config (if exists, keep)
│
├── include/                     # Headers (keep)
├── src/                         # Source code (keep)
├── org/                         # Legacy code (keep)
├── tests/                       # Tests (keep)
├── tools/                       # Tools (keep)
├── build/                       # Build output (.gitignored)
│
├── resources/                   # Required runtime resources (NEW)
│   ├── templates/               # Template files (Atomic_*.pdb)
│   ├── config/                  # Configuration files
│   └── test_sets/               # Test set definitions
│
├── data/                        # Test data and generated outputs
│   ├── pdb/                     # Test PDB files (test inputs)
│   ├── json/                    # Generated modern JSON
│   └── json_legacy/             # Generated legacy JSON
│
├── scripts/                     # Scripts (organized)
│   ├── README.md                # Scripts documentation
│   ├── compare_json.py          # Main comparison tool
│   ├── rebuild_json.py          # JSON regeneration tool
│   └── archive/                 # Archived scripts
│
└── docs/                        # Documentation (organized)
    ├── README.md                # Docs index
    ├── TESTING_GUIDE.md         # Testing guide
    ├── BUILD_INSTRUCTIONS.md    # Build instructions
    ├── QUICK_START.md           # Quick start
    ├── [other active docs]      # Active documentation
    ├── legacy/                  # Legacy code docs
    ├── modernization/           # Modernization plan
    └── archive/                 # Archived docs
        ├── status_reports/      # Status/session files
        ├── comparison_reports/  # Comparison results
        ├── investigation/       # Investigation docs
        └── pdb_analysis/        # PDB analysis docs
```

### 8.2 Update Documentation References

**Action**: Update documentation that references old `data/` paths to use new `resources/` paths.

**Files to update**:
- `README.md` - Update data directory structure section
- `docs/DATA_STRUCTURE.md` - Update to reflect new structure
- `docs/BUILD_INSTRUCTIONS.md` - Update any path references
- `docs/TESTING_GUIDE.md` - Update test set references

**Changes needed**:
- Replace `data/templates/` → `resources/templates/`
- Replace `data/config/` → `resources/config/`
- Replace `data/test_sets/` → `resources/test_sets/`
- Keep `data/pdb/`, `data/json/`, `data/json_legacy/` references (these stay in data/)

### 8.3 Create Cleanup Summary Document

**Action**: Create `CLEANUP_EXECUTION_LOG.md` to document what was cleaned and why.

---

## Execution Order

1. **Phase 1**: Remove temporary files (safest, immediate impact)
2. **Phase 2**: Consolidate root markdown files (major cleanup)
3. **Phase 3**: Reorganize resources directory (move required resources)
4. **Phase 4**: Clean up data directory (organize test/generated data)
5. **Phase 5**: Organize scripts (improve usability)
6. **Phase 6**: Consolidate documentation (final polish)
7. **Phase 7**: Update .gitignore (prevent future clutter)
8. **Phase 8**: Create final structure summary and update docs

---

## Estimated Impact

### Files to Remove
- **~380** `.inp` files from root
- **~5** `.out` files from root
- **~60** redundant markdown files from root
- **~40** one-off scripts
- **Multiple** temporary directories

### Files to Archive
- **~60** status/session markdown files
- **~40** analysis scripts
- **~20** investigation/comparison files

### Net Result
- **Root directory**: From 400+ files to ~6 essential files
- **Clear organization**: Everything in appropriate directories
- **Improved maintainability**: Easy to find what you need
- **Better .gitignore**: Prevents future clutter

---

## Safety Considerations

### Before Starting
1. **Create a backup branch**: `git checkout -b backup-before-cleanup`
2. **Commit current state**: Ensure all work is saved
3. **Test git status**: Understand what will be changed

### During Execution
1. **One phase at a time**: Execute phases sequentially
2. **Test after each phase**: Verify nothing critical was removed
3. **Document changes**: Update CLEANUP_EXECUTION_LOG.md

### After Execution
1. **Verify build still works**: `make clean && make release`
2. **Verify tests still pass**: `make test`
3. **Review git status**: Check what will be committed
4. **Update README.md**: Reflect new structure if needed

---

## Rollback Plan

If something goes wrong:
1. **Revert git changes**: `git reset --hard backup-before-cleanup`
2. **Restore from backup**: If git revert doesn't work
3. **Selective restore**: Use git to restore specific files if needed

---

## Next Steps

1. **Review this plan** with team/stakeholders
2. **Create backup branch** (`git checkout -b backup-before-cleanup`)
3. **Execute Phase 1** (safest, most impactful)
4. **Review results** and proceed with remaining phases
5. **Update documentation** to reflect new structure

---

## Questions to Resolve

Before executing, clarify:
1. ~~Are `.par` files needed?~~ → **RESOLVED**: They are generated output, should be removed
2. Are any `.inp` files actually needed for reference?
3. Should archived scripts be removed entirely or kept in archive?
4. Are there any critical files not accounted for?
5. Should we update CMake to install resources directory? (For packaging/distribution)
6. Should `data/pdb/` files be kept or removed? (These are test inputs, not required resources)

## Resources Directory Rationale

### Why `resources/` instead of `data/`?

1. **Clear separation**: 
   - `resources/` = Required runtime resources (templates, config)
   - `data/` = Test inputs and generated outputs

2. **Standard practice**: 
   - Many C++ projects use `resources/` for required files
   - Makes it clear what's needed to run vs what's test/generated data

3. **Packaging/distribution**: 
   - Easy to identify what needs to be included in distribution
   - Can be installed to standard locations (e.g., `/usr/share/x3dna/resources/`)

4. **Backward compatibility**: 
   - Code will check `resources/` first, then fall back to `data/` if needed
   - Allows gradual migration

