# Cleanup Verification Checklist

**Date**: 2025-11-27  
**Purpose**: Verify all cleanup changes work correctly

---

## ✅ File Organization Verification

### Root Directory
- [x] Only essential markdown files remain (README.md + cleanup docs)
- [x] No `.inp` files in root
- [x] No `.out` files in root
- [x] No temporary directories (`temp_compare/`, `base_frame_calc/`)
- [x] No temporary data files (`.pdb`, `.dat`, `.scr`, `.par`)

### Resources Directory
- [x] `resources/templates/` exists with template files
- [x] `resources/config/` exists with config files
- [x] `resources/test_sets/` exists with test set JSON files
- [x] All required runtime resources are accessible

### Data Directory
- [x] Only contains test inputs (`pdb/`) and generated outputs (`json/`, `json_legacy/`)
- [x] No redundant test directories
- [x] Clean structure

### Scripts Directory
- [x] Only core scripts in main directory (`compare_json.py`, `rebuild_json.py`)
- [x] One-off scripts archived in `scripts/archive/`
- [x] Scripts README updated

### Documentation
- [x] Status/session files archived to `docs/archive/status_reports/`
- [x] Comparison reports archived to `docs/archive/comparison_reports/`
- [x] Investigation docs archived to `docs/archive/investigation/`
- [x] Active documentation organized in `docs/`

---

## 🔍 Code Compatibility Verification

### Template Loading
- [ ] **Test**: Code can find templates in `resources/templates/`
- [ ] **Test**: Code falls back to `data/templates/` if needed (backward compatibility)
- [ ] **Verify**: `StandardBaseTemplates` checks `resources/templates/` first

### Test Set Loading
- [ ] **Test**: Scripts can find test sets in `resources/test_sets/`
- [ ] **Test**: Scripts fall back to `data/test_sets/` if needed (backward compatibility)
- [ ] **Verify**: `compare_json.py` and `rebuild_json.py` use new paths

### Config Loading
- [ ] **Test**: Config files accessible from `resources/config/`
- [ ] **Verify**: Any code that loads config uses correct path

---

## 🧪 Build and Test Verification

### Build System
- [ ] **Test**: `make clean && make release` succeeds
- [ ] **Verify**: CMake configuration works with new structure
- [ ] **Check**: No broken paths in build system

### Tests
- [ ] **Test**: `make test` passes
- [ ] **Verify**: Unit tests can find required resources
- [ ] **Check**: Integration tests work correctly

### Runtime
- [ ] **Test**: Application can load templates from `resources/templates/`
- [ ] **Test**: Scripts can access test sets from `resources/test_sets/`
- [ ] **Verify**: All file paths resolve correctly

---

## 📝 Documentation Verification

### Updated References
- [x] `README.md` updated with new directory structure
- [x] `docs/DATA_STRUCTURE.md` updated with resources/ directory
- [x] `scripts/README.md` updated with archive information
- [ ] **Check**: All path references in docs are correct

### Archive Organization
- [x] Status reports archived correctly
- [x] Comparison reports archived correctly
- [x] Investigation docs archived correctly
- [x] PDB analysis docs archived correctly

---

## 🚀 Next Steps (After Verification)

Once verification is complete:

1. **Test build**: `make clean && make release`
2. **Run tests**: `make test`
3. **Test scripts**: Run `compare_json.py` with a test set
4. **Verify templates**: Run a test that uses template loading
5. **Update documentation**: Fix any broken links or paths found during testing

---

## 📊 Cleanup Statistics

- **Files removed**: ~380 temporary files
- **Files archived**: 85+ documentation files
- **Files organized**: 15+ resource files moved
- **Scripts organized**: 34 archived, 2 core kept
- **Root directory**: From 400+ files to 4 markdown files

---

## ✅ Completion Checklist

- [x] Phase 1: Remove temporary files
- [x] Phase 2: Archive status/session files
- [x] Phase 3: Create resources/ directory
- [x] Phase 4: Clean up data/ directory
- [x] Phase 5: Organize scripts directory
- [x] Phase 6: Consolidate documentation
- [x] Phase 7: Update .gitignore
- [x] Phase 8: Update documentation references
- [ ] **Verification**: Build and test verification (TODO)

---

**Note**: Use this checklist to verify all cleanup changes work correctly before committing.

