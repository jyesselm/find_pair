# Cleanup Plan - Current State

**Date**: December 2, 2025  
**Status**: Many files still need cleanup

---

## Current Issues

**Root directory has**:
- 20+ markdown files (should be ~3-4)
- Test .inp files (1AQ4.inp, temp_7EH2.inp, test_*.inp)
- Parameter .par files
- temp/ directory
- ref_frames_modern.dat

---

## Cleanup Actions

### 1. Archive Validation Documents

**Move to** `docs/validation/`:

```bash
mkdir -p docs/validation

# Validation docs
mv VALIDATION_COMPLETE.md docs/validation/
mv VALIDATION_TODO.md docs/validation/
mv VALIDATION_IN_PROGRESS.md docs/validation/
mv VALIDATION_BUGS_FOUND.md docs/validation/
mv WEEK1_VALIDATION_COMPLETE.md docs/validation/
mv INDEX_MATCHING_PROOF.md docs/validation/

# Investigation docs
mv INVESTIGATION_7EH2.md docs/validation/
mv RESOLUTION_7EH2.md docs/validation/
```

### 2. Archive Old Planning Docs

**Move to** `docs/archive/planning/`:

```bash
mkdir -p docs/archive/planning

mv CLEANUP_AND_RESTRUCTURE_PLAN.md docs/archive/planning/
mv CLEANUP_FINAL_SUMMARY.md docs/archive/planning/
mv CLEANUP_VERIFICATION.md docs/archive/planning/
mv README_CLEANUP.md docs/archive/planning/
mv REPOSITORY_CLEANUP_PLAN.md docs/archive/planning/
mv IMPLEMENTATION_PLAN.md docs/archive/planning/
```

### 3. Archive Session/Status Docs

**Move to** `docs/archive/sessions/`:

```bash
mkdir -p docs/archive/sessions

mv SESSION_SUMMARY.md docs/archive/sessions/
mv NEXT_SESSION_QUICK_START.md docs/archive/sessions/
mv START_HERE.md docs/archive/sessions/
mv IMMEDIATE_ACTIONS.md docs/archive/sessions/
mv BEFORE_AFTER_COMPARISON.md docs/archive/sessions/
```

### 4. Keep Essential Root Files

**KEEP in root** (active/current):
- `README.md` - Main project README
- `NEXT_STEPS.md` - What to do next
- `RUN_FULL_VALIDATION.sh` - Active script
- Build files (Makefile, CMakeLists.txt, etc.)

### 5. Remove Temporary Files

```bash
# Test input files
rm -f 1AQ4.inp
rm -f temp_7EH2.inp  
rm -f test_1TTT.inp
rm -f test_9CF3.inp

# Parameter files (should be in data/ if needed)
rm -f auxiliary.par bp_helical.par bp_step.par cf_7methods.par

# Temporary data files
rm -f ref_frames_modern.dat

# Temp directory
rm -rf temp/
```

### 6. Update .gitignore

```bash
# Add to .gitignore
echo "" >> .gitignore
echo "# Temporary files" >> .gitignore
echo "*.inp" >> .gitignore
echo "*.out" >> .gitignore
echo "*.par" >> .gitignore
echo "temp/" >> .gitignore
echo "ref_frames_*.dat" >> .gitignore
```

---

## Execute Cleanup

Run these commands:

```bash
cd /Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/cpp/find_pair_2

# Create directories
mkdir -p docs/validation
mkdir -p docs/archive/planning
mkdir -p docs/archive/sessions

# Move validation docs
mv VALIDATION_COMPLETE.md INDEX_MATCHING_PROOF.md \
   INVESTIGATION_7EH2.md RESOLUTION_7EH2.md \
   VALIDATION_TODO.md VALIDATION_IN_PROGRESS.md \
   VALIDATION_BUGS_FOUND.md WEEK1_VALIDATION_COMPLETE.md \
   docs/validation/

# Move planning docs
mv CLEANUP_*.md REPOSITORY_CLEANUP_PLAN.md \
   README_CLEANUP.md IMPLEMENTATION_PLAN.md \
   docs/archive/planning/

# Move session docs
mv SESSION_SUMMARY.md NEXT_SESSION_QUICK_START.md \
   START_HERE.md IMMEDIATE_ACTIONS.md \
   BEFORE_AFTER_COMPARISON.md \
   docs/archive/sessions/

# Remove temporary files
rm -f *.inp *.par ref_frames_modern.dat
rm -rf temp/

# Update .gitignore
cat >> .gitignore << 'EOF'

# Temporary test files
*.inp
*.out
*.par
temp/
ref_frames_*.dat
EOF

# Commit cleanup
git add -A
git commit -m "Clean up root directory - archive validation and planning docs"
```

---

## Result

After cleanup, root should have:
- `README.md`
- `NEXT_STEPS.md`
- `Makefile`, `CMakeLists.txt`, `pyproject.toml`
- `.gitignore`
- Build/config files only

All documentation organized in `docs/`:
- `docs/validation/` - Validation documents
- `docs/archive/planning/` - Planning documents
- `docs/archive/sessions/` - Session notes
- `docs/` - Active documentation

---

## Run Cleanup Now

Just execute:

```bash
./scripts/cleanup_root.sh  # If script exists
# Or run the commands above manually
```

