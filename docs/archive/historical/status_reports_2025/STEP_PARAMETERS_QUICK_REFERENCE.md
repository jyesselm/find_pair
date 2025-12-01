# Step Parameters - Quick Reference Guide

**Date**: 2025-11-29  
**Quick reference for step parameter generation and comparison**

---

## 🚀 Quick Start

### Generate Modern Step Parameters
```bash
# Step 1: Generate input file
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb /tmp/<PDB_ID>.inp

# Step 2: Generate step parameters
./build/analyze_app /tmp/<PDB_ID>.inp

# Output: data/json/bpstep_params/<PDB_ID>.json
#         data/json/helical_params/<PDB_ID>.json
```

### Generate Legacy Step Parameters
```bash
org/build/bin/find_pair_analyze data/pdb/<PDB_ID>.pdb

# Output: data/json_legacy/bpstep_params/<PDB_ID>.json
#         data/json_legacy/helical_params/<PDB_ID>.json
```

### Compare Step Parameters
```bash
# Single PDB
python3 scripts/compare_json.py steps <PDB_ID>

# Test set
python3 scripts/compare_json.py steps --test-set 10

# Verbose output
python3 scripts/compare_json.py steps <PDB_ID> --verbose
```

---

## 📁 File Locations

### Modern JSON Files
```
data/json/
├── bpstep_params/
│   └── <PDB_ID>.json    # Step parameters
└── helical_params/
    └── <PDB_ID>.json    # Helical parameters
```

### Legacy JSON Files
```
data/json_legacy/
├── bpstep_params/
│   └── <PDB_ID>.json    # Step parameters
└── helical_params/
    └── <PDB_ID>.json    # Helical parameters
```

---

## 📊 Step Parameters

### bpstep_params
- **Shift**: Translation along x-axis (Å)
- **Slide**: Translation along y-axis (Å)
- **Rise**: Translation along z-axis (Å)
- **Tilt**: Rotation around x-axis (degrees)
- **Roll**: Rotation around y-axis (degrees)
- **Twist**: Rotation around z-axis (degrees)

### helical_params
- **x_displacement**: X-displacement (Å)
- **y_displacement**: Y-displacement (Å)
- **rise**: Helical rise (Å)
- **inclination**: Inclination angle (degrees)
- **tip**: Tip angle (degrees)
- **twist**: Helical twist (degrees)

---

## 🔍 Verification

### Check if Files Exist
```bash
ls data/json/bpstep_params/<PDB_ID>.json
ls data/json_legacy/bpstep_params/<PDB_ID>.json
```

### Count Records
```bash
python3 -c "import json; print(len(json.load(open('data/json/bpstep_params/<PDB_ID>.json'))))"
```

### Compare Values
```bash
python3 scripts/compare_json.py steps <PDB_ID> --verbose
```

---

## ⚠️ Common Issues

### Issue: Modern analyze calculates 0 step parameters

**Cause**: Input file uses atom indices (from legacy find_pair)

**Solution**: Use modern find_pair to generate input file
```bash
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb /tmp/<PDB_ID>.inp
./build/analyze_app /tmp/<PDB_ID>.inp
```

### Issue: Legacy step parameters not found

**Cause**: Legacy code wasn't generating step parameter JSON files

**Solution**: Already fixed! Legacy now generates step parameters correctly.

### Issue: Different record counts

**Expected**: Legacy may have 2× more records (processes multiple duplexes)
- Legacy: 48 records (2 duplexes × 24 each)
- Modern: 20 records (single set)

**This is normal** - both are correct implementations.

---

## 📝 JSON Format

### Modern Format
```json
{
  "type": "bpstep_params",
  "bp_idx1": 3,
  "bp_idx2": 4,
  "shift": -0.326662,
  "slide": -2.096079,
  "rise": 2.910300,
  "tilt": 3.171240,
  "roll": 11.777801,
  "twist": 28.807861,
  "midstep_frame": {...}
}
```

### Legacy Format
```json
{
  "type": "bpstep_params",
  "bp_idx1": 3,
  "bp_idx2": 4,
  "params": {
    "Shift": -0.326662,
    "Slide": -2.096079,
    "Rise": 2.910300,
    "Tilt": 3.171240,
    "Roll": 11.777801,
    "Twist": 28.807861
  },
  "mst_org": [...],
  "mst_orien": [[...], [...], [...]]
}
```

**Note**: Comparison script handles both formats automatically.

---

## ✅ Verification Checklist

- [ ] Modern step parameters generated
- [ ] Legacy step parameters generated
- [ ] Files exist in correct directories
- [ ] Record counts are reasonable
- [ ] Comparison script runs without errors
- [ ] Values match when base pair indices align

---

## 📚 Related Documentation

- [STEP_PARAMETERS_IMPLEMENTATION.md](STEP_PARAMETERS_IMPLEMENTATION.md) - Detailed implementation
- [STEP_PARAMETERS_STATUS.md](STEP_PARAMETERS_STATUS.md) - Current status
- [STEP_PARAMETERS_DIFFERENCES.md](STEP_PARAMETERS_DIFFERENCES.md) - Modern vs Legacy differences
- [STEP_PARAMETERS_SUMMARY.md](STEP_PARAMETERS_SUMMARY.md) - Complete summary
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing workflows

---

## 🎯 Key Points

1. ✅ **Both implementations work** - Modern and legacy generate step parameters correctly
2. ✅ **Values match exactly** - When base pair indices align, all 6 parameters match
3. ⚠️ **Different counts expected** - Legacy processes multiple duplexes, modern processes once
4. ⚠️ **Different pairs expected** - Legacy and modern find_pair select different pairs
5. ✅ **Calculations are correct** - Differences are in pair selection, not calculations

---

*Quick Reference - 2025-11-29*

