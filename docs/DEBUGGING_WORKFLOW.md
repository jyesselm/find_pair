# Complete Debugging Workflow Guide

**Created**: 2025-01-XX  
**Purpose**: Step-by-step guide for using the debugging infrastructure

---

## Quick Start

### Automated Workflow (Recommended)

Use the automated script that does everything:

```bash
./scripts/debug_minimal_case.sh <PDB_ID> [num_pairs]
```

**Example**:
```bash
# Debug a minimal case with 2 consecutive base pairs
./scripts/debug_minimal_case.sh 1AQ4 2

# Debug with 3 consecutive pairs
./scripts/debug_minimal_case.sh 1AQ4 3
```

**What it does**:
1. Extracts minimal test case from PDB
2. Generates legacy JSON output
3. Generates modern JSON output
4. Compares all step-by-step results
5. Shows summary of matches/differences

---

## Manual Workflow (For More Control)

### Step 1: Extract Minimal Test Case

```bash
python3 scripts/extract_minimal_pairs.py \
    --pdb data/pdb/<PDB_ID>.pdb \
    --legacy-json data/json_legacy/find_bestpair_selection/<PDB_ID>.json \
    --output-dir data/pdb/minimal \
    --num-pairs 2 \
    --max-fragments 1
```

**Output**: Creates minimal PDB in `data/pdb/minimal/`

### Step 2: Generate Legacy JSON

```bash
./org/build/bin/find_pair_original data/pdb/minimal/<PDB_ID>_minimal_pairs_1_2.pdb
```

**Output**: Creates JSON files in `data/json_legacy/`:
- `best_partner_candidates/<PDB_ID>_minimal_pairs_1_2.json`
- `mutual_best_decisions/<PDB_ID>_minimal_pairs_1_2.json`
- `iteration_states/<PDB_ID>_minimal_pairs_1_2.json`

### Step 3: Generate Modern JSON

```bash
./build/generate_modern_json \
    data/pdb/minimal/<PDB_ID>_minimal_pairs_1_2.pdb \
    data/json/<PDB_ID>_minimal_pairs_1_2.json \
    --fix-indices
```

**Output**: Creates JSON files in `data/json/<PDB_ID>_minimal_pairs_1_2.json/`

### Step 4: Compare Results

**Compare best partner finding for a specific residue**:
```bash
python3 scripts/compare_best_partner.py <MINIMAL_PDB_ID> <res_i> \
    --legacy-dir data/json_legacy \
    --modern-dir "data/json/<MINIMAL_PDB_ID>.json" \
    --verbose
```

**Compare iteration states**:
```bash
python3 scripts/compare_iteration.py <MINIMAL_PDB_ID> \
    --legacy-dir data/json_legacy \
    --modern-dir "data/json/<MINIMAL_PDB_ID>.json" \
    --verbose
```

**Compare mutual best decisions**:
```bash
python3 scripts/compare_mutual_best.py <MINIMAL_PDB_ID> \
    --legacy-dir data/json_legacy \
    --modern-dir "data/json/<MINIMAL_PDB_ID>.json" \
    --verbose
```

---

## Interpreting Results

### ✅ All Match

If all comparisons show matches:
- ✅ Best partner matches
- ✅ Iteration states match
- ✅ Mutual best decisions match

**Conclusion**: The minimal test case matches perfectly. This confirms the algorithm works correctly for this case.

**Next steps**:
- Try another minimal case (different pair positions)
- Try with more pairs (3, 4, etc.)
- Test with full PDB

### ❌ Differences Found

If differences are found, the comparison tools will show:

**Best Partner Differences**:
- Which candidates differ
- Score differences
- Eligibility differences

**Iteration State Differences**:
- Which iterations differ
- Which pairs found in each iteration
- Matched residue differences

**Mutual Best Decision Differences**:
- Which pairs have different decisions
- Whether pairs are selected differently

---

## Debugging Strategy

### Strategy 1: Find the Divergence Point

1. **Start with iteration states**: 
   - Which iteration first shows differences?
   - This tells you how far into the algorithm the difference occurs

2. **Check mutual best decisions for that iteration**:
   - Which pairs were considered?
   - Which pairs were selected?
   - This narrows down to specific pair decisions

3. **Check best partner candidates for problematic pairs**:
   - Which candidates were considered?
   - Which candidate was selected?
   - Why was it selected (score, eligibility)?

4. **Deep dive into validation**:
   - Use existing validation comparison tools
   - Check quality scores, H-bonds, etc.

### Strategy 2: Progressive Testing

1. **Start small**: Test with 1-2 pairs
2. **Gradually increase**: Add more pairs
3. **Find breaking point**: Identify when differences first appear
4. **Isolate**: Focus on the minimal case that shows the difference

---

## Example: Debugging a Difference

**Scenario**: Full PDB shows different pair selection

```bash
# 1. Extract minimal case
python3 scripts/extract_minimal_pairs.py \
    --pdb data/pdb/DIFFERENT_PDB.pdb \
    --legacy-json data/json_legacy/find_bestpair_selection/DIFFERENT_PDB.json \
    --output-dir data/pdb/minimal \
    --num-pairs 2

# 2. Run automated debugging
./scripts/debug_minimal_case.sh DIFFERENT_PDB 2

# 3. Analyze differences
# If differences found, investigate:
# - Check iteration states to see which iteration differs
# - Check mutual best decisions for that iteration
# - Check best partner candidates for residues in problematic pairs
```

---

## Files and Directories

### Minimal Test Cases
- `data/pdb/minimal/<PDB_ID>_minimal_pairs_<N>_<M>.pdb`
  - Minimal PDB fragments with N consecutive pairs

### Legacy JSON Output
- `data/json_legacy/best_partner_candidates/<PDB_ID>.json`
- `data/json_legacy/mutual_best_decisions/<PDB_ID>.json`
- `data/json_legacy/iteration_states/<PDB_ID>.json`

### Modern JSON Output
- `data/json/<PDB_ID>.json/best_partner_candidates/<PDB_ID>.json`
- `data/json/<PDB_ID>.json/mutual_best_decisions/<PDB_ID>.json`
- `data/json/<PDB_ID>.json/iteration_states/<PDB_ID>.json`

---

## Troubleshooting

### Issue: Modern JSON files not found

**Cause**: Modern code might not have generated them

**Solution**: 
- Ensure `generate_modern_json` was run with `--fix-indices`
- Check that JSON writer is working (look for `[JSON_WRITER]` messages)

### Issue: Comparison shows "Missing in modern"

**Cause**: Modern code might not output that JSON type yet

**Solution**:
- Check that modern code has been built with latest changes
- Verify JSON writer methods are being called

### Issue: Path errors in comparison scripts

**Cause**: Modern JSON uses nested directory structure

**Solution**: 
- Use full path: `data/json/<PDB_ID>.json` as `--modern-dir`
- Scripts should handle this automatically

---

## Tips

1. **Start small**: Always start with 2-pair minimal cases
2. **Use verbose**: Add `--verbose` to see detailed comparisons
3. **Check JSON files**: Look at the actual JSON to understand structure
4. **Iterate**: Fix one difference, then test again
5. **Document**: Keep notes on what you find

---

*This workflow provides a systematic approach to finding and fixing differences between legacy and modern code.*

