# Step-by-Step: Comparing Residue Reference Frames

**Goal**: Verify that modern code calculates the same reference frames as legacy

**Why This Matters**: Reference frames are the foundation for all geometric calculations. If these don't match, nothing else will.

---

## What is a Reference Frame?

A reference frame is a **local coordinate system** for each nucleotide base, consisting of:

1. **Origin** (x, y, z): A point in 3D space (usually near the base center)
2. **Orientation**: Three orthogonal axes (x-axis, y-axis, z-axis)
   - Stored as a 3x3 rotation matrix
   - Or as 9 values: [x1, x2, x3, y1, y2, y3, z1, z2, z3]

### In Legacy JSON

```json
{
  "residue_idx": 44,
  "chain_id": "A",
  "residue_seq": 44,
  "orien": [0.123, -0.456, 0.789, ...],  // 9 values (rotation matrix)
  "org": [12.34, 56.78, 90.12]            // 3 values (origin)
}
```

### In Modern JSON

```json
{
  "legacy_residue_idx": 44,
  "chain_id": "A", 
  "residue_seq": 44,
  "orien": [0.123, -0.456, 0.789, ...],  // 9 values
  "org": [12.34, 56.78, 90.12]            // 3 values
}
```

---

## Step 1: Pick a Simple Test Case

### Criteria for Good Test PDB

- ✅ Small (< 30 nucleotides)
- ✅ Standard Watson-Crick pairs
- ✅ No modified nucleotides
- ✅ Already validated (in PASS list)

### Recommended: 1EHZ

```bash
# Check if it's validated
grep 1EHZ data/index_validation_status.csv
# Should show: 1EHZ,76,76,76,0,PASS,All indices match perfectly
```

**Why 1EHZ?**
- Well-studied structure
- 76 nucleotides (medium size)
- Standard duplex
- Good reference in documentation

---

## Step 2: Generate Modern JSON for Test PDB

```bash
# Generate with --only-paired to match legacy
./build/generate_modern_json data/pdb/1EHZ.pdb data/json --only-paired

# This creates (among others):
# - data/json/base_frame_calc/1EHZ.json
# - data/json/ls_fitting/1EHZ.json  
# - data/json/frame_calc/1EHZ.json
```

**Check it worked**:
```bash
# Should have same count as legacy
echo "Modern:" && cat data/json/base_frame_calc/1EHZ.json | jq '. | length'
echo "Legacy:" && cat data/json_legacy/base_frame_calc/1EHZ.json | jq '. | length'
```

---

## Step 3: Compare Frame Origins (org)

The origin is simpler than orientation - just 3 numbers (x, y, z).

### Quick Check

```python
import json
import numpy as np

# Load both
with open('data/json/base_frame_calc/1EHZ.json') as f:
    modern = json.load(f)
with open('data/json_legacy/base_frame_calc/1EHZ.json') as f:
    legacy = json.load(f)

# Match by legacy_residue_idx
legacy_by_idx = {r['residue_idx']: r for r in legacy}
modern_by_idx = {r['legacy_residue_idx']: r for r in modern}

# Compare origins
print("Comparing frame origins (org):")
print(f"{'Idx':<6} | {'Chain':<5} | {'Seq':<4} | {'Max Diff':<10} | Status")
print("-" * 50)

for idx in sorted(modern_by_idx.keys())[:10]:  # First 10
    if idx not in legacy_by_idx:
        continue
    
    mod_org = modern_by_idx[idx].get('org', [])
    leg_org = legacy_by_idx[idx].get('org', [])
    
    if len(mod_org) == 3 and len(leg_org) == 3:
        diff = max(abs(mod_org[i] - leg_org[i]) for i in range(3))
        status = "✅" if diff < 0.01 else "⚠️" if diff < 0.1 else "❌"
        
        chain = modern_by_idx[idx].get('chain_id', '?')
        seq = modern_by_idx[idx].get('residue_seq', 0)
        
        print(f"{idx:<6} | {chain:<5} | {seq:<4} | {diff:<10.6f} | {status}")
```

### What to Look For

- **Difference < 0.001 Å**: ✅ Perfect match (floating point precision)
- **Difference < 0.01 Å**: ✅ Excellent match (acceptable)
- **Difference < 0.1 Å**: ⚠️ Small difference (investigate)
- **Difference > 0.1 Å**: ❌ Significant difference (bug likely)

---

## Step 4: Compare Frame Orientations (orien)

Orientation is a 3x3 rotation matrix stored as 9 values.

### Quick Check

```python
# Compare orientations
print("\nComparing frame orientations (orien):")
print(f"{'Idx':<6} | {'Chain':<5} | {'Seq':<4} | {'Max Diff':<10} | Status")
print("-" * 50)

for idx in sorted(modern_by_idx.keys())[:10]:
    if idx not in legacy_by_idx:
        continue
    
    mod_orien = modern_by_idx[idx].get('orien', [])
    leg_orien = legacy_by_idx[idx].get('orien', [])
    
    if len(mod_orien) == 9 and len(leg_orien) == 9:
        diff = max(abs(mod_orien[i] - leg_orien[i]) for i in range(9))
        status = "✅" if diff < 0.0001 else "⚠️" if diff < 0.001 else "❌"
        
        chain = modern_by_idx[idx].get('chain_id', '?')
        seq = modern_by_idx[idx].get('residue_seq', 0)
        
        print(f"{idx:<6} | {chain:<5} | {seq:<4} | {diff:<10.8f} | {status}")
```

### What to Look For

Orientation differences are more sensitive:
- **< 0.0001**: ✅ Perfect match
- **< 0.001**: ✅ Very good (acceptable numerical difference)
- **< 0.01**: ⚠️ Small rotation difference
- **> 0.01**: ❌ Significant rotation difference (investigate)

---

## Step 5: Use Existing Compare Tools

### Option A: Use compare_frames.py

```bash
python3 scripts/compare_frames.py \
    data/json_legacy/base_frame_calc/1EHZ.json \
    data/json/base_frame_calc/1EHZ.json \
    --tolerance 0.001 \
    --show-details
```

### Option B: Use compare_ref_frames.py

If you have legacy ref_frames.dat file:

```bash
# Generate legacy ref_frames.dat
cd org && ./build/bin/find_pair_original ../data/pdb/1EHZ.pdb 1ehz.inp

# This creates ref_frames.dat in org/
python3 scripts/compare_ref_frames.py \
    org/ref_frames.dat \
    ref_frames_modern.dat \
    --by-residue \
    --legacy-inp org/1ehz.inp
```

---

## Step 6: Interpret Results

### If Everything Matches ✅

**Great!** Reference frames are correct. Move to next step:
- Compare base pair parameters
- Compare step parameters

### If Small Differences (< 0.001) ⚠️

**Acceptable!** Due to:
- Floating point precision
- Different compilers
- Numerical method differences

**Action**: Document and accept

### If Large Differences (> 0.01) ❌

**Investigate!** Possible causes:
- Template matching differences
- Atom selection differences  
- Least-squares fitting differences
- Bug in modern or legacy code

**Action**: Debug specific residue

---

## Debugging a Specific Residue

If you find a difference in residue frame:

```bash
# 1. Extract that residue to minimal PDB
python3 scripts/extract_residues.py 1EHZ --residue 44 \
    --output minimal_res44.pdb

# 2. Run legacy with debug
cd org
./build/bin/find_pair_original ../minimal_res44.pdb test.inp

# Check debug output for frame calculation

# 3. Run modern with debug (add debug prints to code)
# 4. Compare step-by-step
```

---

## Quick Start Commands

### Just Do This First

```bash
# 1. Pick test PDB
PDB="1EHZ"

# 2. Generate modern
./build/generate_modern_json data/pdb/${PDB}.pdb data/json --only-paired

# 3. Quick compare (Python one-liner)
python3 << EOF
import json

with open('data/json/base_frame_calc/${PDB}.json') as f:
    modern = json.load(f)
with open('data/json_legacy/base_frame_calc/${PDB}.json') as f:
    legacy = json.load(f)

print(f"Modern: {len(modern)} residues")
print(f"Legacy: {len(legacy)} residues")

# Match and compare
matched = 0
for m in modern:
    m_idx = m.get('legacy_residue_idx', m.get('residue_idx'))
    for l in legacy:
        l_idx = l.get('residue_idx')
        if m_idx == l_idx:
            matched += 1
            # Compare org
            if 'org' in m and 'org' in l:
                diff = max(abs(m['org'][i] - l['org'][i]) for i in range(3))
                if diff > 0.01:
                    print(f"⚠️  Idx {m_idx}: origin diff = {diff:.4f}")
            break

print(f"Matched: {matched}/{len(modern)}")
EOF
```

---

## Expected Timeline

- **Step 1-2**: 2 minutes (setup and generate)
- **Step 3**: 5 minutes (compare origins)
- **Step 4**: 5 minutes (compare orientations)
- **Step 5**: 10 minutes (use comparison tools)
- **Step 6**: Depends on results

**Total**: 20-30 minutes for initial frame comparison

---

## Files to Read

- `scripts/compare_frames.py` - Frame comparison tool
- `scripts/compare_ref_frames.py` - Alternative comparison
- `docs/legacy/ARCHIVED_LEGACY_KNOWLEDGE.md` - Legacy frame calculation details

---

## What Success Looks Like

After this step, you'll know:

✅ Do reference frame **origins** match?  
✅ Do reference frame **orientations** match?  
✅ What's the typical difference (if any)?  
✅ Are there systematic issues or just numerical precision?  

Then you can confidently move to comparing base pair parameters!

---

**Start here**: Generate modern JSON for 1EHZ and compare! 🚀

