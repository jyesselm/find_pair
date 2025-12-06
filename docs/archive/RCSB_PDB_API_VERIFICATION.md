# RCSB PDB API Verification - Lessons Learned

**Date**: December 5, 2025  
**Issue**: Incorrectly classified several modified nucleotides as "legacy bugs"  
**Solution**: Used RCSB PDB GraphQL API to verify actual base types  

---

## The Problem

I incorrectly assumed that nucleotides with C7 substituents (like CM0) were thymine derivatives because thymine has a C5 methyl group. This led to misclassifying several nucleotides and claiming legacy code had bugs when it was actually correct!

---

## The Correct Approach

**Always verify using the authoritative source: RCSB PDB Chemical Component Dictionary**

### GraphQL API Query:
```graphql
query molecule ($id: String!) {
    chem_comp(comp_id:$id){
        chem_comp {
            id
            name
            formula
            type
        }
        rcsb_chem_comp_synonyms {
          name
          type
        }
    }
}
```

### How to Determine Base Type:
Look at the chemical name for key terms:
- **URIDINE** → URACIL base → code: `u`
- **CYTIDINE** → CYTOSINE base → code: `c`
- **GUANOSINE** → GUANINE base → code: `g`
- **THYMIDINE** → THYMINE base → code: `t`
- **ADENOSINE** → ADENINE base → code: `a`
- **INOSINE** → INOSINE base → code: `i`

---

## Verification Results

### CM0 - I WAS WRONG ❌
- **API Result**: "5-(CARBOXYMETHOXY) **URIDINE**-5'-MONOPHOSPHATE"
- **Actual Base**: URACIL (not Thymine!)
- **My Error**: Saw C7 group, assumed thymine
- **Truth**: The modification is on position 5, base is still uracil
- **Correction**: Changed from THYMINE to URACIL in registry
- **Legacy Status**: Legacy was CORRECT! ✅

### JSP - NOT A STANDARD NUCLEOTIDE ❌
- **API Result**: "(1R)-1-(4-amino-1-methyl-2-oxo-1,2-dihydropyrimidin-5-yl)..."
- **Actual Base**: Unknown/Complex pyrimidine derivative
- **My Error**: Assumed it was a standard nucleotide type
- **Truth**: Complex modified pyrimidine, should use auto-detection
- **Correction**: Removed from registry
- **Legacy Status**: Let auto-detection handle it

### NCA - NOT EVEN A NUCLEOTIDE! ❌❌❌
- **API Result**: "NICOTINAMIDE"
- **Actual Compound**: Vitamin B3 (niacin)
- **My Error**: Completely misidentified as a nucleotide
- **Truth**: It's a vitamin, not a nucleotide at all!
- **Correction**: Removed from registry
- **Legacy Status**: Should never have been added

### 9DG - I WAS CORRECT ✅
- **API Result**: "9-DEAZAGUANINE"
- **Actual Base**: GUANINE
- **My Classification**: GUANINE → Correct!
- **Legacy**: Used Uracil template ❌
- **Modern**: Uses Guanine template ✅
- **Status**: This IS a real legacy bug - modern is more correct

---

## Template Assignment Rules (IMPORTANT!)

### For Modified Nucleotides:
- Lowercase one-letter code → `Atomic.x.pdb` (e.g., `Atomic.u.pdb`)
- Uppercase one-letter code → `Atomic_X.pdb` (e.g., `Atomic_U.pdb`)

### Examples:
- Standard Uracil (`U`): `Atomic_U.pdb`
- Modified Uracil (`u`): `Atomic.u.pdb`
- CM0 (modified uracil): code=`u` → `Atomic.u.pdb` ✅

---

## Lessons Learned

### ❌ DON'T:
1. **Don't assume** chemical structure from substituent groups
2. **Don't guess** nucleotide types based on atom names
3. **Don't call something a "legacy bug"** without verification
4. **Don't add to registry** without checking PDB database

### ✅ DO:
1. **Always use RCSB PDB API** to verify chemical identity
2. **Look for -URIDINE, -CYTIDINE, etc.** in the chemical name
3. **Trust the PDB database** - it's the authoritative source
4. **Verify legacy behavior** before claiming it's wrong
5. **Use lowercase codes** for modified nucleotides

---

## Corrected Registry

### Before (WRONG):
```json
"CM0": {"code": "t", "type": "THYMINE", ...},  // WRONG!
"JSP": {"code": "t", "type": "THYMINE", ...},  // WRONG!
"NCA": {"code": "t", "type": "THYMINE", ...}   // WRONG!
```

### After (CORRECT):
```json
"CM0": {"code": "u", "type": "URACIL", ...},  // ✅ Fixed
// JSP removed - not standard nucleotide
// NCA removed - not a nucleotide
"9DG": {"code": "g", "type": "GUANINE", ...}  // ✅ Was already correct
```

---

## Impact

### Registry Changes:
- **85 → 83 nucleotides** (removed 2 non-nucleotides)
- **1 corrected** (CM0: THYMINE → URACIL)
- **1 confirmed** (9DG: correctly GUANINE)

### "Legacy Bugs" Status:
- **Before**: Claimed 4 legacy bugs (CM0, JSP, NCA, 9DG)
- **After**: Only 1 real legacy bug (9DG)
- **Accuracy**: 75% of my "bug" claims were wrong!

### Template Matching:
- **CM0**: Now matches legacy ✅
- **9DG**: Modern is more correct than legacy ✅
- **JSP/NCA**: Removed from registry, use auto-detection

---

## Tool Created

**File**: `scripts/check_nucleotide_types.py`

Queries RCSB PDB GraphQL API to verify:
- Chemical name
- Formula  
- Base type classification
- One-letter code assignment

**Usage**:
```bash
python3 scripts/check_nucleotide_types.py
```

This tool should be used to verify ANY modified nucleotide before adding to the registry!

---

## Summary

**Key Takeaway**: The RCSB PDB Chemical Component Dictionary is the **ONLY** authoritative source for nucleotide classification. Chemical structure analysis (looking at atoms, substituents, etc.) can be misleading.

**Result**: 
- Fixed 1 misclassification (CM0)
- Removed 2 non-nucleotides (JSP, NCA)  
- Confirmed 1 real legacy bug (9DG)
- Created verification tool for future additions

**Lesson**: Trust the database, not assumptions! 🎓

---

**Created**: December 5, 2025  
**Author**: AI Assistant (learning from mistakes!)  
**Teacher**: User (who taught me to use RCSB PDB API properly)

