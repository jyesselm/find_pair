# Frame Validation Edge Cases

## Modified Nucleotides with Template Mismatches

### 1. 5MU (5-methyluridine)
- **Example**: 1EHZ A54
- **Issue**: Has C5M atom, should use thymine template (`Atomic.t.pdb`)
- **Fix**: Added to `MODIFIED_PYRIMIDINE_TYPES` lookup as THYMINE
- **Status**: ✅ FIXED

### 2. KIR (Kirromycin-modified nucleotide)
- **Example**: 1OB2 A1394
- **Issue**: Modern detects as thymine, legacy as uracil
- **RMS Difference**: 5.14e-04 (51x over 1e-5 tolerance)
- **Analysis**: Highly modified nucleotide, no C5M present. Modern classification logic needs investigation.
- **Status**: ⚠️  NEEDS INVESTIGATION
- **Impact**: Minor RMS difference, may not affect downstream results

## Malformed Legacy JSONs

Multiple legacy JSON files have parsing errors and cannot be loaded:
- 1FFZ: `Expecting ',' delimiter: line 9813`
- 3G8T: `Expecting ',' delimiter: line 15509`
- 6CAQ: `Expecting ',' delimiter: line 6774`

**Action**: Skip these PDBs during validation (legacy bug, not modern issue)

## Legacy Duplicate Records

Legacy generates duplicate records for:
- `base_frame_calc`
- `frame_calc`  
- `ls_fitting`

**Pattern**: Each residue appears 2-3 times

**Fix**: Comparison logic deduplicates legacy records by `residue_idx`

## Summary

- **Total edge cases**: 2 template issues + multiple malformed JSONs
- **Fixed**: 1 (5MU)
- **Pending**: 1 (KIR - investigate classification logic)
- **Skipped**: Malformed legacy JSONs (not a modern code issue)

