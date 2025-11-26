# bp_type_id Implementation Requirements

## Summary

The `bp_type_id` calculation is currently simplified and requires full implementation from Stage 6 (base pair step parameters).

## Current Implementation (Simplified)

**Location**: `src/x3dna/algorithms/base_pair_finder.cpp::calculate_bp_type_id`

**Current Logic**:
1. Initializes `bp_type_id = -1` (matches legacy)
2. Checks direction vector conditions: `dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0`
3. If conditions met:
   - Checks if base pair string matches WC_LIST (AT, TA, AU, UA, GC, CG, IC, CI)
   - Sets `bp_type_id = 2` if Watson-Crick (in WC_LIST)
   - Sets `bp_type_id = 1` if Wobble (BasePairType::WOBBLE)
4. **FIXED**: Preserves `bp_type_id = -1` for valid pairs (legacy behavior)
5. Only sets to `0` for invalid pairs

**Missing**: Full `check_wc_wobble_pair` logic with geometric parameters

## Legacy Implementation (Full)

**Location**: `org/src/cmn_fncs.c::check_wc_wobble_pair`

**Legacy Logic**:
```c
void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch, double opening)
{
    if (fabs(stretch) > 2.0 || fabs(opening) > 60)
        return;
    if (dval_in_range(fabs(shear), 1.8, 2.8))
        *bpid = 1;  // Wobble
    if (fabs(shear) <= 1.8 && num_strmatch(bp, WC, 1, 8))
        *bpid = 2;  // Watson-Crick
}
```

**Called from**: `calculate_more_bppars` (line ~4618 in `cmn_fncs.c`)
```c
if (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0) {
    check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
    // pars[1] = shear, pars[2] = stretch, pars[6] = opening
    if (*bpid == 2)
        rtn_val[5] -= 2.0;  // Quality score adjustment
}
```

## Required Parameters

### Base Pair Step Parameters (bpstep_par)
The full implementation requires `bpstep_par` which provides:
- **shear** (`pars[1]`): Horizontal displacement
- **stretch** (`pars[2]`): Vertical displacement  
- **opening** (`pars[6]`): Rotation angle

These parameters are calculated in Stage 6 (ParameterCalculator) and are not yet available.

## Impact on Quality Score

When `bp_type_id == 2` (Watson-Crick):
- Legacy subtracts `2.0` from `rtn_val[5]` (quality_score)
- Modern also subtracts `2.0` from quality_score (line 361 in `base_pair_finder.cpp`)

**Problem**: Without full `check_wc_wobble_pair` logic:
- Some pairs get `bp_type_id = 0` instead of `2`
- Missing the `-2.0` quality_score adjustment
- Affects pair selection when quality_scores are close

## Implementation Plan

### Step 1: Wait for Stage 6
- Implement `ParameterCalculator` to calculate `bpstep_par`
- Extract `shear`, `stretch`, and `opening` from `bpstep_par`

### Step 2: Implement Full `check_wc_wobble_pair`
```cpp
void check_wc_wobble_pair(int& bp_type_id, const std::string& bp_str,
                         double shear, double stretch, double opening) {
    // Check stretch and opening limits
    if (std::abs(stretch) > 2.0 || std::abs(opening) > 60.0) {
        return;  // Keep current bp_type_id
    }
    
    // Check for Wobble: fabs(shear) in [1.8, 2.8]
    if (std::abs(shear) >= 1.8 && std::abs(shear) <= 2.8) {
        bp_type_id = 1;  // Wobble
    }
    
    // Check for Watson-Crick: fabs(shear) <= 1.8 && base pair in WC_LIST
    if (std::abs(shear) <= 1.8) {
        // WC_LIST: "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
        if (bp_str == "AT" || bp_str == "TA" || bp_str == "AU" || bp_str == "UA" ||
            bp_str == "GC" || bp_str == "CG" || bp_str == "IC" || bp_str == "CI") {
            bp_type_id = 2;  // Watson-Crick
        }
    }
}
```

### Step 3: Update `calculate_bp_type_id`
```cpp
int BasePairFinder::calculate_bp_type_id(
    const Residue* res1, const Residue* res2,
    const ValidationResult& result,
    double /* quality_score */,
    const BasePairStepParameters& bpstep_par  // NEW: from Stage 6
) const {
    int bp_type_id = -1;
    
    // Check direction vector conditions
    if (result.dir_x > 0.0 && result.dir_y < 0.0 && result.dir_z < 0.0) {
        // Get base pair string
        char base1 = res1->one_letter_code();
        char base2 = res2->one_letter_code();
        std::string bp_str = std::string(1, std::toupper(base1)) + 
                            std::string(1, std::toupper(base2));
        
        // Get geometric parameters from bpstep_par
        double shear = bpstep_par.shear();
        double stretch = bpstep_par.stretch();
        double opening = bpstep_par.opening();
        
        // Call full check_wc_wobble_pair logic
        check_wc_wobble_pair(bp_type_id, bp_str, shear, stretch, opening);
    }
    
    // Default to 0 if still -1
    if (bp_type_id == -1 && result.is_valid) {
        bp_type_id = 0;
    } else if (!result.is_valid) {
        bp_type_id = 0;
    }
    
    return bp_type_id;
}
```

## Blocking Dependencies

1. **Stage 6: ParameterCalculator**
   - Must implement `BasePairStepParameters` calculation
   - Must calculate `shear`, `stretch`, and `opening` for each base pair
   - Must integrate with `BasePairFinder`

2. **Base Pair Step Parameters Class**
   - Already exists: `include/x3dna/core/parameters.hpp::BasePairStepParameters`
   - Need to ensure it has `shear()`, `stretch()`, and `opening()` accessors

## Current Workaround

The simplified version works for many cases but may miss some Watson-Crick pairs that should get `bp_type_id = 2` and the `-2.0` quality score adjustment.

## Test Cases

Once Stage 6 is implemented, test with:
- Pairs that should be `bp_type_id = 2` (Watson-Crick with `fabs(shear) <= 1.8`)
- Pairs that should be `bp_type_id = 1` (Wobble with `fabs(shear) in [1.8, 2.8]`)
- Pairs that should be `bp_type_id = 0` (other cases)

Verify that quality score adjustments match legacy exactly.

