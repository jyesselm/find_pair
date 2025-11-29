# Protocol Implementation vs Legacy Comparison

**Date**: Current  
**Purpose**: Compare FindPairProtocol implementation with legacy find_pair workflow

## Legacy Workflow (from `org/src/find_pair.c`)

### Main Function: `duplex()`

```c
static void duplex(long num, long num_residue, char *bseq, long **seidx, long *RY,
                   char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                   char **Miscs, double **xyz, struct_args *args, char *parfile,
                   miscPars *misc_pars)
```

### Step-by-Step Legacy Workflow

1. **Initialize Arrays**
   ```c
   orien = dmatrix(1, num_residue, 1, 9);   // Frame rotation matrices
   org = dmatrix(1, num_residue, 1, 3);      // Frame origins
   NC1xyz = dmatrix(1, num_residue, 1, 7);  // N, C1' coordinates
   o3_p = dmatrix(1, num_residue, 1, 8);    // O3' coordinates
   ```

2. **Initialize Water/HTM Handling**
   ```c
   init_htm_water(args->waters, num, num_residue, idx, htm_water);
   identify_htw(...);  // Identify HTM waters
   ```

3. **Calculate Base Frames** (`base_info()`)
   ```c
   base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
             ResSeq, Miscs, xyz, orien, org, NC1xyz, o3_p);
   ```
   - Calculates reference frames for all residues
   - Fills `orien` (rotation matrices) and `org` (origins)
   - Calculates `NC1xyz` and `o3_p` coordinates

4. **Get Ring Atoms** (`ring_oidx()`)
   ```c
   ring_atom = lmatrix(1, num_residue, 1, 19);
   ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);
   ```
   - Identifies ring atoms for each residue
   - Used in pair finding and validation

5. **Find Pairs** (conditional)
   ```c
   if (args->pairs) {
       all_pairs(...);  // Exhaustive search
       goto ALL_PAIRS;
   }
   base_pairs = lmatrix(1, num_residue, 1, nout_p1);
   num_bp = find_bestpair(nout, base_pairs, num_residue, bseq, seidx, RY, AtomName,
                          xyz, idx, orien, org, NC1xyz, ring_atom, misc_pars);
   ```
   - If `-p` flag: use `all_pairs()` (exhaustive)
   - Otherwise: use `find_bestpair()` (greedy mutual best match)

6. **Record JSON** (after find_bestpair)
   ```c
   if (num_bp > 0) {
       json_writer_record_find_bestpair_selection(num_bp, base_pairs);
   }
   ```

7. **Handle No Pairs Case**
   ```c
   if (!num_bp) {
       no_basepairs(args->pdbfile, args->outfile, parfile);
       goto NO_BASE_PAIR;
   }
   ```

8. **Reorder Pairs and Detect Helices** (`re_ordering()`)
   ```c
   bp_idx = lvector(1, num_bp);
   helix_marker = lvector(1, num_bp);
   helix_idx = lmatrix(1, num_bp, 1, 7);
   re_ordering(num_bp, base_pairs, bp_idx, helix_marker, helix_idx, misc_pars,
               &num_helix, o3_p, bseq, seidx, ResName, ChainID, ResSeq, Miscs);
   ```
   - Reorders pairs (5' to 3')
   - Detects helices
   - Marks helix boundaries

9. **Write Output Files**
   ```c
   write_bestpairs(...);  // Write .inp file
   write_helix(...);      // Write helix regions
   // Additional output based on options (curves, etc.)
   ```

## Modern Implementation (FindPairProtocol)

### Current Implementation

```cpp
void FindPairProtocol::execute(core::Structure& structure) {
    // Step 1: Calculate frames for all residues
    calculate_frames(structure);

    // Step 2: Find base pairs
    find_pairs(structure);

    // Step 3: Detect helices (when available)
    // detect_helices(structure);

    // Step 4: Reorder pairs if needed
    // reorder_pairs(structure);
}
```

### Comparison Table

| Legacy Step | Modern Implementation | Status |
|------------|----------------------|--------|
| 1. Initialize arrays | Handled by Structure/Residue classes | ✅ Complete |
| 2. Initialize water/HTM | Not yet implemented | ⏳ Missing |
| 3. Calculate frames (`base_info`) | `calculate_frames()` → `BaseFrameCalculator` | ✅ Complete |
| 4. Get ring atoms (`ring_oidx`) | Handled internally by algorithms | ✅ Complete |
| 5. Find pairs (`find_bestpair`/`all_pairs`) | `find_pairs()` → `BasePairFinder` | ✅ Complete |
| 6. Record JSON | `find_pairs_with_recording()` | ✅ Complete |
| 7. Handle no pairs | Not yet implemented | ⏳ Missing |
| 8. Reorder & detect helices (`re_ordering`) | Placeholders only | ⏳ Missing |
| 9. Write output files | Not in protocol (handled at app level) | ✅ By design |

## Detailed Comparison

### Frame Calculation

**Legacy**:
```c
base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
          ResSeq, Miscs, xyz, orien, org, NC1xyz, o3_p);
```
- Calculates frames for all residues
- Stores in arrays: `orien[1..num_residue][1..9]`, `org[1..num_residue][1..3]`

**Modern**:
```cpp
void FindPairProtocol::calculate_frames(core::Structure& structure) {
    frame_calculator_.set_is_rna(is_rna);
    frame_calculator_.calculate_all_frames(structure);
}
```
- Uses `BaseFrameCalculator::calculate_all_frames()`
- Stores frames on `Residue` objects (OOP design)
- ✅ **Matches legacy behavior**

### Pair Finding

**Legacy**:
```c
if (args->pairs) {
    all_pairs(...);  // Exhaustive
} else {
    num_bp = find_bestpair(...);  // Greedy mutual best
}
```

**Modern**:
```cpp
void FindPairProtocol::find_pairs(core::Structure& structure) {
    if (find_all_pairs_) {
        pair_finder_.set_strategy(PairFindingStrategy::ALL_PAIRS);
    } else {
        pair_finder_.set_strategy(PairFindingStrategy::BEST_PAIR);
    }
    base_pairs_ = pair_finder_.find_pairs(structure);
}
```
- ✅ **Matches legacy behavior**
- Strategy selection matches legacy options

### JSON Recording

**Legacy**:
```c
if (num_bp > 0) {
    json_writer_record_find_bestpair_selection(num_bp, base_pairs);
}
```

**Modern**:
```cpp
if (json_writer_) {
    base_pairs_ = pair_finder_.find_pairs_with_recording(structure, json_writer_);
}
```
- ✅ **Matches legacy behavior**
- Recording happens during pair finding (more efficient)

### Missing Features

#### 1. Water/HTM Handling ⏳
**Legacy**: `init_htm_water()`, `identify_htw()`
**Modern**: Not yet implemented
**Impact**: Low (optional feature)

#### 2. No Pairs Handling ⏳
**Legacy**: `no_basepairs()` writes special output
**Modern**: Not yet implemented
**Impact**: Medium (error handling)

#### 3. Helix Detection & Reordering ⏳
**Legacy**: `re_ordering()` does both
**Modern**: Placeholders only
**Impact**: High (needed for complete workflow)

**Legacy `re_ordering()` does**:
- Reorders pairs (5' to 3')
- Detects helices
- Marks helix boundaries
- Handles circular structures

## Differences in Design

### Data Storage

**Legacy**: Uses 1-based arrays
```c
orien[1..num_residue][1..9]  // 1-based indexing
```

**Modern**: Uses OOP with 0-based indexing
```cpp
residue.reference_frame()  // Stored on Residue object
```

**Impact**: ✅ Handled correctly (conversion at boundaries)

### Iteration Order

**Legacy**: Iterates in PDB file order (1..num_residue)
**Modern**: Can use legacy order via `get_residues_in_legacy_order()`

**Impact**: ✅ Legacy mode support ready

### Error Handling

**Legacy**: Uses `goto` for error handling
**Modern**: Uses exceptions and return values

**Impact**: ✅ Better design, but may need to match legacy behavior

## Recommendations

### High Priority

1. **Implement Helix Detection** (`re_ordering` equivalent)
   - Reorder pairs (5' to 3')
   - Detect helix boundaries
   - Handle circular structures

2. **Implement No Pairs Handling**
   - Match legacy `no_basepairs()` behavior
   - Write appropriate output files

### Medium Priority

3. **Water/HTM Handling**
   - Implement if needed for exact compatibility
   - May be optional depending on use case

### Low Priority

4. **Output File Writing**
   - Currently handled at application level
   - May want to add to protocol for completeness

## Summary

✅ **What Matches**:
- Frame calculation workflow
- Pair finding workflow
- Strategy selection (best pair vs all pairs)
- JSON recording

⏳ **What's Missing**:
- Helix detection and reordering
- No pairs handling
- Water/HTM handling (optional)

🎯 **Overall**: Core workflow matches legacy. Missing features are primarily helix detection and error handling.

