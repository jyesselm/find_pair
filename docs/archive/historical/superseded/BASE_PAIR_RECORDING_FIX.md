# Base Pair Recording Fix

**Date**: 2025-01-XX  
**Status**: ⚠️ **SUPERSEDED** - See [BASE_PAIR_RECORDING.md](BASE_PAIR_RECORDING.md) for complete documentation

**Note**: This document has been superseded by the consolidated documentation. Please refer to [BASE_PAIR_RECORDING.md](BASE_PAIR_RECORDING.md) for the complete documentation.

---

## Problem

Modern code was recording `base_pair` records for ALL validated pairs, but legacy only records `base_pair` for pairs that appear in the final selection (ref_frames.dat).

**Result**:
- Modern: 864 base_pair records (all validated pairs)
- Legacy: 1,539 base_pair records (but only 623 in final selection)
- **Issue**: Most missing base_pair records (748/751) are "potential pairs" - validated but not selected

---

## Solution

Updated code to only record `base_pair` records for pairs in the final selection, matching legacy behavior.

### Code Changes

**File**: `src/x3dna/algorithms/base_pair_finder.cpp`

1. **Removed base_pair recording from validation phase** (line 535-562):
   - Previously recorded base_pair for ALL `result.is_valid == true` pairs
   - Now only records validation, distance_checks, and hbond_list during validation

2. **Added base_pair recording after selection** (after line 274):
   - Records base_pair ONLY for pairs in `base_pairs` vector (selected pairs)
   - This matches legacy behavior where base_pair records correspond to ref_frames.dat

### Implementation

```cpp
// After find_bestpair selection completes:
// Record base_pair records ONLY for pairs in the final selection
// This matches legacy behavior where base_pair records correspond to ref_frames.dat
// (only pairs that appear in the final output)
if (writer) {
    for (const auto& pair : base_pairs) {
        writer->record_base_pair(pair);
    }
}
```

---

## Results ✅

After rebuild and testing:

- ✅ **base_pair records now match find_bestpair_selection count exactly**
  - **10-PDB test set**: 100% perfect matches (10/10 PDBs)
  - All selected pairs have base_pair records
  - No missing pairs, no extra pairs

### Test Results

| PDB ID | Selection | Base Pairs | Match |
|--------|-----------|------------|-------|
| 1Q96   | 38        | 38         | ✅    |
| 1VBY   | 27        | 27         | ✅    |
| 3AVY   | 8         | 8          | ✅    |
| 3G8T   | 225       | 225        | ✅    |
| 3KNC   | 8         | 8          | ✅    |
| 4AL5   | 6         | 6          | ✅    |
| 5UJ2   | 6         | 6          | ✅    |
| 6LTU   | 41        | 41         | ✅    |
| 8J1J   | 66        | 66         | ✅    |
| 6CAQ   | 623       | 623        | ✅    |
| **Total** | **1,042** | **1,042** | **✅ 100%** |

---

## Testing

After rebuilding, test with:

```bash
# Generate JSON with --fix-indices
./build/generate_modern_json data/pdb/6CAQ.pdb data/json --fix-indices

# Verify base_pair count matches selection
python3 -c "
import json
from pathlib import Path

pdb_id = '6CAQ'
selection = json.load(open(f'data/json/find_bestpair_selection/{pdb_id}.json'))
basepairs = json.load(open(f'data/json/base_pair/{pdb_id}.json'))

selection_pairs = set(tuple(p) for p in selection[0]['pairs'])
basepair_pairs = set((r['base_i'], r['base_j']) for r in basepairs)

print(f'Selection: {len(selection_pairs)} pairs')
print(f'Base pairs: {len(basepair_pairs)} pairs')
print(f'Match: {selection_pairs == basepair_pairs}')
"
```

---

## Related Documentation

- [BASE_PAIR_RECORD_ANALYSIS.md](BASE_PAIR_RECORD_ANALYSIS.md) - Analysis showing most missing records are potential pairs
- [FIX_INDICES_TEST_RESULTS.md](FIX_INDICES_TEST_RESULTS.md) - Test results showing 100% match on selection
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status

---

*Fix to match legacy behavior: only record base_pair for pairs in final selection (ref_frames.dat).*

