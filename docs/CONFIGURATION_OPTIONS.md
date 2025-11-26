# Configuration Options

## BaseFrameCalculator Configuration

The `BaseFrameCalculator` class supports several configuration options for controlling frame calculation behavior.

### 1. RNA Detection (`set_is_rna()`)

Controls whether the structure is treated as RNA, which includes C1' atom in ring atom matching.

**Default:** `false` (DNA mode)

**Usage:**
```cpp
BaseFrameCalculator calculator("data/templates");
calculator.set_is_rna(true);  // Enable RNA mode
```

**Effect:**
- **RNA mode (true):** Includes `" C1'"` as the first atom in ring atom matching
  - Purines: `[" C1'", " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "]` (10 atoms)
  - Pyrimidines: `[" C1'", " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "]` (7 atoms)
- **DNA mode (false):** Does not include C1'
  - Purines: `[" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "]` (9 atoms)
  - Pyrimidines: `[" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "]` (6 atoms)

### 2. Legacy Compatibility Mode (`set_legacy_mode()`)

Controls whether to exclude C4 atom from ring atom matching, matching the behavior of legacy code.

**Default:** `false` (correct behavior, includes C4)

**Usage:**
```cpp
BaseFrameCalculator calculator("data/templates");
calculator.set_legacy_mode(true);  // Enable legacy compatibility mode
```

**Effect:**
- **Legacy mode (true):** Excludes `" C4 "` from ring atom matching
  - Purines: `[" N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "]` (8 atoms)
  - Pyrimidines: `[" N3 ", " C2 ", " N1 ", " C6 ", " C5 "]` (5 atoms)
- **Normal mode (false):** Includes C4 when present (correct behavior)
  - Purines: `[" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "]` (9 atoms)
  - Pyrimidines: `[" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "]` (6 atoms)

**Note:** Legacy mode replicates a bug in the legacy code where C4 atom was systematically excluded from matched atoms, even though it was present in both the PDB and standard templates. This mode is provided for compatibility with legacy results but is **not recommended** for new analyses.

### 3. Combined Configuration

Both options can be combined:

```cpp
BaseFrameCalculator calculator("data/templates");
calculator.set_is_rna(true);        // RNA mode
calculator.set_legacy_mode(true);   // Legacy compatibility mode
```

**Result for RNA + Legacy mode:**
- Purines: `[" C1'", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "]` (9 atoms, C4 excluded)
- Pyrimidines: `[" C1'", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "]` (6 atoms, C4 excluded)

## Command-Line Tool Configuration

The `generate_modern_json` tool supports the `--legacy` flag to enable legacy compatibility mode:

```bash
# Normal mode (includes C4, correct behavior)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S.json

# Legacy mode (excludes C4, matches legacy behavior)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S.json --legacy
```

**Usage:**
```
Usage: generate_modern_json <input_pdb_file> <output_json_file> [--legacy]
Options:
  --legacy    Exclude C4 atom from matching (matches legacy behavior)
```

The tool automatically detects RNA structures by checking for O2' atoms, so the RNA flag is set automatically.

## Recommendation

**For new analyses:** Use default settings (RNA auto-detection enabled, legacy mode disabled). This provides:
- Complete atom matching (includes C4 when present)
- Better fitting accuracy (more atoms for least-squares fitting)
- RNA support (includes C1' for RNA structures)

**For legacy compatibility:** Enable legacy mode only when:
- Comparing results with legacy code output
- Reproducing legacy calculations
- Debugging legacy behavior

## Impact on Results

### RMS Fit Differences

Excluding C4 atom affects least-squares fitting results:

**Example (G T:5 from 1H4S.pdb):**
- **Normal mode:** RMS fit = 0.012142 (9 atoms including C4)
- **Legacy mode:** RMS fit = 0.008534 (8 atoms excluding C4)

The difference is expected since fewer atoms are used for fitting in legacy mode.

### Matched Atoms Differences

**Normal mode (default):**
```json
{
  "matched_atoms": [" C1'", " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "],
  "num_matched_atoms": 10
}
```

**Legacy mode:**
```json
{
  "matched_atoms": [" C1'", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "],
  "num_matched_atoms": 9
}
```

## Implementation Details

### RingAtomMatcher

The `RingAtomMatcher::match()` and `RingAtomMatcher::get_ring_atom_names()` functions accept an `exclude_c4` parameter that controls C4 atom exclusion:

```cpp
static MatchedAtoms match(const core::Residue& residue,
                          const core::Structure& standard_template,
                          bool is_rna = false,
                          bool exclude_c4 = false);

static std::vector<std::string> get_ring_atom_names(core::ResidueType residue_type,
                                                     bool is_rna = false,
                                                     bool exclude_c4 = false);
```

### BaseFrameCalculator

The `BaseFrameCalculator` internally passes the legacy mode flag to `RingAtomMatcher`:

```cpp
FrameCalculationResult BaseFrameCalculator::calculate_frame_impl(const core::Residue& residue) const {
    // ...
    MatchedAtoms matched = RingAtomMatcher::match(residue, standard_template, is_rna_, legacy_mode_);
    // ...
}
```

## See Also

- [LEGACY_VS_MODERN_MATCHED_ATOMS.md](LEGACY_VS_MODERN_MATCHED_ATOMS.md) - Detailed explanation of differences between legacy and modern matched atoms recording

