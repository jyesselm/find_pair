# Integration Test Checklist by Stage

This checklist tracks integration testing requirements at each stage where comparison with legacy JSON is possible.

## Stage 2: Core Objects ✅

**When to Run**: After implementing PDB parser and Structure classes

**Test File**: `tests/integration/test_core_objects_integration.cpp`

**Tests to Implement**:
- [ ] `ParsePDBMatchesLegacy` - Compare parsed Structure with legacy `pdb_atoms` records
- [ ] `JSONRoundTrip` - Test Structure → JSON → Structure round-trip
- [ ] `AtomDataMatches` - Compare atom names, coordinates, residue info
- [ ] `ResidueDataMatches` - Compare residue sequences, types

**Comparison Points**:
- Atom counts
- Residue counts
- Atom names and coordinates
- Residue names and sequence numbers
- Chain IDs

**Run Command**: `./tests/integration/test_core_objects_integration`

---

## Stage 3: I/O Layer ✅

**When to Run**: After implementing JSON reader/writer

**Test File**: `tests/integration/test_io_integration.cpp`

**Tests to Implement**:
- [ ] `PDBToJSONMatchesLegacy` - Compare our JSON output with legacy JSON
- [ ] `JSONToStructureMatches` - Load Structure from legacy JSON
- [ ] `RoundTripPDBJSON` - PDB → JSON → Structure → JSON
- [ ] `AllRecordTypesPresent` - Verify all calculation types are written

**Comparison Points**:
- `pdb_atoms` records
- `residue_indices` records
- JSON structure and format

**Run Command**: `./tests/integration/test_io_integration`

---

## Stage 4: Frame Calculation ✅

**When to Run**: After implementing BaseFrameCalculator

**Test File**: `tests/integration/test_frame_calculation_integration.cpp`

**Tests to Implement**:
- [ ] `FramesMatchLegacy` - Compare calculated frames with legacy `ref_frame` records
- [ ] `FrameCalculationDetailsMatch` - Compare `base_frame_calc` records
- [ ] `LSFittingMatches` - Compare `ls_fitting` records (rotation, translation, RMS)
- [ ] `FrameCalcMatches` - Compare `frame_calc` records (matched coordinates)

**Comparison Points**:
- Rotation matrices (`orien`) - tolerance 0.001
- Origins (`org`) - tolerance 0.001
- RMS fit values - tolerance 0.001
- Matched atom lists
- Template file paths

**Run Command**: `./tests/integration/test_frame_calculation_integration`

---

## Stage 5: Base Pair Finding ✅

**When to Run**: After implementing BasePairFinder

**Test File**: `tests/integration/test_base_pair_integration.cpp`

**Tests to Implement**:
- [ ] `BasePairsMatchLegacy` - Compare found base pairs with legacy `base_pair` records
- [ ] `PairValidationMatches` - Compare `pair_validation` records
- [ ] `HBondListMatches` - Compare `hbond_list` records
- [ ] `DistanceChecksMatch` - Compare `distance_checks` records
- [ ] `RingAtomsMatch` - Compare `ring_atoms` records

**Comparison Points**:
- Base pair indices (`base_i`, `base_j`)
- Base pair types
- Reference frames for each pair
- Hydrogen bond counts and details
- Validation check results
- Distance/angle measurements

**Run Command**: `./tests/integration/test_base_pair_integration`

---

## Stage 6: Parameter Calculation ✅

**When to Run**: After implementing ParameterCalculator

**Test File**: `tests/integration/test_parameter_integration.cpp`

**Tests to Implement**:
- [ ] `StepParametersMatchLegacy` - Compare `bpstep_params` records
- [ ] `HelicalParametersMatchLegacy` - Compare `helical_params` records
- [ ] `MidstepFramesMatch` - Compare midstep frames
- [ ] `AllParametersMatch` - Verify all 6 parameters for all steps

**Comparison Points**:
- **CRITICAL**: All 6 step parameters (Shift, Slide, Rise, Tilt, Roll, Twist) - tolerance 0.001
- All 6 helical parameters - tolerance 0.001
- Midstep frames (`mst_orien`, `mst_org`) - tolerance 0.001
- Helical midstep frames (`mst_orienH`, `mst_orgH`) - tolerance 0.001

**Run Command**: `./tests/integration/test_parameter_integration`

**Note**: This is the final output - must match exactly!

---

## Stage 7: Protocols ✅

**When to Run**: After implementing protocols

**Test File**: `tests/integration/test_protocol_integration.cpp`

**Tests to Implement**:
- [ ] `FindPairProtocolMatchesLegacy` - Complete find_pair workflow
- [ ] `AnalyzeProtocolMatchesLegacy` - Complete analyze workflow
- [ ] `AllCalculationTypesMatch` - Compare all record types
- [ ] `EndToEndWorkflow` - Full pipeline test

**Comparison Points**:
- All calculation types in legacy JSON
- Complete JSON structure
- All numerical values
- Record ordering (if relevant)

**Run Command**: `./tests/integration/test_protocol_integration`

---

## Running All Integration Tests

```bash
# Build all integration tests
cd build
cmake ..
make test_integration

# Run all integration tests
ctest -R integration

# Run with verbose output
ctest -R integration -VV

# Run specific test
./tests/integration/test_frame_calculation_integration
```

## Current Test Data

**PDB/JSON Pairs Available**:
- ✅ `data/pdb/100D.pdb` → `data/json_legacy/100D.json`
- ✅ `data/pdb/157D.pdb` → `data/json_legacy/157D.json`

**To Add More Test Data**:
1. Generate JSON for a PDB file using original code
2. Place in `data/json_legacy/`
3. Integration tests will automatically discover and test it

## Success Criteria

For each stage:
- [ ] Integration test file created
- [ ] Tests run on all discovered PDB/JSON pairs
- [ ] All comparisons pass within tolerance
- [ ] No test failures
- [ ] Test execution time acceptable

---

*This checklist ensures comprehensive testing against legacy JSON at every stage.*

