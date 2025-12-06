# Legacy / Modern 100% Match Plan

## 1. Close Stage-by-Stage Validation Gaps

- **stabilize-json-generation**: Fix the duplicate-record bug hitting `pdb_atoms` (seen on 7EH2) by auditing `tools/generate_modern_json.cpp` + `JsonWriter` so each record type is emitted once, remove `--fix-indices`/legacy reads per `docs/SEPARATION_OF_CONCERNS.md`, and rerun 7EH2 to confirm counts.
- **stage-harness**: Implement the staged workflow described in `docs/CLEAN_SLATE_VALIDATION_PLAN.md` (`scripts/batch_validation_workflow.py`, stage result JSONs) so we can generate/validate in batches, delete passing JSON, and log failures (focus on stages 2–6 where modern JSON already exists for ~22 PDBs).
- **stage3-6-validation**: For each remaining record type (`distance_checks`, `hbond_list`, `pair_validation`, `find_bestpair_selection`, `base_pair`), run the batch validator across available PDBs, investigate mismatches (quality score ties, hbonds, etc.), and only proceed when 100% match is recorded in `data/validation_results/`.
- **stage7-8-generation**: Extend `FindPairProtocol`/`AnalyzeProtocol` to emit `bpstep_params` and `helical_params`, generate modern JSON for the 22-PDB test set, and add pytest/CLI comparisons (`compare_json.py steps|helical`) so phases 7–8 join the staged pipeline.
- **full-scale-runs**: After all stages pass on the 22-PDB slice, scale the batch runner to `test_set_100` then the 3,602 fast PDBs, addressing timeouts via either higher limits or hotspots identified in `docs/NEXT_STEPS.md` (e.g., ribosomal structures) and updating `docs/VALIDATION_PROGRESS.md` after each milestone.

## 2. Streamline Files, Data, and Tests

- **on-demand-json**: Enforce the clean-slate policy (only generate the current stage, delete passing JSON) by wiring `generate_modern_json --stage` into the batch runner, archiving existing `data/json*/` per `docs/START_HERE.md`, and documenting space savings so the repo no longer carries stale gigabytes of outputs.
- **test-suite-thinning**: Collapse redundant Python scripts (`scripts/test_ls_fitting*.py`, manual batch runners) into pytest suites under `tests_python/` using shared fixtures from `conftest.py`, while C++ `tests/` focus on unit/integration coverage; document the unified workflow in `docs/TESTING_GUIDE.md`.
- **comparison-tooling**: Execute the `docs/REFACTOR_JSON_COMPARISON.md` plan—normalize file discovery/regeneration inside `x3dna_json_compare`, expose a single `compare`/`batch` CLI, and delete bespoke JSON diff helpers once parity is confirmed.
- **inventory-tracking**: Keep `data/validation_results/stage*_results.json` as the sole source of truth plus lightweight CSV/markdown summaries (retire ad-hoc spreadsheets/scripts) so progress reporting is automatic.

## 3. Improve Modern Code Reuse & Extensibility

- **protocol-separation**: Finish the `FrameJsonRecorder`/recorder pattern from `docs/REFACTOR_FRAME_CALCULATORS.md` by removing recording logic from `BasePairFinder` and `JsonWriter`, adding `PairJsonRecorder`/`ParameterJsonRecorder`, and letting protocols orchestrate algorithms + recorders.
- **lean-generate-modern-json**: Apply `docs/REFACTOR_GENERATE_MODERN_JSON.md`—strip legacy-dependent flags, fold RNA detection & stage setup into helpers, and shrink the tool to “parse PDB → run stage protocol → exit”, making it safe for CI/regeneration.
- **comparison-library**: Restructure `x3dna_json_compare/` per `REFACTOR_JSON_COMPARISON.md` so comparisons, regeneration, and batching share typed models and a single executor (helps future automation and external reuse).
- **future-proofing**: As recorders/protocols stabilize, codify RAII + short-function rules in headers (`include/x3dna/...`) and add extension points (e.g., dependency injection for template paths, hooks for alternative JSON sinks) so new consumers (CLI, services) can reuse the same pipelines without forked code.
