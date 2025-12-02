# Implementation Plan: Index Matching and Unified Comparison

**Date**: December 2, 2025  
**Priority**: CRITICAL - Block all other work until complete

---

## Phase 1: ResidueTracker Implementation (C++)

### 1.1 Create Header File

**File**: `include/x3dna/residue_tracker.hpp`

```cpp
#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <optional>

namespace x3dna {

struct ResidueRecord {
    // Tracking indices
    int read_index;           // Order read from PDB (0-based)
    int legacy_index;         // Index from legacy JSON (1-based, -1 if not set)
    int modern_index;         // Final index after filtering (0-based, -1 if filtered)
    bool filtered;            // Was this residue filtered out?
    std::string filter_reason;  // Why was it filtered (empty if not filtered)
    
    // PDB properties for matching
    std::string chain_id;
    int residue_seq;
    std::string insertion;
    std::string residue_name;
    
    ResidueRecord(int read_idx, const std::string& chain, int seq, 
                  const std::string& ins, const std::string& name)
        : read_index(read_idx), legacy_index(-1), modern_index(-1),
          filtered(false), filter_reason(""),
          chain_id(chain), residue_seq(seq), insertion(ins), residue_name(name) {}
};

class ResidueTracker {
public:
    ResidueTracker() = default;
    
    // Add residue as it's read from PDB (in order)
    void add_residue(const std::string& chain_id, int residue_seq,
                    const std::string& insertion, const std::string& residue_name);
    
    // Mark residue as filtered (won't appear in final list)
    void mark_filtered(int read_index, const std::string& reason);
    
    // Assign final modern index to a residue
    void assign_modern_index(int read_index, int modern_index);
    
    // Load legacy indices from JSON file
    bool load_legacy_indices(const std::string& legacy_json_path);
    
    // Validate that our indices match legacy exactly
    struct ValidationResult {
        bool success;
        int num_residues_read;
        int num_legacy;
        int num_modern;
        int num_filtered;
        int num_matched;
        int num_unmatched;
        std::vector<std::string> errors;
    };
    ValidationResult validate() const;
    
    // Export mapping to JSON for debugging
    void export_mapping(const std::string& output_path) const;
    
    // Get legacy index for a modern index
    std::optional<int> get_legacy_index(int modern_index) const;
    
    // Get modern index for a legacy index (1-based)
    std::optional<int> get_modern_index(int legacy_index) const;
    
    // Get all residues (for debugging)
    const std::vector<ResidueRecord>& get_residues() const { return residues_; }
    
private:
    std::vector<ResidueRecord> residues_;
    
    // Helper to find residue by PDB properties
    int find_by_pdb_props(const std::string& chain, int seq, 
                         const std::string& ins) const;
};

} // namespace x3dna
```

### 1.2 Implement ResidueTracker

**File**: `src/x3dna/residue_tracker.cpp`

```cpp
#include "x3dna/residue_tracker.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <sstream>
#include <iostream>

using json = nlohmann::json;

namespace x3dna {

void ResidueTracker::add_residue(const std::string& chain_id, int residue_seq,
                                 const std::string& insertion, 
                                 const std::string& residue_name) {
    int read_idx = static_cast<int>(residues_.size());
    residues_.emplace_back(read_idx, chain_id, residue_seq, insertion, residue_name);
}

void ResidueTracker::mark_filtered(int read_index, const std::string& reason) {
    if (read_index >= 0 && read_index < static_cast<int>(residues_.size())) {
        residues_[read_index].filtered = true;
        residues_[read_index].filter_reason = reason;
    }
}

void ResidueTracker::assign_modern_index(int read_index, int modern_index) {
    if (read_index >= 0 && read_index < static_cast<int>(residues_.size())) {
        residues_[read_index].modern_index = modern_index;
    }
}

bool ResidueTracker::load_legacy_indices(const std::string& legacy_json_path) {
    // Try to open file
    std::ifstream file(legacy_json_path);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open legacy JSON: " << legacy_json_path << std::endl;
        return false;
    }
    
    // Parse JSON
    json j;
    try {
        file >> j;
    } catch (const json::exception& e) {
        std::cerr << "Error parsing legacy JSON: " << e.what() << std::endl;
        return false;
    }
    
    // Load base_frame_calc records
    if (!j.contains("base_frame_calc")) {
        std::cerr << "No base_frame_calc in legacy JSON" << std::endl;
        return false;
    }
    
    for (const auto& record : j["base_frame_calc"]) {
        // Match by PDB properties
        std::string chain = record["chain_id"];
        int seq = record["residue_seq"];
        std::string ins = record.value("insertion", "");
        int legacy_idx = record["residue_idx"];  // 1-based
        
        // Find matching residue
        int read_idx = find_by_pdb_props(chain, seq, ins);
        if (read_idx >= 0) {
            residues_[read_idx].legacy_index = legacy_idx;
        }
    }
    
    return true;
}

int ResidueTracker::find_by_pdb_props(const std::string& chain, int seq,
                                      const std::string& ins) const {
    for (size_t i = 0; i < residues_.size(); ++i) {
        const auto& r = residues_[i];
        if (r.chain_id == chain && r.residue_seq == seq && r.insertion == ins) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

ResidueTracker::ValidationResult ResidueTracker::validate() const {
    ValidationResult result;
    result.success = true;
    result.num_residues_read = static_cast<int>(residues_.size());
    result.num_filtered = 0;
    result.num_legacy = 0;
    result.num_modern = 0;
    result.num_matched = 0;
    result.num_unmatched = 0;
    
    // Count residues by category
    for (const auto& r : residues_) {
        if (r.filtered) {
            result.num_filtered++;
        }
        if (r.legacy_index >= 0) {
            result.num_legacy++;
        }
        if (r.modern_index >= 0) {
            result.num_modern++;
        }
        if (r.legacy_index >= 0 && r.modern_index >= 0) {
            result.num_matched++;
        }
    }
    
    // Check: num_modern should equal num_legacy for perfect match
    if (result.num_modern != result.num_legacy) {
        result.success = false;
        std::ostringstream oss;
        oss << "Count mismatch: modern=" << result.num_modern 
            << " legacy=" << result.num_legacy;
        result.errors.push_back(oss.str());
    }
    
    // Check: all non-filtered residues should have both indices
    for (const auto& r : residues_) {
        if (!r.filtered) {
            if (r.modern_index < 0) {
                result.success = false;
                std::ostringstream oss;
                oss << "Non-filtered residue " << r.chain_id << r.residue_seq
                    << r.insertion << " has no modern index";
                result.errors.push_back(oss.str());
            }
            if (r.legacy_index < 0) {
                result.success = false;
                std::ostringstream oss;
                oss << "Non-filtered residue " << r.chain_id << r.residue_seq
                    << r.insertion << " has no legacy index";
                result.errors.push_back(oss.str());
            }
        }
    }
    
    result.num_unmatched = result.num_modern - result.num_matched;
    
    return result;
}

void ResidueTracker::export_mapping(const std::string& output_path) const {
    json j = json::array();
    
    for (const auto& r : residues_) {
        json record;
        record["read_index"] = r.read_index;
        record["legacy_index"] = r.legacy_index;
        record["modern_index"] = r.modern_index;
        record["filtered"] = r.filtered;
        record["filter_reason"] = r.filter_reason;
        record["chain_id"] = r.chain_id;
        record["residue_seq"] = r.residue_seq;
        record["insertion"] = r.insertion;
        record["residue_name"] = r.residue_name;
        j.push_back(record);
    }
    
    std::ofstream file(output_path);
    file << j.dump(2);
}

std::optional<int> ResidueTracker::get_legacy_index(int modern_index) const {
    for (const auto& r : residues_) {
        if (r.modern_index == modern_index) {
            return r.legacy_index;
        }
    }
    return std::nullopt;
}

std::optional<int> ResidueTracker::get_modern_index(int legacy_index) const {
    for (const auto& r : residues_) {
        if (r.legacy_index == legacy_index) {
            return r.modern_index;
        }
    }
    return std::nullopt;
}

} // namespace x3dna
```

### 1.3 Integration Points

**Where to integrate ResidueTracker:**

1. **PDB Parsing** (`src/x3dna/pdb_parser.cpp`):
   - Create tracker at start of parse
   - Add residue for each one read
   - Mark filtered residues
   - Assign modern indices

2. **JSON Generation** (`tools/generate_modern_json.cpp`):
   - Load legacy indices
   - Validate match
   - Export mapping
   - Abort if validation fails

3. **Testing** (`tests/test_residue_tracker.cpp`):
   - Test with known PDBs
   - Validate all test_set_10

---

## Phase 2: Unified Comparison Script (Python)

### 2.1 Create unified_compare.py

**File**: `scripts/unified_compare.py`

```python
#!/usr/bin/env python3
"""
Unified comparison framework for X3DNA legacy vs modern code.

This is the SINGLE entry point for all comparisons.
Do not create new comparison scripts - extend this one.
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set
from dataclasses import dataclass, asdict
from datetime import datetime
import subprocess

# Add x3dna_json_compare to path
sys.path.insert(0, str(Path(__file__).parent.parent / "x3dna_json_compare"))

from x3dna_json_compare.json_comparison import JSONComparer
from x3dna_json_compare.json_file_finder import find_json_files

@dataclass
class ComparisonStatus:
    """Track status of comparison stages for a PDB"""
    pdb_id: str
    index_match: str = "PENDING"  # PASS, FAIL, PENDING, SKIP
    atoms: str = "PENDING"
    frames: str = "PENDING"
    hbonds: str = "PENDING"
    pairs: str = "PENDING"
    step_params: str = "PENDING"
    last_updated: str = ""
    notes: str = ""
    
    def to_dict(self):
        return asdict(self)
    
    @staticmethod
    def from_dict(d: Dict) -> 'ComparisonStatus':
        return ComparisonStatus(**d)


class StatusTracker:
    """Manages comparison_status.csv"""
    
    def __init__(self, csv_path: Path):
        self.csv_path = csv_path
        self.statuses: Dict[str, ComparisonStatus] = {}
        self.load()
    
    def load(self):
        """Load status from CSV"""
        if not self.csv_path.exists():
            return
        
        with open(self.csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pdb_id = row['pdb_id']
                self.statuses[pdb_id] = ComparisonStatus.from_dict(row)
    
    def save(self):
        """Save status to CSV"""
        self.csv_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(self.csv_path, 'w', newline='') as f:
            fieldnames = list(ComparisonStatus.__dataclass_fields__.keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for status in sorted(self.statuses.values(), key=lambda s: s.pdb_id):
                writer.writerow(status.to_dict())
    
    def get(self, pdb_id: str) -> ComparisonStatus:
        """Get status for PDB, create if doesn't exist"""
        if pdb_id not in self.statuses:
            self.statuses[pdb_id] = ComparisonStatus(pdb_id=pdb_id)
        return self.statuses[pdb_id]
    
    def update(self, pdb_id: str, stage: str, status: str, notes: str = ""):
        """Update status for a stage"""
        s = self.get(pdb_id)
        setattr(s, stage, status)
        s.last_updated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if notes:
            s.notes = notes
        self.save()
    
    def get_summary(self) -> Dict[str, int]:
        """Get summary statistics"""
        stages = ['index_match', 'atoms', 'frames', 'hbonds', 'pairs', 'step_params']
        summary = {stage: {'PASS': 0, 'FAIL': 0, 'PENDING': 0, 'SKIP': 0} 
                  for stage in stages}
        
        for status in self.statuses.values():
            for stage in stages:
                val = getattr(status, stage)
                summary[stage][val] = summary[stage].get(val, 0) + 1
        
        return summary


class UnifiedComparer:
    """Main comparison class"""
    
    def __init__(self, project_root: Path, cache: bool = False):
        self.project_root = project_root
        self.data_dir = project_root / "data"
        self.status_csv = self.data_dir / "comparison_status.csv"
        self.tracker = StatusTracker(self.status_csv)
        self.cache = cache
        
        # Paths
        self.legacy_json_dir = self.data_dir / "json_legacy"
        self.modern_json_dir = self.data_dir / "json"
        self.index_mapping_dir = self.data_dir / "index_mapping"
        self.build_dir = project_root / "build"
        
    def validate_indices(self, pdb_ids: List[str]) -> Dict[str, bool]:
        """Validate residue index matching for PDBs"""
        results = {}
        
        for pdb_id in pdb_ids:
            mapping_file = self.index_mapping_dir / f"{pdb_id}.json"
            
            if not mapping_file.exists():
                print(f"{pdb_id}: No index mapping file - FAIL")
                self.tracker.update(pdb_id, 'index_match', 'FAIL', 
                                  'No index mapping file')
                results[pdb_id] = False
                continue
            
            # Load and validate mapping
            with open(mapping_file) as f:
                mapping = json.load(f)
            
            # Check for issues
            num_read = len(mapping)
            num_legacy = sum(1 for r in mapping if r['legacy_index'] >= 0)
            num_modern = sum(1 for r in mapping if r['modern_index'] >= 0)
            num_filtered = sum(1 for r in mapping if r['filtered'])
            
            if num_modern != num_legacy:
                msg = f"Count mismatch: modern={num_modern} legacy={num_legacy}"
                print(f"{pdb_id}: {msg} - FAIL")
                self.tracker.update(pdb_id, 'index_match', 'FAIL', msg)
                results[pdb_id] = False
                continue
            
            # Check all non-filtered have both indices
            issues = []
            for r in mapping:
                if not r['filtered']:
                    if r['modern_index'] < 0:
                        issues.append(f"No modern index for {r['chain_id']}{r['residue_seq']}")
                    if r['legacy_index'] < 0:
                        issues.append(f"No legacy index for {r['chain_id']}{r['residue_seq']}")
            
            if issues:
                msg = "; ".join(issues[:3])  # First 3 issues
                print(f"{pdb_id}: {msg} - FAIL")
                self.tracker.update(pdb_id, 'index_match', 'FAIL', msg)
                results[pdb_id] = False
                continue
            
            # PASS
            msg = f"{num_modern} residues, {num_filtered} filtered"
            print(f"{pdb_id}: {msg} - PASS")
            self.tracker.update(pdb_id, 'index_match', 'PASS', msg)
            results[pdb_id] = True
        
        return results
    
    def compare_json_stage(self, pdb_id: str, stage: str, verbose: bool = False):
        """Compare a specific JSON stage"""
        # Check index validation first
        status = self.tracker.get(pdb_id)
        if status.index_match != 'PASS':
            print(f"SKIP {pdb_id} {stage}: Index validation not passed")
            self.tracker.update(pdb_id, stage, 'SKIP', 
                              'Index validation failed')
            return
        
        # TODO: Implement actual comparison using x3dna_json_compare
        # This is a placeholder
        pass
    
    def batch_compare(self, pdb_ids: List[str], stage: Optional[str] = None,
                     verbose: bool = False, resume: bool = False):
        """Batch comparison with status tracking"""
        stages = ['index_match', 'atoms', 'frames', 'hbonds', 'pairs', 'step_params']
        
        if stage:
            # Single stage
            if stage == 'index_match':
                self.validate_indices(pdb_ids)
            else:
                for pdb_id in pdb_ids:
                    if resume:
                        status = self.tracker.get(pdb_id)
                        if getattr(status, stage) == 'PASS':
                            print(f"SKIP {pdb_id} {stage}: Already passed")
                            continue
                    self.compare_json_stage(pdb_id, stage, verbose)
        else:
            # All stages
            for s in stages:
                print(f"\n=== Stage: {s} ===")
                if s == 'index_match':
                    self.validate_indices(pdb_ids)
                else:
                    for pdb_id in pdb_ids:
                        if resume:
                            status = self.tracker.get(pdb_id)
                            if getattr(status, s) == 'PASS':
                                continue
                        self.compare_json_stage(pdb_id, s, verbose)
    
    def print_status(self):
        """Print status summary"""
        summary = self.tracker.get_summary()
        
        print("\n=== Comparison Status Summary ===\n")
        print(f"{'Stage':<15} {'PASS':<8} {'FAIL':<8} {'PENDING':<8} {'SKIP':<8}")
        print("-" * 55)
        
        for stage, counts in summary.items():
            print(f"{stage:<15} {counts['PASS']:<8} {counts['FAIL']:<8} "
                  f"{counts['PENDING']:<8} {counts['SKIP']:<8}")
        
        print(f"\nTotal PDBs: {len(self.tracker.statuses)}")
        print(f"Status file: {self.status_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Unified comparison framework for X3DNA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate indices
  %(prog)s validate-indices 1H4S 2BNA
  %(prog)s validate-indices --test-set 10
  
  # Compare single stage
  %(prog)s batch --stage atoms --test-set 10
  
  # Compare all stages
  %(prog)s batch --test-set 100
  
  # Resume from last run
  %(prog)s batch --test-set 100 --resume
  
  # View status
  %(prog)s status
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # validate-indices command
    validate_parser = subparsers.add_parser('validate-indices',
                                           help='Validate residue index matching')
    validate_parser.add_argument('pdb_ids', nargs='*', help='PDB IDs to validate')
    validate_parser.add_argument('--test-set', type=int, choices=[10, 50, 100, 500, 1000],
                                help='Use test set')
    
    # batch command
    batch_parser = subparsers.add_parser('batch', help='Batch comparison')
    batch_parser.add_argument('--stage', choices=['index_match', 'atoms', 'frames', 
                                                  'hbonds', 'pairs', 'step_params'],
                             help='Compare specific stage only')
    batch_parser.add_argument('--test-set', type=int, choices=[10, 50, 100, 500, 1000],
                             help='Use test set')
    batch_parser.add_argument('--resume', action='store_true',
                             help='Skip already passed PDBs')
    batch_parser.add_argument('--verbose', action='store_true',
                             help='Verbose output')
    batch_parser.add_argument('--cache', action='store_true',
                             help='Enable caching (disabled by default)')
    
    # status command
    status_parser = subparsers.add_parser('status', help='Show status summary')
    
    args = parser.parse_args()
    
    # Get project root
    project_root = Path(__file__).parent.parent
    
    # Create comparer
    comparer = UnifiedComparer(project_root, cache=getattr(args, 'cache', False))
    
    # Execute command
    if args.command == 'validate-indices':
        if args.test_set:
            # TODO: Load test set
            pdb_ids = []  # placeholder
        else:
            pdb_ids = args.pdb_ids
        comparer.validate_indices(pdb_ids)
    
    elif args.command == 'batch':
        if args.test_set:
            # TODO: Load test set
            pdb_ids = []  # placeholder
        else:
            print("Error: --test-set required for batch command")
            sys.exit(1)
        comparer.batch_compare(pdb_ids, args.stage, args.verbose, args.resume)
    
    elif args.command == 'status':
        comparer.print_status()


if __name__ == '__main__':
    main()
```

---

## Phase 3: Action Checklist

### Immediate Actions (Do First)

- [ ] 1. Review and approve CLEANUP_AND_RESTRUCTURE_PLAN.md
- [ ] 2. Create feature branch: `git checkout -b fix-index-matching`
- [ ] 3. Implement ResidueTracker class (C++)
- [ ] 4. Add ResidueTracker to CMakeLists.txt
- [ ] 5. Integrate ResidueTracker into PDB parsing
- [ ] 6. Test on single PDB (1H4S)
- [ ] 7. Generate index mapping for 1H4S
- [ ] 8. Validate 1H4S mapping against legacy JSON
- [ ] 9. Fix any issues found
- [ ] 10. Test on test_set_10

### Week 1: Index Validation

- [ ] 11. Generate index mappings for all test_set_10
- [ ] 12. Validate all test_set_10 mappings
- [ ] 13. Document any failures
- [ ] 14. Fix root causes of failures
- [ ] 15. Achieve 100% pass on test_set_10
- [ ] 16. Generate mappings for test_set_100
- [ ] 17. Validate test_set_100
- [ ] 18. Create index_validation_status.csv
- [ ] 19. Commit index validation changes
- [ ] 20. Merge to main

### Week 2: Unified Framework

- [ ] 21. Create scripts/unified_compare.py
- [ ] 22. Implement StatusTracker class
- [ ] 23. Implement validate_indices command
- [ ] 24. Test validate_indices on test_set_10
- [ ] 25. Implement batch comparison skeleton
- [ ] 26. Implement status command
- [ ] 27. Test on test_set_10
- [ ] 28. Add comparison stages (atoms, frames, etc.)
- [ ] 29. Test full workflow
- [ ] 30. Commit unified framework

### Week 3: Tool Consolidation

- [ ] 31. Archive tools/ to tools/archive/
- [ ] 32. Update CMakeLists.txt to only build essential tools
- [ ] 33. Test build
- [ ] 34. Archive scripts/ to scripts/archive/
- [ ] 35. Keep only: unified_compare.py, rebuild_json.py, download_pdbs.py
- [ ] 36. Test scripts still work
- [ ] 37. Update Makefile
- [ ] 38. Commit consolidation

### Week 4: Documentation Cleanup

- [ ] 39. Create docs/archive/ if needed
- [ ] 40. Move old docs to archive
- [ ] 41. Keep only 6-7 core docs
- [ ] 42. Update README.md
- [ ] 43. Update CODE_FLOW.md
- [ ] 44. Simplify TESTING_GUIDE.md
- [ ] 45. Create docs/README.md
- [ ] 46. Test that docs make sense
- [ ] 47. Remove redundant files
- [ ] 48. Commit documentation changes

### Week 5: Full Validation

- [ ] 49. Run unified_compare.py on test_set_100
- [ ] 50. Investigate all failures
- [ ] 51. Fix root causes
- [ ] 52. Re-run until 100% on atoms stage
- [ ] 53. Re-run until 100% on frames stage
- [ ] 54. Re-run until 100% on hbonds stage
- [ ] 55. Analyze pairs stage results
- [ ] 56. Document path to 100% match
- [ ] 57. Create final report
- [ ] 58. Merge to main
- [ ] 59. Tag release: v1.0-index-matching
- [ ] 60. Celebrate! 🎉

---

## Success Criteria

### Phase 1: Index Matching (Week 1)
✅ 100% index match validation for test_set_10
✅ Index mapping JSON files for all test PDBs
✅ Clear documentation of any mismatches
✅ Automated validation in build process

### Phase 2: Unified Framework (Week 2)
✅ Single entry point for all comparisons
✅ CSV status tracking working
✅ Can validate indices for any PDB
✅ Can resume from failures

### Phase 3: Consolidation (Week 3)
✅ 36 tools → 1 tool + archive
✅ 60 scripts → 3 scripts + archive
✅ Build still works
✅ All tests pass

### Phase 4: Documentation (Week 4)
✅ 50 docs → 7 core docs
✅ README is clear and concise
✅ Easy for new developer to understand

### Phase 5: Validation (Week 5)
✅ 100% pass on atoms for test_set_100
✅ 100% pass on frames for test_set_100
✅ 100% pass on hbonds for test_set_100
✅ Clear path to fixing remaining pairs

---

## Notes

- **BLOCK ALL OTHER WORK** until index validation is 100%
- **NO NEW TOOLS** - extend unified_compare.py instead
- **CSV IS TRUTH** - Always check comparison_status.csv first
- **NO CACHING BY DEFAULT** - Only enable with --cache flag
- **VALIDATE EARLY** - Run validation after every change

---

## Questions & Answers

**Q: What if a PDB doesn't validate?**
A: Document in index_validation_status.csv, investigate root cause, fix before proceeding.

**Q: What if we need a new comparison type?**
A: Add it as a stage to unified_compare.py, don't create new script.

**Q: What about the cluster tools?**
A: Keep them but integrate into unified framework as batch mode.

**Q: Can we delete archived files?**
A: Keep them in archive/ for now, delete after confirming not needed (6 months).

**Q: What about x3dna_json_compare module?**
A: Keep it - it's the comparison library used by unified_compare.py.


