================================================================================
REMAINING 0.5% DIFFERENCES - ANALYSIS & RECOMMENDATIONS
================================================================================

SUMMARY
--------------------------------------------------------------------------------
Total PDBs with differences: 13
Total missing pairs: 46
Total extra pairs: 52

PATTERN ANALYSIS
--------------------------------------------------------------------------------
Same residue, different partner: 52 cases
  This suggests quality score differences or tie-breaking issues

Examples:
  3CF5: Missing (3239, 3680), Extra (3238, 3680) (shared: 3680)
  3CF5: Missing (3239, 3680), Extra (3239, 3679) (shared: 3239)
  3CF5: Missing (4202, 4206), Extra (4203, 4206) (shared: 4206)
  3CF5: Missing (5130, 5131), Extra (3684, 5130) (shared: 5130)
  3CF5: Missing (5858, 5873), Extra (5856, 5873) (shared: 5873)
  3CF5: Missing (3297, 3332), Extra (3297, 3333) (shared: 3297)
  3CF5: Missing (3786, 3788), Extra (3786, 3789) (shared: 3786)
  3CF5: Missing (3236, 3238), Extra (3238, 3680) (shared: 3238)
  3CF5: Missing (3236, 3238), Extra (3236, 3237) (shared: 3236)
  3CF5: Missing (4807, 5085), Extra (5081, 5085) (shared: 5085)
  ... and 42 more

Adjacent/close residues: 19 pairs
  These are likely tie-breaking cases where quality scores are very close

TOP PDBs WITH DIFFERENCES
--------------------------------------------------------------------------------
3CF5: 9 missing, 14 extra (total: 23)
3CME: 7 missing, 8 extra (total: 15)
6G5I: 6 missing, 5 extra (total: 11)
4JV5: 5 missing, 5 extra (total: 10)
6CAP: 4 missing, 5 extra (total: 9)
5IWA: 3 missing, 4 extra (total: 7)
8T2T: 3 missing, 3 extra (total: 6)
6J6G: 3 missing, 2 extra (total: 5)
7XHT: 2 missing, 3 extra (total: 5)
3IWN: 1 missing, 1 extra (total: 2)

================================================================================
RECOMMENDATIONS
================================================================================

1. QUALITY SCORE INVESTIGATION
   - Use C++ tool to compare quality scores for specific pairs:
     build/compare_quality_scores <pdb_id> <res1> <res2>
   - Focus on pairs with same residue, different partner
   - Check adjusted quality score calculation (base + adjust_pairQuality - bp_type_id adjustment)

2. TIE-BREAKING INVESTIGATION
   - When quality scores are equal, check iteration order
   - Verify that modern code uses same iteration order as legacy
   - Check that matched residues are excluded correctly

3. H-BOND ADJUSTMENT INVESTIGATION
   - Verify adjust_pairQuality calculation matches legacy
   - Check H-bond counting (good H-bonds: distance in [2.5, 3.5])
   - Formula: adjust = (num_good_hb >= 2 ? -3.0 : -num_good_hb)

4. BP_TYPE_ID INVESTIGATION
   - Verify bp_type_id calculation matches legacy
   - Check final score adjustment: -2.0 if bp_type_id == 2

5. SPECIFIC PAIRS TO INVESTIGATE
   - Focus on top PDBs with most differences:
     3CF5: 9 missing, 14 extra
       Missing: [(3239, 3680), (4202, 4206), (5130, 5131)]
       Extra: [(3684, 5130), (3238, 3680), (3786, 3789)]
     3CME: 7 missing, 8 extra
       Missing: [(4851, 4981), (4045, 4132), (6112, 6144)]
       Extra: [(4726, 4740), (4046, 4132), (3806, 4206)]
     6G5I: 6 missing, 5 extra
       Missing: [(857, 904), (1153, 1157), (219, 223)]
       Extra: [(858, 904), (1157, 1307), (113, 223)]
     4JV5: 5 missing, 5 extra
       Missing: [(169, 198), (9, 16), (17, 887)]
       Extra: [(9, 887), (1042, 1167), (170, 198)]
     6CAP: 4 missing, 5 extra
       Missing: [(499, 512), (1123, 1125), (495, 498)]
       Extra: [(495, 512), (499, 508), (1105, 1125)]
