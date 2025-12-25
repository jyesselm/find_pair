# Plan: Template-Based Pair Classification with RMSD + H-Bond Scoring

## Goal

For each candidate base pair:
1. Apply ALL LW class templates (cWW, tSH, tWH, cSW, etc.)
2. Compute RMSD for each template alignment
3. Compute H-bond score based on expected H-bond patterns
4. Determine best classification and confidence (gap to second-best)

## Why This Matters

- **cWW pairs**: Should have lowest RMSD to cWW template AND best H-bond match
- **Other pairs (tSH, tWH, cSW)**: Need this approach since N1N9 distance alone doesn't discriminate
- **Confidence metric**: Gap between best and second-best template helps identify uncertain classifications

## Data We Have

### Templates
- `basepair-idealized/` - Idealized templates for each LW class + sequence
- `basepair-exemplars/` - Exemplar templates from real structures
- Format: `{sequence}_{LW}.pdb` (e.g., `GC_cWW.pdb`, `GA_tSH.pdb`)

### H-Bond Data
- `data/json/all_hbond_list/` - Modern H-bond output with:
  - donor_atom, acceptor_atom, distance
  - context (base_base, base_backbone, etc.)
  - classification (standard, non_standard)

### Expected H-Bond Patterns (per LW class)
```
cWW (Watson-Crick):
  GC: O6-N4, N1-N3, N2-O2 (3 H-bonds)
  AU: N6-O4, N1-N3 (2 H-bonds)
  GU: O6-N3, N1-O2 (2 H-bonds, wobble)

tSH (trans Sugar-Hoogsteen):
  GA: N2-N7, N3-N6 (2 H-bonds)

tWH (trans Watson-Hoogsteen):
  GG: O6-N7, N1-N2 (2 H-bonds)

cSW (cis Sugar-Watson):
  GC: N2-O2, O2'-N3 (2 H-bonds)
```

## Implementation Steps

### Phase 1: Template Alignment Engine

**File**: `prototypes/pair_identification/template_aligner.py`

```python
class TemplateAligner:
    """Align target pair to template and compute RMSD."""

    def __init__(self, template_dir: Path):
        self.templates = self._load_templates(template_dir)

    def align_to_template(
        self,
        target_res1: Residue,
        target_res2: Residue,
        template_path: Path
    ) -> Tuple[float, np.ndarray]:
        """
        Align target pair to template using Kabsch algorithm.

        Returns:
            rmsd: RMSD after optimal superposition
            rotation: Rotation matrix applied
        """
        # 1. Extract common atoms (C1', N1/N9, base ring atoms)
        # 2. Apply Kabsch algorithm
        # 3. Return RMSD

    def classify_pair(
        self,
        target_res1: Residue,
        target_res2: Residue,
        sequence: str,
        lw_classes: List[str]
    ) -> ClassificationResult:
        """
        Try all templates and return best match.

        Returns:
            best_lw: Best matching LW class
            best_rmsd: RMSD to best template
            second_best_lw: Second best LW class
            second_best_rmsd: RMSD to second best
            confidence: (second_best_rmsd - best_rmsd)
        """
```

### Phase 2: H-Bond Pattern Scorer

**File**: `prototypes/pair_identification/hbond_scorer.py`

```python
# Expected H-bond patterns for each LW class
HBOND_PATTERNS = {
    "cWW": {
        "GC": [("O6", "N4"), ("N1", "N3"), ("N2", "O2")],
        "AU": [("N6", "O4"), ("N1", "N3")],
        "GU": [("O6", "N3"), ("N1", "O2")],
        # ...
    },
    "tSH": {
        "GA": [("N2", "N7"), ("N3", "N6")],
        "AG": [("N7", "N2"), ("N6", "N3")],
        # ...
    },
    "tWH": {
        "GG": [("O6", "N7"), ("N1", "N2")],
        # ...
    },
    # ... other LW classes
}

class HBondScorer:
    """Score H-bond pattern match for each LW class."""

    def score_hbonds(
        self,
        observed_hbonds: List[HBond],
        lw_class: str,
        sequence: str
    ) -> float:
        """
        Score how well observed H-bonds match expected pattern.

        Returns:
            score: 0.0 to 1.0 (1.0 = perfect match)
        """
        expected = HBOND_PATTERNS[lw_class][sequence]
        matched = 0
        for donor, acceptor in expected:
            for hb in observed_hbonds:
                if hb.donor_atom == donor and hb.acceptor_atom == acceptor:
                    matched += 1
                    break
        return matched / len(expected) if expected else 0.0
```

### Phase 3: Combined Classifier

**File**: `prototypes/pair_identification/lw_classifier.py`

```python
@dataclass
class ClassificationResult:
    sequence: str
    best_lw: str
    best_rmsd: float
    best_hbond_score: float
    best_combined_score: float

    second_lw: str
    second_rmsd: float
    second_hbond_score: float
    second_combined_score: float

    confidence: float  # gap between best and second

class LWClassifier:
    """Classify base pairs using template RMSD + H-bond scoring."""

    def __init__(
        self,
        template_dir: Path,
        rmsd_weight: float = 0.5,
        hbond_weight: float = 0.5
    ):
        self.aligner = TemplateAligner(template_dir)
        self.hbond_scorer = HBondScorer()
        self.rmsd_weight = rmsd_weight
        self.hbond_weight = hbond_weight

    def classify(
        self,
        res1: Residue,
        res2: Residue,
        hbonds: List[HBond],
        lw_classes: List[str] = ["cWW", "tSH", "tWH", "cSW", "cWS"]
    ) -> ClassificationResult:
        """
        Classify pair by trying all LW templates.

        For each LW class:
        1. Compute RMSD to template
        2. Compute H-bond pattern score
        3. Combine: score = rmsd_weight * (1/rmsd) + hbond_weight * hbond_score

        Return best and second-best with confidence.
        """
        results = []

        for lw in lw_classes:
            rmsd = self.aligner.align_to_template(res1, res2, lw)
            hbond_score = self.hbond_scorer.score_hbonds(hbonds, lw, sequence)

            # Combine scores (lower RMSD is better, higher H-bond is better)
            combined = self.hbond_weight * hbond_score - self.rmsd_weight * rmsd

            results.append((lw, rmsd, hbond_score, combined))

        # Sort by combined score (descending)
        results.sort(key=lambda x: x[3], reverse=True)

        best = results[0]
        second = results[1] if len(results) > 1 else None

        return ClassificationResult(
            best_lw=best[0],
            best_rmsd=best[1],
            best_hbond_score=best[2],
            best_combined_score=best[3],
            second_lw=second[0] if second else None,
            second_rmsd=second[1] if second else None,
            # ...
            confidence=best[3] - second[3] if second else 1.0
        )
```

### Phase 4: Validation Against DSSR

**File**: `prototypes/pair_identification/validate_classification.py`

```python
def validate_classifier(
    pdb_id: str,
    classifier: LWClassifier,
    dssr_path: Path,
    pdb_path: Path,
    hbond_path: Path
) -> Dict:
    """
    Validate classifier against DSSR ground truth.

    For each DSSR pair:
    1. Run classifier
    2. Check if best_lw matches DSSR LW
    3. Compute accuracy, confusion matrix
    """
    # Load DSSR pairs with LW classifications
    # Run classifier on each
    # Compare results

    return {
        "accuracy": correct / total,
        "confusion_matrix": {...},
        "low_confidence_pairs": [...],  # Where classifier was uncertain
    }
```

## Scoring Formula

```
Combined Score = (hbond_weight × hbond_score) + (rmsd_weight × rmsd_score)

where:
  hbond_score = matched_hbonds / expected_hbonds  (0 to 1)
  rmsd_score = max(0, 1 - rmsd / rmsd_cutoff)     (0 to 1, cutoff ~2.0Å)

Default weights:
  hbond_weight = 0.6  (H-bonds are more discriminative)
  rmsd_weight = 0.4   (RMSD helps break ties)
```

## Template Selection Strategy

For each sequence (e.g., GC, GA), try templates in this order:
1. **Idealized templates** (`basepair-idealized/`) - Perfect geometry
2. **Exemplar templates** (`basepair-exemplars/`) - Real structures

If idealized not available, fall back to exemplar.

## Expected Output

For each pair:
```json
{
  "res_id1": "A-G-1",
  "res_id2": "A-C-72",
  "sequence": "GC",
  "classification": {
    "best": {
      "lw_class": "cWW",
      "rmsd": 0.45,
      "hbond_score": 1.0,
      "combined_score": 0.82
    },
    "second": {
      "lw_class": "cSW",
      "rmsd": 1.89,
      "hbond_score": 0.33,
      "combined_score": 0.12
    },
    "confidence": 0.70
  }
}
```

## Success Criteria

1. **cWW accuracy**: ≥95% match with DSSR cWW classification
2. **Other LW classes**: ≥85% accuracy on tSH, tWH, cSW
3. **Confidence calibration**: Low confidence pairs should have higher error rate
4. **Speed**: <100ms per pair classification

## Files to Create

1. `template_aligner.py` - Kabsch alignment + RMSD calculation
2. `hbond_scorer.py` - H-bond pattern matching
3. `lw_classifier.py` - Combined classifier
4. `validate_classification.py` - Validation against DSSR
5. `hbond_patterns.json` - Expected H-bond patterns for all LW classes

## Next Steps

1. [ ] Extract H-bond patterns for all LW classes from DSSR data
2. [ ] Implement template_aligner.py with Kabsch algorithm
3. [ ] Implement hbond_scorer.py with pattern matching
4. [ ] Implement lw_classifier.py combining both scores
5. [ ] Validate on cWW pairs first (should get ~95%+)
6. [ ] Extend to tSH, tWH, cSW and validate
7. [ ] Tune weights based on validation results
