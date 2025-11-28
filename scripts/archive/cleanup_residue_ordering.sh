#!/bin/bash
# Clean up residue ordering JSON files, keeping only essential files

cd "$(dirname "$0")/.."

RESIDUE_ORDERING_DIR="data/residue_ordering"

echo "🧹 Cleaning up residue ordering directory..."

# Count files before cleanup
BEFORE=$(find "$RESIDUE_ORDERING_DIR" -name "*_modern.json" -o -name "*_legacy.json" 2>/dev/null | wc -l | tr -d ' ')

# Delete individual modern and legacy JSON files
echo "  Deleting individual JSON files..."
find "$RESIDUE_ORDERING_DIR" -name "*_modern.json" -delete
find "$RESIDUE_ORDERING_DIR" -name "*_legacy.json" -delete

# Keep summary files, test files, and the log
echo "  Keeping summary files, test files, and logs..."

# Count files after cleanup
AFTER=$(find "$RESIDUE_ORDERING_DIR" -type f 2>/dev/null | wc -l | tr -d ' ')

echo ""
echo "✅ Cleanup complete!"
echo "   Deleted: $BEFORE individual JSON files"
echo "   Remaining: $AFTER files (summary, tests, logs)"
echo ""
echo "   To regenerate JSON files, run:"
echo "   python3 scripts/generate_and_compare_residue_ordering_batch.py 1000"

