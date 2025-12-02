#!/bin/bash
# Test ref_frames with a single base pair to isolate the issue

set -e

PDB_FILE="${1:-data/pdb/100D.pdb}"
TEST_DIR="/tmp/test_ref_frames_single"
mkdir -p "$TEST_DIR"

echo "=== Testing ref_frames with single pair ==="
echo "PDB: $PDB_FILE"
echo "Test dir: $TEST_DIR"
echo ""

# Step 1: Generate legacy ref_frames
echo "Step 1: Generating legacy ref_frames..."
cd "$(dirname "$0")/.."
org/build/bin/find_pair_original "$PDB_FILE" > "$TEST_DIR/legacy_find_pair.log" 2>&1 || true

if [ -f "ref_frames.dat" ]; then
    mv ref_frames.dat "$TEST_DIR/legacy_ref_frames.dat"
    echo "  ✓ Legacy ref_frames.dat created"
else
    echo "  ✗ Legacy ref_frames.dat not found"
    exit 1
fi

# Step 2: Generate modern ref_frames (without legacy ordering)
echo ""
echo "Step 2: Generating modern ref_frames (no legacy ordering)..."
./build/find_pair_app "$PDB_FILE" "$TEST_DIR/modern.inp" > "$TEST_DIR/modern_find_pair.log" 2>&1

if [ -f "ref_frames_modern.dat" ]; then
    mv ref_frames_modern.dat "$TEST_DIR/modern_ref_frames_no_ordering.dat"
    echo "  ✓ Modern ref_frames_modern.dat created"
else
    echo "  ✗ Modern ref_frames_modern.dat not found"
    exit 1
fi

# Step 3: Generate modern ref_frames (with legacy ordering)
echo ""
echo "Step 3: Generating modern ref_frames (with legacy ordering)..."
# First, we need a legacy .inp file - use the modern one for now
if [ -f "$TEST_DIR/modern.inp" ]; then
    ./build/find_pair_app --legacy-inp="$TEST_DIR/modern.inp" "$PDB_FILE" "$TEST_DIR/modern2.inp" > "$TEST_DIR/modern_find_pair_ordered.log" 2>&1
    
    if [ -f "ref_frames_modern.dat" ]; then
        mv ref_frames_modern.dat "$TEST_DIR/modern_ref_frames_with_ordering.dat"
        echo "  ✓ Modern ref_frames_modern.dat (with ordering) created"
    fi
fi

# Step 4: Compare
echo ""
echo "Step 4: Comparing ref_frames..."
echo ""

echo "=== Comparison: Legacy vs Modern (no ordering) ==="
python3 scripts/debug_ref_frames.py \
    "$TEST_DIR/legacy_ref_frames.dat" \
    "$TEST_DIR/modern_ref_frames_no_ordering.dat" \
    --verbose

if [ -f "$TEST_DIR/modern_ref_frames_with_ordering.dat" ]; then
    echo ""
    echo "=== Comparison: Legacy vs Modern (with ordering) ==="
    python3 scripts/debug_ref_frames.py \
        "$TEST_DIR/legacy_ref_frames.dat" \
        "$TEST_DIR/modern_ref_frames_with_ordering.dat" \
        --verbose
fi

echo ""
echo "=== Test complete ==="
echo "Results in: $TEST_DIR"

