#!/bin/bash
# Script to run both original and our code, then compare debug outputs

PDB_FILE="${1:-../legacy_archive/fixtures/pdb/1EHZ.pdb}"

if [ ! -f "$PDB_FILE" ]; then
    echo "Error: PDB file not found: $PDB_FILE"
    exit 1
fi

echo "=========================================="
echo "Comparing Original vs Our Code"
echo "=========================================="
echo "PDB file: $PDB_FILE"
echo ""

# Run original
echo "1. Running original x3dna..."
cd org
./build/bin/find_pair_original "$PDB_FILE" > /tmp/orig_output.inp 2> /tmp/orig_debug.txt
echo "   Original debug: $(wc -l < /tmp/orig_debug.txt) lines"
cd ..

# Run ours (if built)
if [ -f "build/debug/find_pair" ]; then
    echo "2. Running our code..."
    ./build/debug/find_pair "$PDB_FILE" > /tmp/ours_output.inp 2> /tmp/ours_debug.txt
    echo "   Our debug: $(wc -l < /tmp/ours_debug.txt) lines"
else
    echo "2. Our code not built. Build it first:"
    echo "   ./scripts/build.sh debug"
    exit 1
fi

echo ""
echo "3. Extracting first pair comparison..."
echo ""

# Extract first pair from original
echo "=== ORIGINAL (first pair) ==="
grep -A 20 "\[DEBUG calculate_more_bppars\] i=1 j=" /tmp/orig_debug.txt | head -25

echo ""
echo "=== OUR CODE (first pair) ==="
grep -A 20 "\[DEBUG calculate_more_bppars\] i=1 j=" /tmp/ours_debug.txt | head -25

echo ""
echo "4. Parameter comparison:"
echo ""

# Extract parameters
orig_params=$(grep "\[DEBUG\] Result:" /tmp/orig_debug.txt | head -1 | sed 's/.*Result: //')
ours_params=$(grep "\[DEBUG\] Result:" /tmp/ours_debug.txt | head -1 | sed 's/.*Result: //')

if [ -n "$orig_params" ] && [ -n "$ours_params" ]; then
    echo "Original: $orig_params"
    echo "Our:      $ours_params"
    
    # Use Python for comparison
    python3 << EOF
orig = "$orig_params".split()
ours = "$ours_params".split()
params = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist']

print("\nDifferences:")
print("-" * 60)
for i, param in enumerate(params):
    if i < len(orig) and i < len(ours):
        try:
            orig_val = float(orig[i])
            ours_val = float(ours[i])
            diff = abs(orig_val - ours_val)
            match = "✓" if diff < 0.001 else "✗"
            print(f"{param:6s}: Diff={diff:10.6f} {match}")
        except:
            pass
EOF
fi

echo ""
echo "Debug files saved:"
echo "  Original: /tmp/orig_debug.txt"
echo "  Our code: /tmp/ours_debug.txt"

