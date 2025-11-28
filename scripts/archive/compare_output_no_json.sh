#!/bin/bash
# Compare legacy vs modern output without JSON generation

set -e

PDB_DIR="data/pdb"
LEGACY_BIN="org/build/bin/find_pair"
MODERN_BIN="build/find_pair_app"
LEGACY_ANALYZE="org/build/bin/analyze"
MODERN_ANALYZE="build/analyze_app"

TEMP_DIR=$(mktemp -d)
echo "Using temp directory: $TEMP_DIR"

if [ $# -eq 0 ]; then
    echo "Usage: $0 <pdb1> [pdb2] [pdb3] ..."
    echo "Example: $0 6V9Q 7EH2"
    exit 1
fi

for PDB in "$@"; do
    PDB_FILE="$PDB_DIR/${PDB}.pdb"
    
    if [ ! -f "$PDB_FILE" ]; then
        echo "ERROR: PDB file not found: $PDB_FILE"
        continue
    fi
    
    echo ""
    echo "========================================"
    echo "Comparing: $PDB"
    echo "========================================"
    
    LEGACY_INP="$TEMP_DIR/${PDB}_legacy.inp"
    MODERN_INP="$TEMP_DIR/${PDB}_modern.inp"
    
    # Step 1: Run find_pair (legacy with --no-json)
    echo ""
    echo "--- Step 1: Running find_pair ---"
    echo "Legacy:"
    if [ -f "$LEGACY_BIN" ]; then
        cd org && ./build/bin/find_pair --no-json "../$PDB_FILE" "../$LEGACY_INP" 2>&1 | tail -5
        cd ..
    else
        echo "  Legacy binary not found: $LEGACY_BIN"
    fi
    
    echo "Modern:"
    if [ -f "$MODERN_BIN" ]; then
        $MODERN_BIN "$PDB_FILE" "$MODERN_INP" 2>&1 | tail -5
    else
        echo "  Modern binary not found: $MODERN_BIN"
    fi
    
    # Step 2: Compare .inp files
    echo ""
    echo "--- Step 2: Comparing .inp files ---"
    if [ -f "$LEGACY_INP" ] && [ -f "$MODERN_INP" ]; then
        if diff -q "$LEGACY_INP" "$MODERN_INP" > /dev/null; then
            echo "✅ .inp files are identical"
        else
            echo "⚠️  .inp files differ:"
            diff "$LEGACY_INP" "$MODERN_INP" | head -20
        fi
    else
        echo "⚠️  Missing .inp files (legacy: $([ -f "$LEGACY_INP" ] && echo '✓' || echo '✗'), modern: $([ -f "$MODERN_INP" ] && echo '✓' || echo '✗'))"
    fi
    
    # Step 3: Run analyze and compare parameters
    echo ""
    echo "--- Step 3: Running analyze ---"
    
    LEGACY_STEP="$TEMP_DIR/${PDB}_legacy_step.par"
    LEGACY_HELIX="$TEMP_DIR/${PDB}_legacy_helix.par"
    MODERN_STEP="$TEMP_DIR/${PDB}_modern_step.par"
    MODERN_HELIX="$TEMP_DIR/${PDB}_modern_helix.par"
    
    if [ -f "$LEGACY_INP" ] && [ -f "$LEGACY_ANALYZE" ]; then
        echo "Legacy analyze:"
        cd org && ./build/bin/analyze --no-json "../$LEGACY_INP" 2>&1 | tail -3
        if [ -f "../$LEGACY_INP" ]; then
            cd .. && cp "$(dirname $LEGACY_INP)/bp_step.par" "$LEGACY_STEP" 2>/dev/null || true
            cp "$(dirname $LEGACY_INP)/bp_helical.par" "$LEGACY_HELIX" 2>/dev/null || true
        fi
        cd ..
    fi
    
    if [ -f "$MODERN_INP" ] && [ -f "$MODERN_ANALYZE" ]; then
        echo "Modern analyze:"
        $MODERN_ANALYZE "$MODERN_INP" > "$TEMP_DIR/${PDB}_modern_analyze.out" 2>&1
        # Extract step and helical parameters from output
        grep -A 100 "Step Parameters" "$TEMP_DIR/${PDB}_modern_analyze.out" | grep -v "^==" | head -50 > "$MODERN_STEP" 2>/dev/null || true
        grep -A 100 "Helical Parameters" "$TEMP_DIR/${PDB}_modern_analyze.out" | grep -v "^==" | head -50 > "$MODERN_HELIX" 2>/dev/null || true
    fi
    
    # Step 4: Compare parameters
    echo ""
    echo "--- Step 4: Comparing Step Parameters ---"
    if [ -f "$LEGACY_STEP" ] && [ -f "$MODERN_STEP" ]; then
        echo "Comparing step parameters..."
        # Would need custom comparison here
        echo "  Legacy: $([ -f "$LEGACY_STEP" ] && wc -l < "$LEGACY_STEP" || echo 0) lines"
        echo "  Modern: $([ -f "$MODERN_STEP" ] && wc -l < "$MODERN_STEP" || echo 0) lines"
    else
        echo "⚠️  Missing parameter files"
    fi
    
    echo ""
    echo "--- Step 4: Comparing Helical Parameters ---"
    if [ -f "$LEGACY_HELIX" ] && [ -f "$MODERN_HELIX" ]; then
        echo "Comparing helical parameters..."
        echo "  Legacy: $([ -f "$LEGACY_HELIX" ] && wc -l < "$LEGACY_HELIX" || echo 0) lines"
        echo "  Modern: $([ -f "$MODERN_HELIX" ] && wc -l < "$MODERN_HELIX" || echo 0) lines"
    else
        echo "⚠️  Missing parameter files"
    fi
    
    echo ""
done

echo ""
echo "========================================"
echo "Comparison complete. Files in: $TEMP_DIR"
echo "========================================"
