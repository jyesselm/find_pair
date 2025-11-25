#!/bin/bash
# Script to add debug output to C++ code for validation mismatches
# This will help trace exactly what's happening in good_hb_atoms()

echo "To add debug output, we need to:"
echo "1. Add conditional debug logging to good_hb_atoms()"
echo "2. Log atom names, patterns, idx values, and return values"
echo "3. Enable via environment variable or compile flag"
echo ""
echo "Suggested approach:"
echo "- Add #ifdef DEBUG_GOOD_HB_ATOMS blocks"
echo "- Or use std::cerr with conditional compilation"
echo "- Log to a file for specific pairs: 1VBY(20,21) and 3AVY(1204,1223)"

