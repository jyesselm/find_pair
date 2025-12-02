#!/bin/bash
# Quick script to check validation progress

CSV="data/index_validation_status.csv"

if [ ! -f "$CSV" ]; then
    echo "CSV file not found: $CSV"
    exit 1
fi

TOTAL=$(tail -n +2 "$CSV" | wc -l | tr -d ' ')
PASS=$(grep ",PASS," "$CSV" | wc -l | tr -d ' ')
FAIL=$(grep ",FAIL," "$CSV" | wc -l | tr -d ' ')
SKIP=$(grep ",SKIP," "$CSV" | wc -l | tr -d ' ')

echo "=== Index Validation Progress ==="
echo "Total tested: $TOTAL"
echo "✅ PASS: $PASS"
echo "❌ FAIL: $FAIL"
echo "⏭️  SKIP: $SKIP"
echo ""

if [ "$FAIL" -gt 0 ]; then
    echo "Failed PDBs:"
    grep ",FAIL," "$CSV" | cut -d',' -f1
fi

echo ""
echo "Last 5 tested:"
tail -5 "$CSV" | cut -d',' -f1,6 | sed 's/,/ -> /'

