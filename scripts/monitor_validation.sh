#!/bin/bash
# Monitor ls_fitting validation progress

LOG_FILE="data/validation_results/ls_fitting_fixed_validation.log"

while true; do
    clear
    echo "=== LS_FITTING VALIDATION MONITOR ==="
    echo "Time: $(date)"
    echo ""
    
    if [ ! -f "$LOG_FILE" ]; then
        echo "Waiting for log file to be created..."
        sleep 5
        continue
    fi
    
    # Show last progress line
    echo "Latest progress:"
    tail -20 "$LOG_FILE" | grep -E "Progress|SUMMARY" | tail -5
    echo ""
    
    # Count lines in log
    lines=$(wc -l < "$LOG_FILE")
    echo "Log file has $lines lines"
    
    # Check if done
    if grep -q "SUMMARY" "$LOG_FILE"; then
        echo ""
        echo "=== VALIDATION COMPLETE ==="
        tail -30 "$LOG_FILE"
        break
    fi
    
    sleep 10
done

