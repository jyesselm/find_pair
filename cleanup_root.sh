#!/bin/bash
#
# Clean up root directory - archive docs and remove temp files
#

set -e

echo "🧹 Starting root directory cleanup..."

# Create archive directories
mkdir -p docs/validation
mkdir -p docs/archive/planning
mkdir -p docs/archive/sessions

# Move validation docs
echo "📦 Archiving validation documents..."
for file in VALIDATION_COMPLETE.md INDEX_MATCHING_PROOF.md \
            INVESTIGATION_7EH2.md RESOLUTION_7EH2.md \
            VALIDATION_TODO.md VALIDATION_IN_PROGRESS.md \
            VALIDATION_BUGS_FOUND.md WEEK1_VALIDATION_COMPLETE.md; do
    if [ -f "$file" ]; then
        mv "$file" docs/validation/
        echo "  ✅ $file → docs/validation/"
    fi
done

# Move planning docs
echo "📦 Archiving planning documents..."
for file in CLEANUP_AND_RESTRUCTURE_PLAN.md CLEANUP_FINAL_SUMMARY.md \
            CLEANUP_VERIFICATION.md REPOSITORY_CLEANUP_PLAN.md \
            README_CLEANUP.md IMPLEMENTATION_PLAN.md; do
    if [ -f "$file" ]; then
        mv "$file" docs/archive/planning/
        echo "  ✅ $file → docs/archive/planning/"
    fi
done

# Move session docs
echo "📦 Archiving session documents..."
for file in SESSION_SUMMARY.md NEXT_SESSION_QUICK_START.md \
            START_HERE.md IMMEDIATE_ACTIONS.md \
            BEFORE_AFTER_COMPARISON.md; do
    if [ -f "$file" ]; then
        mv "$file" docs/archive/sessions/
        echo "  ✅ $file → docs/archive/sessions/"
    fi
done

# Remove temporary files
echo "🗑️  Removing temporary files..."
rm -f *.inp && echo "  ✅ Removed .inp files"
rm -f *.par && echo "  ✅ Removed .par files"
rm -f ref_frames_modern.dat && echo "  ✅ Removed ref_frames_modern.dat"
rm -rf temp/ && echo "  ✅ Removed temp/ directory"

# Update .gitignore
echo "📝 Updating .gitignore..."
if ! grep -q "# Temporary test files" .gitignore 2>/dev/null; then
    cat >> .gitignore << 'EOF'

# Temporary test files
*.inp
*.out
*.par
temp/
ref_frames_*.dat
bestpairs.pdb
hel_regions.pdb
bp_order.dat
col_*.scr
EOF
    echo "  ✅ Updated .gitignore"
else
    echo "  ℹ️  .gitignore already updated"
fi

echo ""
echo "✅ Cleanup complete!"
echo ""
echo "Root directory now has:"
ls -1 *.md 2>/dev/null | wc -l | xargs echo "  Markdown files:"
echo ""
echo "Documentation organized in docs/:"
echo "  - docs/validation/ (validation reports)"
echo "  - docs/archive/planning/ (cleanup and planning)"
echo "  - docs/archive/sessions/ (session notes)"
echo ""
echo "Ready to commit:"
echo "  git add -A"
echo "  git commit -m 'Clean up root directory'"

