#!/usr/bin/env python3
"""Entry point for cWW miss annotator CLI.

This script can be run directly from the prototypes/pair_identification directory.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from cww_miss_annotator.cli import main

if __name__ == "__main__":
    sys.exit(main())
