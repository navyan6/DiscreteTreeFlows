#!/usr/bin/env python3
"""Thin wrapper: temporal HIV Env split only. See prepare_hiv_splits.py."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# Re-invoke combined prep with --skip-geo
from scripts.prepare_hiv_splits import main

if __name__ == "__main__":
    sys.argv = [sys.argv[0], "--skip-geo", *sys.argv[1:]]
    main()
