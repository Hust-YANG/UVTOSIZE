#!/usr/bin/env python3
"""
WSGI entry point for production deployment with Gunicorn.

Usage:
    gunicorn wsgi:app -w 4 -b 127.0.0.1:8765
"""

import sys
import os
from pathlib import Path

# Ensure the skill scripts are importable
SKILL_DIR = Path(__file__).resolve().parent.parent / ".claude" / "skills" / "UVTOSIZE" / "scripts"
sys.path.insert(0, str(SKILL_DIR))

from server import app

# Production config
app.config["DEBUG"] = False

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8765))
    app.run(host="127.0.0.1", port=port, debug=False)
