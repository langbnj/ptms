#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

# from blang_mysql import *
from blang import Run

# Start

# # Download
# Run("Download", "download.py")

Run("Main", "main.py")

print("Done!")
