#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

# from blang_mysql import *
from blang import Run

# # Download
# Run("Download", "download.py")

# Start

# Run("Main", "main.py")
Run("Get only manually curated catalytic residues in M-CSA", "main.py input/residues.json")

# This here actually doesn't work: the homologous residues are in residue_sequences, but all except the first get dropped by the script. Same for residue_chains.
# Since the homolog mapping is apparently not great, I'll simply use the manually curated human residues.
# Run("Get all the catalytic residues in M-CSA including residues found by homology (apparently very roughly mapped)", "main.py input/homologues_residues.json")

print("Done!")
