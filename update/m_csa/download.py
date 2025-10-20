#!/usr/bin/env python3
"""
Download: Download source data
"""

# Initialize

# from blang_mysql import *
from blang import Run

# Download

# https://www.ebi.ac.uk/thornton-srv/m-csa/download/

# Flat files haven't been updated since Mar/Apr 2018, need to use API downloads instead
# Run("Download", "wget -NP input 'https://www.ebi.ac.uk/thornton-srv/m-csa/media/flat_files/curated_data.csv'")

# API downloads:
# Manually curated catalytic residues in M-CSA
Run("Download", "wget -O input/residues.json 'https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json'")
# All the catalytic residues in M-CSA including residues found by homology
# "This file includes all the homologues found in Swiss-Prot and the PDB without regard to the residue conservation or function of the target sequence. This means that the residues in the homologous sequences may be performing a different catalytic function, or no catalytic function at all. Residues with 'null' fields could not be aligned to the target sequence."
# >> Not useful
# Run("Download", "wget -O input/homologues_residues.json 'https://www.ebi.ac.uk/thornton-srv/m-csa/api/homologues_residues.json'")

# Show directory
Run("Show directory", "ls -lah input")

print("Done!")
