#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

# from blang_mysql import *
from blang import Run, Waitforjobs

# Start

# All PTM sites in table 'unimod'
Run("Get PTM site contacts for PTMs from all sources", "~/scripts/qsub.sh main.py human all 1000")

# Breakdown by source
Run("Get PTM site contacts for PTMs from source UniProt", "~/scripts/qsub.sh main.py human uniprot 1000")
Run("Get PTM site contacts for PTMs from source Ochoa", "~/scripts/qsub.sh main.py human ochoa 1000")
Run("Get PTM site contacts for PTMs from source PhosphoSitePlus", "~/scripts/qsub.sh main.py human phosphositeplus 1000")
Run("Get PTM site contacts for PTMs from source dbPTM", "~/scripts/qsub.sh main.py human dbptm 1000")

Waitforjobs()

# Combine output
Run("Get header", "head -n 1 output-contacts-human-all-1000.txt > output-contacts-human-1000.txt")
Run("Concatenate output (without headers)", "cat output-contacts-human-*-1000.txt | grep -iPv '^ptm\tacc\tsite\t' >> output-contacts-human-1000.txt")
Run("Import into table", "~/scripts/import.pl output-contacts-human-1000.txt ptm_contact_aas -overwrite")

print("Done!")
