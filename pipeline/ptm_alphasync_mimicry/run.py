#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

# from blang_mysql import *
from blang import Run, Waitforjobs

# Start

# Get mimicry cases for all negatively charged PTMs
# (...and K-ac, which might overlap with K-mal and K-suc, and hence K-ac sites could alternatively perhaps also be K-mal or K-suc sites)
for ptm in ('S-p', 'T-p', 'Y-p', 'K-ac', 'K-mal', 'K-suc'):

    # Negatively charged PTM originating from D or E (-) contacting K or R (+) (einsi_tree_1para)
    Run("Get any mimicry cases", f"~/scripts/qsub.sh main.py human einsi_tree_1para {ptm} DE KR rate4site 4 root 10 1 0.9 1 0.9 -debug")
    Run("Get primate-specific mimicry cases", f"~/scripts/qsub.sh main.py human einsi_tree_1para {ptm} DE KR rate4site 4 primates 10 1 0.9 1 0.9 -debug")

Waitforjobs()

print("Done!")
