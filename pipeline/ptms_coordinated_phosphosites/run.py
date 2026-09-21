#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

import subprocess
import sys

# phospho_reference_asa.py is not run here: its results are already written into
# find_coordinated_phosphoresidues.py as constants. Run it on its own to
# reproduce them.
STEPS = [
    ("Search the PDB for phosphoresidues coordinated by Lys or Arg", "find_coordinated_phosphoresidues.py"),
    ("Collapse the contacts to unique phosphosites", "summarise_phosphosites.py"),
    ("Label kinase and activation-loop sites", "classify_kinase_sites.py"),
    ("Build Supplementary Table 6b", "build_supplementary_table.py"),
]

# Start

for description, script in STEPS:
    print(f"\n{description}:\n >> {script}\n")
    subprocess.run([sys.executable, script], check=True)

print("\nDone!")
