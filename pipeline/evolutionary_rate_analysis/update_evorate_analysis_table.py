#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

table = "evorate_analysis"

(species, resamples) = Args(2, "[species] [resamples]", "human 10000")
outfile = f"output-evorate_analysis-{resamples}-{species}.tsv"

# Remove existing outfile
Run("Remove existing outfile", f"rm -f {outfile}")

# Start
print(f"\nReading 'output/output-pvalues-{resamples}-{species}-*.tsv' and writing to '{outfile}':")

# Get all output/output-pvalues-{resamples}-{species}-*.tsv files
infiles = glob.glob(f"output/output-pvalues-{resamples}-{species}-*.tsv")
for infile in tq(nsort(infiles)):
    
    print(f"\n >> '{infile}'") if Switch('debug') else None

    # e.g. output/output-pvalues-10000-human-all-capra1_ginsi-AlphaFold-coresurf.tsv
    m = rx(r"^output/output-pvalues-(\d+)-(\w+)-(\w+)-([^_]+)_(\w+)-(\w+)-(\w+)\.tsv$", infile)
    if not m:
        Die(f"Couldn't parse '{infile}'")
    else:
        (resamples, species, source, evorate, mafftmode, predictor, mode) = m

    # Read file
    df = pd.read_csv(infile, sep="\t")

    # Add columns parsed from file name at start of data frame
    df.insert(0, "species", species)
    df.insert(1, "predictor", predictor)
    df.insert(2, "evorate", evorate)
    df.insert(3, "mafftmode", mafftmode)
    df.insert(4, "source", source)
    df.insert(5, "mode", mode)

    # Move resamples column towards start
    df.insert(6, "resamples", df.pop("resamples"))

    # Replace *** with 3 in all columns ending in sig
    for col in df.columns:
        if col.endswith("sig"):
            df[col] = df[col].replace("***", 3)
            df[col] = df[col].replace("**", 2)
            df[col] = df[col].replace("*", 1)
            df[col] = df[col].replace("(*)", -1)
            df[col] = df[col].replace("(**)", -2)
            df[col] = df[col].replace("(***)", -3)
            # Replace blank with 0
            df[col] = df[col].replace("", 0)

    # Append to outfile (TSV)
    if not Exists(outfile):
        # Create file
        df.to_csv(outfile, sep="\t", index=False, header=True)
    else:
        # Append to file
        df.to_csv(outfile, sep="\t", index=False, header=False, mode="a")
    

Show(lim=0)

# Import into MySQL
Run(f"Import into MySQL table '{table}'", f"~/scripts/import.pl {outfile} {table} -overwrite -allindices")

print("\nDone!")
