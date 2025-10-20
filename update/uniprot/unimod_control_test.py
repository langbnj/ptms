#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# from scipy.stats import wilcoxon
# from scipy.stats import mannwhitneyu
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

species = "human"
minsites = 1000

Args(0, "", "")

# Start

# Initialize final lists
total = []
ptm = []
ptm_accs = []
control = []
control_accs = []


# Get PTM residues
# print(f"Getting accs that contain ptm IS NOT NULL PTM sites:")
# for acc, in Fetch(Query(f"SELECT m.acc FROM unimod m, uniseq s WHERE m.species='{species}' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='UniProt' GROUP BY m.acc ORDER BY m.acc")):
print(f"Getting accs that contain 'ptm IS NOT NULL' PTM sites with >{minsites} sites, and that have evolutionary rates in table 'evorate_rate4site_einsi_tree_1para':")
# for acc, in Fetch(Query(f"SELECT m.acc FROM unimod m, uniens ue, evorate_rate4site_einsi_tree_1para r WHERE m.species='{species}' AND m.ptm IS NOT NULL AND m.acc=ue.acc AND ue.ensp=r.ensp GROUP BY m.acc ORDER BY m.acc")):
for acc, in Fetch(Query(f"SELECT m.acc FROM unimod m, uniens ue WHERE m.species='{species}' AND m.ptm IN (SELECT ptm FROM (SELECT ptm, COUNT(DISTINCT acc, site) AS c FROM unimod WHERE species='{species}' AND ptm IS NOT NULL GROUP BY ptm ORDER BY ptm) AS t WHERE c > {minsites}) AND m.acc=ue.acc AND ue.ensp IN (SELECT DISTINCT r.ensp FROM evorate_rate4site_einsi_tree_1para r) GROUP BY m.acc ORDER BY m.acc")):
    # print(f" >> {acc}")
    
    # Get sequence
    # seq = FetchOne(Query(f"SELECT DISTINCT seq FROM uniseq WHERE type='UniProt' AND acc='{acc}'"))
    seq = FetchOne(Query(f"SELECT DISTINCT seq FROM uniseq WHERE type IN ('UniProt', 'UniIso') AND acc='{acc}'"))
    # print(f"   >> Sequence: \t{seq}")

    # Get PTM positions in this acc (all)
    ptmpos = sorted(list(FetchSet(Query(f"SELECT m.site FROM unimod m WHERE m.species='{species}' AND m.acc='{acc}' AND m.aa IN ('A', 'C', 'D', 'E', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'Y')"))))
    # print(f"   >> PTM positions: \t{' '.join([str(i) for i in ptmpos])}")
    
    # Get positions of residues of interest (ACDEKMNPQRSTY)
    respos = [i + 1 for i in range(len(seq)) if seq[i] in "ACDEKMNPQRSTY"]
    
    # Remove PTM residues from controlpos set
    controlpos = sorted(list(set(respos) - set(ptmpos)))
    # print(f"   >> Control positions: \t{' '.join([str(i) for i in controlpos])}")
    
    # Add to final lists
    total += [(acc, pos) for pos in respos]
    control += [(acc, pos) for pos in controlpos]
    control_accs += [acc]
    ptm += [(acc, pos) for pos in ptmpos]
    ptm_accs += [acc]


# Print unique control list length
print(f" >> Unique total residues: {Comma(len(set(total)))}")
print(f" >> Unique PTM residues: {Comma(len(set(ptm)))}")
print(f" >> Unique PTM accs: {Comma(len(set(ptm_accs)))}")
print(f" >> Unique control residues: {Comma(len(set(control)))}")
print(f" >> Unique control accs: {Comma(len(set(control_accs)))}")

Show(lim=0)
# ~/pipeline/evolutionary_rate_analysis/tmp/tmp-dataframe-all-human-all-lichtarge_einsi_tree_1para-AlphaFold-coresurf.txt
# u = read_tsv("unimod_control.accsites.txt")
# r = read_tsv("tmp-dataframe-all-human-all-rate4site_einsi_tree_1para-AlphaFold-coresurf.accsites.txt")
# su = set(u.acc + '|' + u.site.astype(str))

evorate_control_df = read_tsv("../../pipeline/evolutionary_rate_analysis/tmp/tmp-dataframe-all-human-all-rate4site_einsi_tree_1para-AlphaFold-coresurf.txt")
evorate_control = [(acc, pos) for acc, pos in zip(evorate_control_df.acc, evorate_control_df.site)]
evorate_control_accs = set(evorate_control_df.acc)

# Retain only controls that are also in evorate_control_accs
# evorate_control = [(acc, pos) for acc, pos in zip(evorate_control_df.acc, evorate_control_df.site) if acc in control_accs]
control_filtered_with_evorates = [(acc, pos) for acc, pos in control if acc in evorate_control_accs]

# Test if evorate_control and control_filtered_with_evorates are identical sets
if set(evorate_control) == set(control_filtered_with_evorates):
    print(" >> evorate_control and control_filtered_with_evorates are identical sets")
else:
    print(" >> evorate_control and control_filtered_with_evorates are NOT identical sets")
d()

print("\nDone!")
