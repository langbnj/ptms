#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# from scipy.stats import wilcoxon
# Import Mann-Whitney U test and t-test
from scipy.stats import mannwhitneyu, ttest_ind
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

# table = "table"
species = "human"

resamples = 10

Args(0, "-noptms: Don't use residues that have any known PTM sites as control residues (doesn't make any conceptual sense since we're looking at enzyme active sites here, though)", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

tmpnoptms = ""
if Switch('noptms'):
    tmpnoptms = "-noptms"

# Start
print(f"\nGetting M-CSA catalytic residues from table 'm_csa':")

# Read results pickles if they exist
if Exists(f"tmp-t{tmpnoptms}.txt") and Exists(f"tmp-res{tmpnoptms}.txt"):
    print(" >> Found existing results pickles")
    # Read results from pickles
    # t = pd.read_pickle("tmp-t.pkl")
    # res = pd.read_pickle("tmp-res.pkl")
    # Read results from TSVs
    t = pd.read_csv(f"tmp-t{tmpnoptms}.txt", sep="\t")
    res = pd.read_csv(f"tmp-res{tmpnoptms}.txt", sep="\t")
    print(" >> Read data frames")
else:
    print(" >> No results pickles found, running")

    # Initialise results data frame (which will be made up of tcsa and tcon)
    t = pd.DataFrame()
    res = pd.DataFrame()
    # for method in ("capra0", "lichtarge", "rate4site"):
    for method in ("rate4site", "lichtarge", "capra0", "norm_rate4site", "norm_lichtarge", "norm_capra0"):
        for mafftmode in ("einsi_tree_1para", "einsi_tree", "ginsi_tree", "linsi_tree", "einsi", "ginsi", "linsi"):
            evorate = f"{method}_{mafftmode}"
            # Get CSA residues
            # tcsa = Panda(f"SELECT DISTINCT '{evorate}' AS evorate, 'csa' AS type, c.acc, c.site, c.aa, r.rate FROM m_csa c, uniens ue, evorate_{evorate} r WHERE c.species='{species}' AND c.acc=ue.acc AND ue.ensp=r.ensp AND c.site=r.site")
            tcsa = Panda(f"SELECT '{evorate}' AS evorate, 'csa' AS type, c.acc, c.site, c.aa, MIN(r.rate) AS rate FROM m_csa c, uniens ue, evorate_{evorate} r WHERE c.species='{species}' AND c.acc=ue.acc AND ue.ensp=r.ensp AND c.site=r.site GROUP BY c.acc, c.site")
            # Add tcsa to results data frame
            t = pd.concat([t, tcsa], ignore_index=True)
            print(f" >> Getting control residues for {evorate}:")
            for (i, e) in tq(tcsa.iterrows(), total=len(tcsa)):
                (evorate, type, acc, site, aa, rate) = e
                
                # Get sequence and disorder prediction
                # seq = FetchOne(Query(f"SELECT seq FROM uniseq WHERE type IN ('UniProt', 'UniIso') AND acc='{acc}'"))
                seq = FetchOne(Query(f"SELECT seq FROM uniseq WHERE type IN ('UniProt') AND acc='{acc}'"))
                query = Query(f"SELECT seq FROM uniseq WHERE type='AlphaFold' AND acc='{acc}'")
                if Numrows(query) == 0:
                    Log("no type='AlphaFold' sequence in uniseq for acc (skipped)", acc)
                    continue
                disseq = FetchOne(query)
                
                # Get type='CoreSurf' for filtering out membrane residues
                query = Query(f"SELECT seq FROM uniseq WHERE type='CoreSurf' AND acc='{acc}'")
                if Numrows(query) == 0:
                    Log("no type='CoreSurf' sequence in uniseq for acc (skipped)", acc)
                    continue
                surfseq = FetchOne(query)
                
                # Check if seq and disseq are the same length
                if len(seq) != len(disseq) or len(seq) != len(surfseq):
                    Log("seq and disseq lengths do not match for acc (skipped)", acc)
                    continue
                
                # Check if seq and surfseq are the same length
                if len(seq) != len(surfseq):
                    Log("seq and surfseq lengths do not match for acc (skipped)", acc)
                    continue
                    
                # Skip if surfseq is 'M' at site
                if surfseq[int(site) - 1] == "M":
                    Log("surfseq is 'M' (membrane) at site for acc (skipped)", acc)
                    continue
                
                # Check aa at site
                if seq[int(site) - 1] != aa:
                    Die(f"aa at site {site} does not match {aa} for acc '{acc}'")
                    # continue
                
                # Get positions in seq that match aa
                aapos = [i + 1 for i in range(len(seq)) if seq[i] == aa]
                
                # Get positions in disseq that are structured (.)
                strpos = [i + 1 for i in range(len(seq)) if disseq[i] == "."]
                
                # Get positions in surfseq that are membrane (M)
                mempos = [i + 1 for i in range(len(seq)) if surfseq[i] == "M"]
                
                # Get intersection of aapos and strpos using set .intersect
                if Switch('structured_controls_only'):
                    controlpos = sorted(list(set(aapos) & set(strpos)))
                else:
                    controlpos = sorted(list(set(aapos)))

                # Remove mempos from controlpos
                if Switch('no_membrane_controls'):
                    controlpos = sorted(list(set(controlpos) - set(mempos)))
                
                # Get M-CSA residue positions (to avoid)
                csapos = FetchSet(Query(f"SELECT site FROM m_csa WHERE species='{species}' AND acc='{acc}'"))
                
                # Remove M-CSA residues from controlpos set
                controlpos = sorted(list(set(controlpos) - set(csapos)))
                
                # Avoid PTM sites
                if Switch('noptms'):
                    # Get positions that have any PTMs in table unimod
                    # ptmpos = FetchSet(Query(f"SELECT DISTINCT site FROM unimod WHERE species='{species}' AND acc='{acc}' AND aa='{aa}' AND ptm IS NOT NULL"))
                    ptmpos = FetchSet(Query(f"SELECT DISTINCT site FROM unimod WHERE species='{species}' AND acc='{acc}' AND aa='{aa}'"))
                    # Remove PTM positions from controlpos set
                    controlpos = sorted(list(set(controlpos) - set(ptmpos)))
                
                # Get matching control residues (same protein, same aa as the original CSA residue, structured, and no other PTM sites (if switch -noptms enabled))
                # tcon = Panda(f"SELECT DISTINCT '{evorate}' AS evorate, 'control' AS type, '{acc}' AS acc, r.site, '{aa}' AS aa, r.rate FROM uniens ue, evorate_{evorate} r WHERE ue.species='{species}' AND ue.acc='{acc}' AND r.site IN ('" + "', '".join(str(x) for x in controlpos) + "') AND ue.ensp=r.ensp")
                tcon = Panda(f"SELECT '{evorate}' AS evorate, 'control' AS type, '{acc}' AS acc, r.site, '{aa}' AS aa, MIN(r.rate) AS rate FROM uniens ue, evorate_{evorate} r WHERE ue.species='{species}' AND ue.acc='{acc}' AND r.site IN ('" + "', '".join(str(x) for x in controlpos) + "') AND ue.ensp=r.ensp GROUP BY ue.acc, r.site")
                if len(tcon) == 0:
                    Log("no control residues for acc|site (skipped)", f"{acc}|{site}")
                    continue
                # Append to results data frame
                t = pd.concat([t, tcon], ignore_index=True)
                
                # Compare tcon to tcsa
                # tcsa[tcsa.acc == acc].rate.mean()
                # tcsa[tcsa.acc == acc].rate.median()
                # mannwhitneyu(tcsa[tcsa.acc == acc].rate, tcon.rate).pvalue
                
                # Calculate normalised ranks
                ranks_df = pd.concat([pd.DataFrame({"type": "csa", "rate": tcsa[tcsa.acc == acc].rate}),
                                      pd.DataFrame({"type": "control", "rate": tcon.rate})])
                ranks_df["nrank"] = (ranks_df.rate.rank() - 1) / (len(ranks_df) - 1)
                mean_nrank_csa = ranks_df[ranks_df.type == "csa"].nrank.mean()
                mean_nrank_control = ranks_df[ranks_df.type == "control"].nrank.mean()
                mean_nrank = mean_nrank_csa - mean_nrank_control
                median_nrank_csa = ranks_df[ranks_df.type == "csa"].nrank.median()
                median_nrank_control = ranks_df[ranks_df.type == "control"].nrank.median()
                median_nrank = median_nrank_csa - median_nrank_control

                res = pd.concat([res, pd.DataFrame({"evorate": evorate, "acc": acc, "site": site, "aa": aa,
                                                    "mean_csa": tcsa[tcsa.acc == acc].rate.mean(),
                                                    "mean_control": tcon.rate.mean(),
                                                    "meandif": tcsa[tcsa.acc == acc].rate.mean() - tcon.rate.mean(),
                                                    "mean_zscore": mean((tcsa[tcsa.acc == acc].rate - tcon.rate.mean()) / tcon.rate.std()),
                                                    "mean_normalised_rank": mean_nrank,
                                                    "median_csa": tcsa[tcsa.acc == acc].rate.median(),
                                                    "median_control": tcon.rate.median(),
                                                    "mediandif": tcsa[tcsa.acc == acc].rate.median() - tcon.rate.median(),
                                                    "median_zscore": median((tcsa[tcsa.acc == acc].rate - tcon.rate.mean()) / tcon.rate.std()),
                                                    "median_normalised_rank": median_nrank,
                                                    "pvalue": mannwhitneyu(tcsa[tcsa.acc == acc].rate, tcon.rate, alternative="less").pvalue,
                                                    "pvalue_twotailed": mannwhitneyu(tcsa[tcsa.acc == acc].rate, tcon.rate).pvalue,
                                                    "pvalue_ttest": ttest_ind(tcsa[tcsa.acc == acc].rate, tcon.rate, alternative="less").pvalue,
                                                    "pvalue_ttest_twotailed": ttest_ind(tcsa[tcsa.acc == acc].rate, tcon.rate).pvalue,
                                                    }, index = [i])])
                # res = pd.concat([res, pd.DataFrame({"evorate": evorate, "acc": acc, "site": site, "aa": aa,
                #                                     "meandif": tcsa[tcsa.acc == acc].rate.mean() - tcon.rate.mean(),
                #                                     "mediandif": tcsa[tcsa.acc == acc].rate.median() - tcon.rate.median()})], ignore_index=True)
                # d()

    # # Save data frames to pickles
    # t.to_pickle("tmp-t.pkl")
    # res.to_pickle("tmp-res.pkl")
    
    # Save data frames to TSVs
    print(" >> Wrote data frames")
    t.to_csv(f"tmp-t{tmpnoptms}.txt", sep="\t")
    res.to_csv(f"tmp-res{tmpnoptms}.txt", sep="\t")

# # Resampling analysis
# Would be ultra slow, especially if not using an evorate|acc-based dict of lists
# for evorate in nsort(res.evorate.unique()):
#     print(f" >> {evorate}")
#     for acc in nsort(res.acc.unique()):
#         print(f"   >> {acc}")
#         for aa in res[res.acc == acc].aa.unique():
#             print(f"\n{aa}:")

# d()
resout = res.groupby('evorate').agg(
                                 mean_pvalue = ('pvalue', 'mean'), 
                                 median_pvalue = ('pvalue', 'median'),
                                 mean_pvalue_twotailed = ('pvalue_twotailed', 'mean'),
                                 median_pvalue_twotailed = ('pvalue_twotailed', 'median'),
                                 mean_pvalue_ttest = ('pvalue_ttest', 'mean'), 
                                 median_pvalue_ttest = ('pvalue_ttest', 'median'),
                                 mean_pvalue_ttest_twotailed = ('pvalue_ttest_twotailed', 'mean'),
                                 median_pvalue_ttest_twotailed = ('pvalue_ttest_twotailed', 'median'),
                                 mean_zscore = ('mean_zscore', 'mean'),
                                 median_zscore = ('median_zscore', 'median'),
                                 mean_normalised_rank = ('mean_normalised_rank', 'mean'),
                                 median_normalised_rank = ('median_normalised_rank', 'median')
                                 ).sort_values('median_pvalue')

print(resout)

# Write resout to TSV file
resout.to_csv(f"output-table{tmpnoptms}.tsv", sep="\t")

Show(lim=0)

print("\nDone!")
