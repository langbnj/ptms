#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
import pandas as pd
import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

table = "m_csa"

# Args(0, "", "")
infile = Args(1, "[input file]", "input/residues.json")

# infile = f"input/residues.json"
tmpfile = f"tmp-{table}.txt"

# Start
print(f"\nReading '{infile}' and restructuring data frame...")

# Load JSON file
t = pd.read_json(infile)

# First entry:
# d()
# t.iloc[0]

# mcsa_id
# A numeric ID for an enzyme. There can be multiple residues (rows) per enzyme.
# >> Keep

# roles_summary
# t["roles_summary"].value_counts()
# roles_summary
# metal ligand                                                                                                                                     937
# electrostatic stabiliser                                                                                                                         922
# electrostatic stabiliser, hydrogen bond donor                                                                                                    369
# proton acceptor, proton donor                                                                                                                    210
# proton shuttle (general acid/base)                                                                                                               205
#                                                                                                                                                 ... 
# hydrogen bond acceptor, hydrogen bond donor, increase electrophilicity, increase nucleophilicity, metal ligand, proton acceptor, proton donor      1
# activator, electrostatic stabiliser, hydrogen bond donor, proton acceptor, proton donor, proton relay                                              1
# activator, increase electrophilicity, polar interaction, steric role                                                                               1
# electrostatic stabiliser, hydrogen bond donor, polar interaction, proton acceptor, proton donor                                                    1
# electrostatic interaction, electrostatic stabiliser, transition state stabiliser                                                                   1
# Name: count, Length: 549, dtype: int64
# >> Simple string
# >> Seems good to keep
# >> Keep

# function_location_abv
# t["function_location_abv"].value_counts()
# function_location_abv
#           4730
# main-N     365
# main-C     110
# ptm         21
# N-term      15
# main         4
# C-term       3
# Name: count, dtype: int64
# >> Not sure what this means, and it's usually blank
# >> Skip

# main_annotation
# t["main_annotation"].value_counts()
# main_annotation
# Acts as a general acid/base.                                                                                                                  124
#                                                                                                                                                81
# Forms part of the magnesium binding site.                                                                                                      76
# Forms part of the zinc binding site.                                                                                                           33
# Binds one of the Zn(II) ions.                                                                                                                  23
#                                                                                                                                              ... 
# Sabilises the cationic transition state with oxocarbenium ion character                                                                         1
# Asp272 is proposed to hydrogen bond with the pyrimidine ring, stabilising the substrate in a linear alignment with the Cys113 nucleophile.      1
# Acts as a general acid/base, activating the nucleophilic cysteine residue.                                                                      1
# Binds the carboxylate group of the substrate, helps maintain the substrate in the correct orientation for the reaction to occur.                1
# Acts as a nucleophile and deprotonates AdoHcy in the last step of catalysis.                                                                    1
# Name: count, Length: 3365, dtype: int64
# >> Simple string
# >> Seems good to keep
# >> Keep

# ptm
# PTMs:
# t["ptm"].value_counts()
# ptm
#        5221
# Llp       4
# Pyr       3
# Cso       3
# Kcx       3
# Csd       2
# Sep       2
# Css       2
# Trq       1
# Fgl       1
# SE7       1
# Tpq       1
# Ocs       1
# Kxc       1
# Nep       1
# Smc       1
# Name: count, dtype: int64
# e.g. https://pdbj.org/emnavi/quick.php?id=Sep (phosphoserine). Only occurs twice, so PTMs are not too relevant. The others are rare modifications where there isn't any mass spec data.

# roles
# >> Too complex
# >> Skip

# residue_chains
# PDB annotation
# >> I first assumed this was the main reference (as shown below for 'residue_sequences', some entries don't have uniprot_ids). However:
# >> Weirdly, some entries don't have residue_chains:
# t["residue_chains"].apply(len).value_counts()
# residue_chains
# 1    5201
# 0      47
# Name: count, dtype: int64
# >> This means PDB can't be the main reference either.
# List of entries without residue_chains:
# t[t["residue_chains"].apply(len) == 0]
#       mcsa_id                                      roles_summary function_location_abv                                    main_annotation ptm                                              roles residue_chains                                  residue_sequences
# 1698      277  activator, covalently attached, electrofuge, e...                        The cysteine "pairs" act as nucleophiles to re...      [{'group_function': 'electrostatic interaction...             []  [{'uniprot_id': 'P00392', 'code': 'Cys', 'is_r...
# 1699      277   electrostatic stabiliser, hydrogen bond acceptor                        Help stabilise and hold in place the reactive ...      [{'group_function': 'electrostatic interaction...             []  [{'uniprot_id': 'P00392', 'code': 'Tyr', 'is_r...
# [...]
# >> All of these have uniprot_ids.
# >> ...so the main reference can either be uniprot_id or PDB.
# >> Let's see if residue_sequences 'is_reference' is used to distinguish:
# t[t["residue_chains"].apply(len) == 0]["residue_sequences"].apply(lambda x: x[0]["is_reference"]).value_counts()
# residue_sequences
# True    47
# Name: count, dtype: int64
# >> No. Even for these 47 that do not have a PDB reference, the residue_sequences 'is_reference' is True.
# Actually, is_reference is always True (as for residue_sequences below)
# t[t["residue_chains"].apply(len) > 0]["residue_chains"].apply(lambda x: x[0]["is_reference"]).value_counts()
# residue_chains
# True    5201
# Name: count, dtype: int64
# >> is_reference is useless (never False).

# residue_sequences
# There is always a single residue_sequence:
# t["residue_sequences"].apply(len).unique()
# array([1])
# 
# Each one is a list of length 1, containing a dict:
# 
# t["residue_sequences"].apply(lambda x: x[0])
# 0       {'uniprot_id': 'P56868', 'code': 'Asp', 'is_re...
# 1       {'uniprot_id': 'P56868', 'code': 'Cys', 'is_re...
# 2       {'uniprot_id': 'P56868', 'code': 'Cys', 'is_re...
# 3       {'uniprot_id': 'P56868', 'code': 'His', 'is_re...
# 4       {'uniprot_id': 'P56868', 'code': 'Ser', 'is_re...
#                               ...                        
# 5243    {'uniprot_id': '', 'code': 'Lys', 'is_referenc...
# 5244    {'uniprot_id': '', 'code': 'Lys', 'is_referenc...
# 5245    {'uniprot_id': '', 'code': 'Gly', 'is_referenc...
# 5246    {'uniprot_id': '', 'code': 'Glu', 'is_referenc...
# 5247    {'uniprot_id': '', 'code': 'Glu', 'is_referenc...
# Name: residue_sequences, Length: 5248, dtype: object
# 
#  t["residue_sequences"].apply(lambda x: x[0])[0]
# {'uniprot_id': 'P56868', 'code': 'Asp', 'is_reference': True, 'resid': 7}
# 
# There are always 4 entries:
# t["residue_sequences"].apply(lambda x: x[0]).apply(len).unique()
#
# t["residue_sequences"].apply(lambda x: x[0].keys()).value_counts()
# residue_sequences
# (uniprot_id, code, is_reference, resid)    5248
# Name: count, dtype: int64
# >> They're always "uniprot_id, code, is_reference, resid"
# >> There are no other sequence IDs than uniprot_id.
# 
# t["residue_sequences"].apply(lambda x: x[0]["uniprot_id"]).value_counts()
# residue_sequences
#           279
# P81701     27
# P0ABI8     23
# P00968     18
# P33247     17
#          ... 
# Q2RSB4      1
# Q56308      1
# P15807      1
# O66186      1
# Q9ZFQ5      1
# >> Sometimes, the uniprot_id is blank.
# Is there a PDB reference in these cases?
# t[t["residue_sequences"].apply(lambda x: x[0]["uniprot_id"]) == ""]
#       mcsa_id                                      roles_summary function_location_abv                                    main_annotation ptm                                              roles                                     residue_chains                                  residue_sequences
# 1333      219                      proton acceptor, proton donor                        His258 operates as the general base and acid t...      [{'group_function': 'proton shuttle (general a...  [{'chain_name': 'A', 'pdb_id': '4kxv', 'assemb...  [{'uniprot_id': '', 'code': 'His', 'is_referen...
# 1336      219                      proton acceptor, proton donor                        Involve in acid-base catalysis by accepting an...      [{'group_function': 'proton shuttle (general a...  [{'chain_name': 'A', 'pdb_id': '4kxv', 'assemb...  [{'uniprot_id': '', 'code': 'Lys', 'is_referen...
# 2257      383                                hydrogen bond donor                        IGP is held by a network of hydrogen bonds for...      [{'group_function': '', 'function_type': 'inte...  [{'chain_name': 'A', 'pdb_id': '1a50', 'assemb...  [{'uniprot_id': '', 'code': 'Ser', 'is_referen...
# 2258      383                                hydrogen bond donor                        IGP is held by a network of hydrogen bonds for...      [{'group_function': '', 'function_type': 'inte...  [{'chain_name': 'A', 'pdb_id': '1a50', 'assemb...  [{'uniprot_id': '', 'code': 'Ser', 'is_referen...
# 2259      383                                hydrogen bond donor                        Tyr102 and Thr183 stabilise Asp60 via hydrogen...      [{'group_function': '', 'function_type': 'inte...  [{'chain_name': 'A', 'pdb_id': '1a50', 'assemb...  [{'uniprot_id': '', 'code': 'Thr', 'is_referen...
# ...       ...                                                ...                   ...                                                ...  ..                                                ...                                                ...                                                ...
# 5243     1004  electrostatic stabiliser, transition state sta...                                 Form a catalytic motif: lysine triangle.      [{'group_function': 'electrostatic interaction...  [{'chain_name': 'A', 'pdb_id': '3lkk', 'assemb...  [{'uniprot_id': '', 'code': 'Lys', 'is_referen...
# 5244     1004  electrostatic stabiliser, transition state sta...                                 Form a catalytic motif: lysine triangle.      [{'group_function': 'electrostatic interaction...  [{'chain_name': 'A', 'pdb_id': '3lkk', 'assemb...  [{'uniprot_id': '', 'code': 'Lys', 'is_referen...
# 5245     1004   hydrogen bond donor, transition state stabiliser                main-N     Gly8 and Lys14 stabilise the transition state.      [{'group_function': 'electrostatic interaction...  [{'chain_name': 'A', 'pdb_id': '3lkk', 'assemb...  [{'uniprot_id': '', 'code': 'Gly', 'is_referen...
# 5246     1005                                    proton acceptor                        Acts as a nucleophile and deprotonates Noradre...      [{'group_function': 'proton shuttle (general a...  [{'chain_name': 'A', 'pdb_id': '1hnn', 'assemb...  [{'uniprot_id': '', 'code': 'Glu', 'is_referen...
# 5247     1005                                    proton acceptor                        Acts as a nucleophile and deprotonates AdoHcy ...      [{'group_function': 'proton shuttle (general a...  [{'chain_name': 'A', 'pdb_id': '1hnn', 'assemb...  [{'uniprot_id': '', 'code': 'Glu', 'is_referen...
# 
# [279 rows x 8 columns]
# >> 279 such cases...let's see if they have a PDB reference:
# t[t["residue_sequences"].apply(lambda x: x[0]["uniprot_id"]) == ""]["residue_chains"].apply(len).value_counts()
# residue_chains
# 1    278
# 0      1
# >> In all except one case, yes, they have a PDB reference.
# >> I have no clue how there can be an entry without either a PDB or a UniProt reference. Here it is:
# t[(t["residue_sequences"].apply(lambda x: x[0]["uniprot_id"]) == "") & (t["residue_chains"].apply(len) == 0)]
#       mcsa_id roles_summary function_location_abv main_annotation ptm                                              roles residue_chains                                  residue_sequences
# 4061      744   nucleophile                                            [{'group_function': 'covalent catalysis', 'fun...             []  [{'uniprot_id': '', 'code': 'Thr', 'is_referen...
# t[(t["residue_sequences"].apply(lambda x: x[0]["uniprot_id"]) == "") & (t["residue_chains"].apply(len) == 0)]["residue_chains"]
# 4061    []
# Name: residue_chains, dtype: object
# t[(t["residue_sequences"].apply(lambda x: x[0]["uniprot_id"]) == "") & (t["residue_chains"].apply(len) == 0)]["residue_sequences"].iloc[0]
# [{'uniprot_id': '', 'code': 'Thr', 'is_reference': True, 'resid': 264}]
# >> Weird. Can ignore this single case, though.
# Also, all of the 279 cases where uniprot_id is blank still have "is_reference" as True:
# t[t["residue_sequences"].apply(lambda x: x[0]["uniprot_id"]) == ""]["residue_sequences"].apply(lambda x: x[0]["is_reference"]).value_counts()
# residue_sequences
# True    279
# Name: count, dtype: int64
# Actually, every single UniProt entry has is_reference as True, so it's not useful:
# t["residue_sequences"].apply(lambda x: x[0]["is_reference"]).value_counts()
# residue_sequences
# True    5248
# Name: count, dtype: int64

# Drop unneeded columns:
t = t.drop(columns=["function_location_abv", "roles"])


# Keep only the first list item in residue_chains (and None if the list is empty):
t["residue_chains"] = t["residue_chains"].apply(lambda x: x[0] if len(x) > 0 else None)

# Keep only the first list item in residue_sequences:
t["residue_sequences"] = t["residue_sequences"].apply(lambda x: x[0])


# Convert residue_sequences to pd.Series and append columns
t = t.join(t["residue_sequences"].apply(pd.Series))

# Convert residue_chains to pd.Series and append columns
t = t.join(t["residue_chains"].apply(pd.Series), rsuffix="_pdb")

# Drop original columns
t = t.drop(columns=["residue_chains", "residue_sequences"])

# Replace '' with NaN for assembly
t["assembly"] = t["assembly"].apply(lambda x: np.nan if x == '' else x)

# Convert float back to int for assembly
t["assembly"] = t["assembly"].apply(lambda x: pd.Series(x, dtype = pd.Int64Dtype()))

# Convert float back to int for assembly, resid and auth_resid (leaving NaNs untouched)
t["resid_pdb"] = t["resid_pdb"].apply(lambda x: pd.Series(x, dtype = pd.Int64Dtype()))
t["auth_resid"] = t["auth_resid"].apply(lambda x: pd.Series(x, dtype = pd.Int64Dtype()))



# Drop unneeded columns
# 'mcsa_id', 'roles_summary', 'main_annotation', 'ptm', 'uniprot_id',
#        'code', 'is_reference', 'resid', 'chain_name', 'pdb_id',
#        'assembly_chain_name', 'assembly', 'code_pdb', 'resid_pdb',
#        'auth_resid', 'is_reference_pdb', 'domain_name', 'domain_cath_id'
t = t.drop(columns=["is_reference", "is_reference_pdb"])

# Rename columns
# "uniprot_id" to acc
t = t.rename(columns={"uniprot_id": "acc"})
# "code" to aa3
t = t.rename(columns={"code": "aa3"})
# "code_pdb" to aa3_pdb
t = t.rename(columns={"code_pdb": "aa3_pdb"})
# "resid" to site
t = t.rename(columns={"resid": "site"})
# "resid_pdb" to site_pdb
t = t.rename(columns={"resid_pdb": "site_pdb"})

# # Add single-letter amino acid column
t["aa"] = t["aa3"].apply(lambda x: ThreeToOne(x))
t["aa_pdb"] = t["aa3_pdb"].apply(lambda x: ThreeToOne(x) if x not in ('X', np.nan) else "")

# Everything except selenocysteine (U) is a standard amino acid in the "aa" column
# t[t["aa"].apply(Aa) == False]["aa"]
# 2629    U
# 3260    U
# 3265    U
# 4514    U
# Name: aa, dtype: object
# >> OK!

# Everything except selenocysteine (U) and unknown (X) is a standard amino acid in the "aa_pdb" column
# t[t["aa_pdb"].apply(lambda x: Aa(x) if not pd.isna(x) else x) == False]["aa_pdb"]
# 323     X
# 2418    X
# 2629    U
# 3260    U
# 3265    U
# 3319    X
# 3682    X
# 5195    X
# >> OK!

# Escape \t and \n in main_annotation and remove \r
# e.g. t["main_annotation"].apply(lambda x: "\r" in x).value_counts()
# Version with homologues doesn't have "main_annotation" column
if "main_annotation" in t.columns:
    t["main_annotation"] = t["main_annotation"].apply(lambda x: x.replace("\t", "\\t").replace("\n", "\\n").replace("\r", ""))

# Filtering 

# Verify acc using table 'uniacc' and get name, canon, and species
print(f"\nFiltering out accs that aren't found in table 'uniacc'...")
removed = 0
for acc in tq(t["acc"].unique()):
    if Numrows(Query(f"SELECT name, canon, species FROM uniacc WHERE acc='{acc}'")) == 0:
        Log(f"m-csa acc not found in table 'uniacc' for acc (skipped)", acc)
        # Remove this acc from data frame
        t = t[t["acc"] != acc]
        removed += 1
    else:
        for name, canon, species in Query(f"SELECT name, canon, species FROM uniacc WHERE acc='{acc}'"):
            # Add name to data frame
            t.loc[t["acc"] == acc, "name"] = name
            # Add canon to data frame
            t.loc[t["acc"] == acc, "canon"] = canon
            # Add species to data frame
            t.loc[t["acc"] == acc, "species"] = species
print(f" >> Removed {removed} accs")

# Verify that residue is correct using table 'uniseq'
# Loop over entire data frame
removed = 0
print(f"\nFiltering out residue mismatches in table 'uniseq'...")
for i, row in tq(t.iterrows(), total=len(t)):
    # Get acc, site, aa
    acc = row["acc"]
    site = row["site"]
    aa = row["aa"]
    # Get residue from table 'uniseq'
    for (residue,) in Query(f"SELECT SUBSTRING(seq, {site}, 1) FROM uniseq WHERE acc='{acc}' AND type IN ('UniProt', 'UniIso')"):
        # If residue is not correct
        if residue != aa:
            Log(f"m-csa residue mismatch in table 'uniseq' for species|acc|site|expected|found (skipped)", f"{species}|{acc}|{int(site)}|{aa}|{residue}")
            Log(f"m-csa residue mismatch in table 'uniseq' for acc|site|expected|found (skipped)", f"{acc}|{int(site)}|{aa}|{residue}")
            Log(f"m-csa residue mismatch in table 'uniseq' for species (skipped)", species)
            # Remove this row from data frame
            t = t.drop(i)
            removed += 1
print(f" >> Removed {removed} rows")

# Reorder data frame, putting acc, site and aa first
# Version with homologues doesn't have "main_annotation" column
if "main_annotation" in t.columns:
    t = t[["name", "acc", "species", "site", "aa", "aa3", "ptm", "mcsa_id", "pdb_id", "site_pdb", "chain_name", "assembly_chain_name", "assembly", "aa_pdb", "aa3_pdb", "domain_name", "domain_cath_id", "roles_summary", "main_annotation"]]
else:
    t = t[["name", "acc", "species", "site", "aa", "aa3", "ptm", "mcsa_id", "pdb_id", "site_pdb", "chain_name", "assembly_chain_name", "assembly", "aa_pdb", "aa3_pdb", "domain_name", "domain_cath_id", "roles_summary"]]

# Write to TSV for import
print(f"\nWriting to '{tmpfile}'...")
t.to_csv(tmpfile, sep="\t", index=False)

# Import into MySQL
Run(f"Import into table '{table}'", f"~/scripts/import.pl {tmpfile} {table} -overwrite -allindices")

# Reorder table by name, acc, species, site
print(f"\nReordering table by mcsa_id, species, name, acc, site...")
Query(f"ALTER TABLE {table} ORDER BY mcsa_id, species, name, acc, site")

print(f"Resetting id column in table '{table}'...")
Query(f"ALTER TABLE {table} DROP COLUMN id, DROP PRIMARY KEY, ADD COLUMN id int unsigned NOT NULL AUTO_INCREMENT FIRST, ADD PRIMARY KEY (id)")

Show(lim=0)

print("\nDone!")
