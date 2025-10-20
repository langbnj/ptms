#!/usr/bin/env python3 -u
"""
Adds 'CoreSurf' core/surface calls to 'uniseq' table based on 'alphasa' table. Mostly for compatibility with older scripts.

Note: Membrane residues will be set to 'M' (membrane) in the core/surface sequence. This causes them to be ignored by some of my analyses.
Note: Membrane residues come from alphasa "membrane" column, and ultimately from UniProt's 'unifeat' table.
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

table = "uniseq"
type = "CoreSurf"
alphamap_type = 'uniprot'
alphamap_version = '2022_04'

column = Args(1, "[column: e.g. dis/dis10] [-humanonly] [-iso] [-plddt70] [-plddt90] [-pae2] [-pae1]\n\n -humanonly: Only human proteins (much faster)\n -iso: In addition to 'UniProt' canonical sequences, also include 'UniIso' isoform sequences", "dis10 -humanonly")

uniseq_types = "('UniProt')"
if Switch('iso'):
    uniseq_types = "('UniProt', 'UniIso')"

# mem = ' '
mem = 'M'

if Switch('plddt70'):
    type = "CoreSurf_pLDDT70"
if Switch('plddt90'):
    type = "CoreSurf_pLDDT90"

if Switch('pae2'):
    type += "_PAE2"
elif Switch('pae1'):
    type += "_PAE1"

if not Switch('debug'):
    print(f"Clearing '{type}' entries from table '{table}'")
    query = Query(f"DELETE FROM {table} WHERE type='{type}'")
    print(f" >> {Comma(Numrows(query))} rows affected")


# Start
print(f"\nGetting 'uniseq' proteins with a matching 'alphaseq' sequence (using table 'alphamap') and inserting their core/surface calls (based on '{column}') into table '{table}' as type '{type}':")

if (Switch('humanonly')):
    # Human only
    # Use table 'alphamap' and best=1
    # # AlphaPTM main protein list (no isoforms, no short peptides <16 aa, no fragmented proteins ≥2700 aa, no BZUX residues):
    # SELECT DISTINCT s.acc FROM uniseq s, alphamap m WHERE s.species='human' AND s.type='UniProt' AND s.acc=s.canon AND s.seq NOT REGEXP '[BZUX]' AND LENGTH(s.seq)>=16 AND LENGTH(s.seq)<=2699 AND s.acc=m.value AND m.type='{alphamap_type}' AND m.version='{alphamap_version}' AND m.best=1 AND m.map IS NOT NULL;
    # AlphaPTM main protein list (no isoforms, no short peptides <16 aa):
    accquery = Query(f"SELECT DISTINCT s.name, s.acc, s.canon, s.species, s.fragments, s.precursor, s.trembl, s.seq, m.map, m.afdb FROM uniseq s, alphamap m WHERE s.species='human' AND s.type IN {uniseq_types} AND s.acc=s.canon AND LENGTH(s.seq)>=16 AND s.acc=m.value AND m.type='{alphamap_type}' AND m.version='{alphamap_version}' AND m.best=1 AND m.map IS NOT NULL ORDER BY s.name, s.acc")
else:
    # All species
    # Use table 'alphamap' and best=1
    accquery = Query(f"SELECT DISTINCT s.name, s.acc, s.canon, s.species, s.fragments, s.precursor, s.trembl, s.seq, m.map, m.afdb FROM uniseq s, alphamap m WHERE s.type IN {uniseq_types} AND s.acc=m.value AND m.type='{alphamap_type}' AND m.version='{alphamap_version}' AND m.map IS NOT NULL AND m.best=1 ORDER BY s.name, s.acc")

for name, acc, canon, species, fragments, precursor, trembl, seq, map, afdb in tq(accquery, total=Numrows(accquery)):

    # Get residue-level core/surface predictions
    # query = Query(f"SELECT site, {column}, membrane, plddt FROM alphasa WHERE acc='{map}' AND afdb='{afdb}' ORDER BY site", Switch('debug'))
    query = Query(f"SELECT a.site, a.{column}, a.membrane, a.plddt, IFNULL(MIN(c.pae), 100) AS minpae FROM alphasa a LEFT OUTER JOIN alphacon c ON a.acc=c.acc AND a.afdb=c.afdb AND (c.site1=a.site OR c.site2=a.site) WHERE a.acc='{map}' AND a.afdb='{afdb}' GROUP BY a.site ORDER BY a.site", Switch('debug'))

    # Verify that it's the correct number of residues (otherwise there must be two different structures in alphasa with different core/surface predictions)
    if (Numrows(query) != len(seq)):
        Die(f"Error: Sequence length for uniseq acc '{acc}' is {len(seq)}, but there are {Numrows(query)} residues for it in table 'alphasa'")

    # Assemble core/surface sequence from individual-residue query
    disseq = ''
    prevsite = 0
    for (site, dis, membrane, plddt, minpae) in query:
        if site != prevsite + 1:
            Die(f"Error: Sites aren't sequential in query results ({site} follows {prevsite})")
        prevsite = site

        # pLDDT filtering
        if Switch('plddt70') and plddt < 70:
            dis = ' '
        if Switch('plddt90') and plddt < 90:
            dis = ' '
        # PAE filtering
        if Switch('pae2') and minpae > 2:
            dis = ' '
        if Switch('pae1') and minpae > 1:
            dis = ' '

        # Membrane filtering
        if membrane == '1' and dis != ' ':
            # # Membrane: set to C (core)
            # disseq += 'C'
            # Membrane: set to M
            disseq += mem
        else:
            disseq += dis

    # Insert
    q = f"INSERT INTO {table} SET name='{name}', acc='{acc}', canon='{canon}', species='{species}', type='{type}', fragments='{fragments}', precursor='{precursor}', trembl='{trembl}', seq='{disseq}'"
    q = q.replace("='None'", "=NULL")
    if not Switch('debug'):
        Query(q)
    else:
        State(q)

    Log(f"Inserted '{type}' core/surface sequence into '{table}' for canon|acc|alphamap acc", f"{canon}|{acc}|{map}")
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for canon|acc", f"{canon}|{acc}")
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for acc|alphamap acc", f"{acc}|{map}")
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for alphamap acc", map)
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for acc", acc)
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for afdb='{afdb}' acc", acc)
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for canon", canon)
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for name", name)
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for species", species)
    Log(f"Inserted '{type}' core/surface sequence into '{table}' for sequence", seq)

Show()

if not Switch('debug'):
    Optimize(table)

print("\nDone!")
