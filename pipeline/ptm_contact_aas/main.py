#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

discol = "dis10"        # Disorder column to use in table 'alphasa' (smoothed across ±10 aa)
surfcol = "surf"        # Core/surface column to use in table 'alphasa' (not smoothed)

(species, source, minsites) = Args(3, "[species] [ptm source in table 'unimod', or 'all' for all] [minimum number of sites for a PTM type]", "human uniprot 1000")

outfile = f"output-contacts-{species}-{source}-{minsites}.txt"

sourcestr = ""
if source != "all":
    # sourcestr = f"AND source='{source}' "
    sourcestr = f"AND m.source='{source}' "

# Start
print(f"\nGetting '{species}' PTMs from source '{source}' with at least {Comma(minsites)} sites:")
with open(outfile, "w") as out:

    # Print header
    out.write("ptm\tacc\talphacc\tafdb\tsite\tsource\tclass\ttype\tdis\tsurf\tplddt\tplddt10\taa1\taa2\tcontact\tatom1\tatom2\tdist\tpae\n")

    # Get PTM types
    for (ptm,) in Query(f"SELECT ptm FROM (SELECT ptm, COUNT(DISTINCT acc, site) AS c FROM unimod m WHERE species='{species}' {sourcestr}AND ptm IS NOT NULL GROUP BY ptm ORDER BY ptm) AS t WHERE c > {minsites} ORDER BY c DESC"):
        print(f" >> {ptm}")

        # Get amino acid
        aa = ptm[0]

        # Get accs with this PTM
        # for (acc,) in tq(Query(f"SELECT DISTINCT acc FROM unimod WHERE species='{species}' {sourcestr}AND ptm='{ptm}' ORDER BY acc")):
        # ...and ensure their uniseq sequence matches that in alphaseq (and therefore alphasa)
        # Note: I verified that adding alphasa to the query doesn't change anything, i.e. alphaseq is enough:
        # SELECT DISTINCT m.acc FROM unimod m, uniseq s, alphaseq a WHERE m.species='human' AND m.ptm='S-p' AND m.acc=s.acc AND s.type='UniProt' AND a.species='human' AND m.acc=a.acc AND s.seq=a.seq ORDER BY m.acc LIMIT 10000000;
        # SELECT DISTINCT m.acc FROM unimod m, uniseq s, alphaseq a, alphasa aa WHERE m.species='human' AND m.ptm='S-p' AND m.acc=s.acc AND s.type='UniProt' AND a.species='human' AND m.acc=a.acc AND s.seq=a.seq AND m.acc=aa.acc ORDER BY m.acc LIMIT 10000000;
        # >> Identical, going with the faster one that only uses alphaseq (~40 secs vs. ~70)
        for (acc, alphacc, afdb) in Fetch(Query(f"SELECT DISTINCT m.acc, am.map, am.afdb FROM unimod m, uniseq s, alphamap am, alphasa a WHERE m.species='{species}' {sourcestr}AND m.ptm='{ptm}' AND m.acc=s.acc AND s.type='UniProt' AND LENGTH(s.seq)>=16 AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND a.membrane IS NULL ORDER BY m.acc")):
            # print(f"   >> {acc}")

            # Get sites for this PTM
            ptmsites = FetchSet(Query(f"SELECT DISTINCT site FROM unimod m WHERE species='{species}' {sourcestr}AND ptm='{ptm}' AND acc='{acc}'"))

            # Get sites for this PTM or for any other (to avoid these sites as controls)
            modsites = FetchSet(Query(f"SELECT DISTINCT site FROM unimod m WHERE species='{species}' {sourcestr}AND acc='{acc}'"))

            # Verify that ptmsites are contained in modsites
            if not ptmsites.issubset(modsites):
                Die("Error: ptmsites are not all contained in modsites")

            # Get alphasa information (dis10, surf) for classification and classify PTM and control residues
            sites = set()
            controlsites = set()
            sitetype = {}
            siteclass = {}
            sitedis = {}
            sitesurf = {}
            for (site, dis, surf, plddt, plddt10) in Query(f"SELECT site, {discol}, {surfcol}, plddt, plddt10 FROM alphasa WHERE acc='{alphacc}' AND afdb='{afdb}' AND aa='{aa}' AND membrane IS NULL"):
                # print(f"\n     >> Residue types")
                # print(f"       >> {site}")
                # print(f"         >> {discol}: {dis}")
                # print(f"         >> {surfcol}: {surf}")
                
                # # Skip if it's a membrane residue
                # if membrane == 1:
                #     Log(f"membrane residue skipped for acc|site (skipped)", f"{acc}|{site}")
                #     continue

                if site in ptmsites:
                    # PTM site
                    type = "mod"
                elif site in modsites:
                    # Other modified site (skip)
                    continue
                else:
                    # Control site
                    type = "control"
                    # Add site to control sites
                    controlsites.add(site)

                # Add site to sites of interest (PTM or control, but not other modified sites)
                sites.add(site)

                # PTM or control?
                sitetype[site] = type

                # Classify site
                # if surf == "C":
                #     type += "_core"
                # elif dis == ".":
                #     type += "_strsurf"
                # else:
                #     type += "_dissurf"
                if surf == "C":
                    myclass = "core"
                elif dis == ".":
                    myclass = "strsurf"
                else:
                    myclass = "dissurf"
                
                # Store site type
                sitedis[site] = dis
                sitesurf[site] = surf
                sitetype[site] = type
                siteclass[site] = myclass
                # print(f"           >> {type}")

                # Insert into table 'unimod_control'
                
            # Verify that ptmsites and controlsites don't overlap
            if not ptmsites.isdisjoint(controlsites):
                Die("Error: ptmsites and controlsites overlap")
            
            # Get contacts
            for site in sites:
                # print(f"\n     >> Contacts")
                # print(f"       >> {site}")
                # print(f"         >> {sitetype[site]}")
                # print(f"           >> {surfcol}: {surf}")

                Log(f"successfully wrote contact to outfile for acc", acc)
                Log(f"successfully wrote contact to outfile for acc|site", f"{acc}|{site}")
                Log(f"successfully wrote contact to outfile for '{sitetype[site]}' acc|site", f"{acc}|{site}")

                # Get contacts
                for (site1, site2, aa1, aa2, type, atom1, atom2, dist, pae) in Query(f"SELECT site1, site2, aa1, aa2, type, atom1, atom2, dist, pae FROM alphacon WHERE acc='{alphacc}' AND afdb='{afdb}' AND (site1='{site}' OR site2='{site}')"):
                    # Write to output file
                    if site1 == site:
                        out.write(f"{ptm}\t{acc}\t{alphacc}\t{afdb}\t{site}\t{source}\t{siteclass[site]}\t{sitetype[site]}\t{sitedis[site]}\t{sitesurf[site]}\t{plddt}\t{plddt10}\t{aa1}\t{aa2}\t{type}\t{atom1}\t{atom2}\t{dist}\t{pae}\n")
                    elif site2 == site:
                        out.write(f"{ptm}\t{acc}\t{alphacc}\t{afdb}\t{site}\t{source}\t{siteclass[site]}\t{sitetype[site]}\t{sitedis[site]}\t{sitesurf[site]}\t{plddt}\t{plddt10}\t{aa2}\t{aa1}\t{type}\t{atom2}\t{atom1}\t{dist}\t{pae}\n")

    

Show(lim=0)

print(f"Wrote to '{outfile}'")

print("\nDone!")
