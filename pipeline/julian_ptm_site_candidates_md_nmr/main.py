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

Args(0, "", "")

# Extra proteins to include (because they have interesting mimicry results, and filamins)
include_list = (
    "HEBP2_HUMAN",
    "FLNA_HUMAN",
    "FLNB_HUMAN",
    "FLNC_HUMAN",
)

outfile = "output/candidate_ptm_sites.tsv"

# Minimum size of the domain the modified residue is found in (to avoid short peptides)
min_domain_size = 40

# Maximum size of the domain the modified residue is found in
# max_domain_size = 100
max_domain_size = 200


# # Maximum number of variants in gnomAD tolerated (to ensure the residue is functionally important)
# max_variant_metric = "ac"
# max_variant = 1
# Allow rare variants (<1e-4, i.e. <0.01%, as in gnomAD PTM SNP figures)
max_variant_metric = "af"
max_variant = 0.0001


# Minimum pLDDT score for the modified residue
min_plddt = 70
# min_plddt = 50
# min_plddt = 0

# Minimum pLDDT score for surrounding residues
min_plddt10 = 70
# min_plddt10 = 50
# min_plddt10 = 0

# UniProt version (for table 'alphamap')
uniprot_release = "2022_04"

# Start
# - The residue must be buried (but close to accessible residues) according to AlphaFold. The criterion I used is that the relative accessible surface area needs to be ≤25% for the residue itself, but >25% within a ±10 aa window.
print(f"""
Getting candidate PTM sites:
- I've removed the requirement for the protein to have a "clinically significant" variant in ClinVar since it was too restrictive. There are likely to be plenty of biologically interesting proteins with sufficient literature on their function that do not yet have known disease variants. I have, however, put proteins with clinically significant disease variants at the top of the table in bold. Also, PTM sites that we previously decided to exclude as candidates are in italics.
- The residue must be buried (but close to accessible residues) according to AlphaFold. The criterion I used is that the relative accessible surface area needs to be ≤25% for the residue itself, but that there must be a residue with RSA >25% within a ±5 aa window. The distance to the nearest "surface" residue is given in the "surface_distance" column. This is actually between 1 and 3 aa for all of the sites in the table.
- The residue must be resolved in at least one PDB structure of ≤{max_domain_size} aa (NMR or X-ray) to ensure the protein expresses well. It also needs to be ≥{min_domain_size} aa, rather than a small fragment.
- The residue's structural position must be predicted reasonably confidently by AlphaFold (pLDDT ≥{min_plddt}).
- The residue's structural surroundings must be predicted reasonably confidently by AlphaFold (pLDDT ±10 aa ≥{min_plddt10}).
- The residue must be in a structured domain of ≤{max_domain_size} aa according to AlphaFold. It also needs to be ≥{min_domain_size} aa to ensure it is an actual independently folded domain.
- The residue must not have any known natural variants, with the exception of {max_variant_metric}<={max_variant} (rare variants), according to gnomAD v4 (~800,000 individuals). This is to ensure the residue is functionally important.
- The residue must be well-conserved across species (to ensure it is functionally important) according to a conservation score (Mayrose et al.'s Rate4Site evolutionary rate score across Ensembl Compara one-to-one orthologs (or best-matching orthologs if none available) ≤ 0.3).
- We do not tolerate cases with missing evolutionary rate data.
- As indicated by the new "category" column, the sites are sorted with clinically significant proteins at the top (especially the two sites we selected for HECBioSim), followed by the evolutionarily interesting HEBP2_HUMAN S37 phosphosite, followed by sites we previously decided to exclude, and finally sites in proteins that do not have clinically significant disease variants (which could still be biologically very interesting). Within these categories, the sites are sorted alphabetically by the protein's UniProt ID and the residue number.
""")

with open(outfile, "w") as out:

    # Write header
    # out.write(f"name\tacc\tsite\tptm\trelasa\tsurf\tsurf10\tdis\tdis10\tsec\tevorate\tpdbs\tclinsig\tclingen\tclinorigin\tclinphen\n")
    # out.write(f"uniprot_name\tgene_symbol\tprotein_name\tuniprot_accession\tuniprot_link\talphafold_link\talphasync_link\tresidue_number\tptm\tptm_source_count\tptm_sources\tproline_upstream\tproline_downstream\tsurface_distance\trelative_accessible_surface_area\tsurface_category\tsurface_category_smoothed_10\tdisorder_category\tdisorder_category_smoothed_10\talphafold_plddt\talphafold_plddt_smoothed_10\tdssp_secondary_structure\tevolutionary_rate\tpdbs\tclinically_significant\tclinical_genetic_test\tclinvar_residues\tclinvar_mutation_origin\tclinical_phenotypes\n")
    out.write(f"UniProt name\tGene symbol\tProtein name\tUniProt accession\tUniProt link\tAlphaFold DB link\tAlphaSync link\tResidue number\tPTM\tPTM sources\tPTM source list\tProline upstream\tProline downstream\tSurface distance\tRelative solvent-accessible surface area\tBuried/surface category\tBuried/surface category smoothed ± 10 aa\tDisorder category\tDisorder category smoothed ± 10 aa\tAlphaFold pLDDT score\tAlphaFold pLDDT score smoothed ± 10 aa\tDSSP secondary structure\tEvolutionary rate\tPDB IDs\tProtein is clinically significant (ClinVar)\tProtein has genetic test (ClinVar)\tClinVar clinically significant residues\tClinVar clinically significant mutation origin\tClinVar clinically significant phenotypes\n")
    # out.write(f"UniProt name\tGene symbol\tProtein name\tUniProt accession\tUniProt link\tAlphaFold DB link\tAlphaSync link\tResidue number\tPTM\tPTM sources\tPTM source list\tProline upstream\tProline downstream\tSurface distance\tRelative solvent-accessible surface area\tBuried/surface category\tBuried/surface category smoothed ± 10 aa\tDisorder category\tDisorder category smoothed ± 10 aa\tAlphaFold pLDDT score\tAlphaFold pLDDT score smoothed ± 10 aa\tDSSP secondary structure\tEvolutionary rate\tPDB IDs\tSite ClinVar mutation origin\tSite ClinVar clinical phenotypes\tProtein is clinically significant (ClinVar)\tProtein has genetic test (ClinVar)\tProtein ClinVar residues\tProtein clinically significant ClinVar residues\n")

    # Ochoa or UniProt PTM sites (reliable - also, phosphosites will be easier to manage experimentally than other PTMs)
    # More complex query (but too slow)
    # for (name, acc, site, ptm) in tq(Query(f"SELECT DISTINCT m.name, m.acc, m.site, m.ptm FROM unimod m, alphasa a, snps_clinvar c, unipdb p WHERE m.acc=a.acc AND a.surf='C' AND a.surf10='S' AND a.dis='.' AND a.dis10='*' AND m.species='human' AND m.source IN ('Ochoa', 'UniProt') AND m.ptm IS NOT NULL AND m.acc=c.acc AND c.clinsig=1 AND p.acc=m.acc AND m.site BETWEEN p.start AND p.stop")):
    # # alphasa: no check (check is too slow)
    # for (name, acc, site, ptm) in tq(Query(f"SELECT DISTINCT name, acc, site, ptm FROM unimod WHERE species='human' AND ptm IS NOT NULL AND source IN ('Ochoa', 'UniProt') ORDER BY name, acc, site, ptm")):
    # # alphasa: must be buried (only just) and structured (only just)
    # for (name, acc, site, ptm) in tq(Query(f"SELECT DISTINCT m.name, m.acc, m.site, m.ptm FROM unimod m, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.source IN ('Ochoa', 'UniProt') AND m.acc=a.acc AND m.site=a.site AND a.surf='C' AND a.surf10='S' AND a.dis='.' AND a.dis10='*' ORDER BY m.name, m.acc, m.site, m.ptm")):
    
    if Switch('debug2'):
        # with pre-specified protein names (much faster, but will be incorrect if filter criteria change)
        # mainquery = Query(f"SELECT DISTINCT m.name, m.acc, m.site, m.ptm, p.symbol, p.fullname FROM unimod m, uniprot p, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.source IN ('Ochoa', 'UniProt') AND p.acc=m.acc AND m.acc=a.acc AND m.site=a.site AND a.surf='C' AND a.surf10='S' AND m.name IN ('EMD_HUMAN', 'GABR2_HUMAN', 'HNRPK_HUMAN', 'HUWE1_HUMAN', 'PAK1_HUMAN', 'PAK2_HUMAN', 'SASH1_HUMAN', 'SNP29_HUMAN', 'TY3H_HUMAN', 'UBA5_HUMAN', 'WASP_HUMAN', 'XIAP_HUMAN', 'ZAP70_HUMAN', 'HEBP2_HUMAN') ORDER BY m.name, m.acc, m.site, m.ptm")
        # mainquery = Query(f"SELECT m.name, m.acc, m.site, m.ptm, p.symbol, p.fullname, COUNT(DISTINCT source) AS sources, GROUP_CONCAT(DISTINCT REPLACE(REPLACE(m.subset, '_small', ' small-scale'), '_large', ' large-scale') ORDER BY subset SEPARATOR ', ') FROM unimod m, uniprot p, alphamap am, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.source IN ('Ochoa', 'UniProt') AND p.acc=m.acc AND m.acc=am.value AND am.type='uniprot' AND am.version='{uniprot_release}' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND a.relasa<0.26 AND a.relasa10>=0.26 AND m.name IN ('EMD_HUMAN', 'GABR2_HUMAN', 'HNRPK_HUMAN', 'HUWE1_HUMAN', 'PAK1_HUMAN', 'PAK2_HUMAN', 'SASH1_HUMAN', 'SNP29_HUMAN', 'TY3H_HUMAN', 'UBA5_HUMAN', 'WASP_HUMAN', 'XIAP_HUMAN', 'ZAP70_HUMAN', 'HEBP2_HUMAN') GROUP BY m.name, m.acc, m.site, m.ptm, p.symbol, p.fullname ORDER BY m.name, m.acc, m.site, m.ptm")
        mainquery = Query(f"SELECT m.name, m.acc, m.site, m.ptm, p.symbol, p.fullname, COUNT(DISTINCT source) AS sources, GROUP_CONCAT(DISTINCT REPLACE(REPLACE(m.subset, '_small', ' small-scale'), '_large', ' large-scale') ORDER BY subset SEPARATOR ', ') FROM unimod m, uniprot p, alphamap am, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND p.acc=m.acc AND m.acc=am.value AND am.type='uniprot' AND am.version='{uniprot_release}' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND a.relasa<0.26 AND a.relasa10>=0.26 AND m.name IN ('EMD_HUMAN', 'GABR2_HUMAN', 'HNRPK_HUMAN', 'HUWE1_HUMAN', 'PAK1_HUMAN', 'PAK2_HUMAN', 'SASH1_HUMAN', 'SNP29_HUMAN', 'TY3H_HUMAN', 'UBA5_HUMAN', 'WASP_HUMAN', 'XIAP_HUMAN', 'ZAP70_HUMAN', 'HEBP2_HUMAN') GROUP BY m.name, m.acc, m.site, m.ptm, p.symbol, p.fullname ORDER BY m.name, m.acc, m.site, m.ptm")
    else:
        # alphasa: must be buried (only just)
        # surf='C' AND surf10='S'
        # mainquery = Query(f"SELECT DISTINCT m.name, m.acc, m.site, m.ptm, p.symbol, p.fullname FROM unimod m, uniprot p, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.source IN ('Ochoa', 'UniProt') AND p.acc=m.acc AND m.acc=a.acc AND m.site=a.site AND a.surf='C' AND a.surf10='S' ORDER BY m.name, m.acc, m.site, m.ptm")
        # surf='C' (and filter below to make sure there's a surface residue within ±5 aa)
        # mainquery = Query(f"SELECT DISTINCT m.name, m.acc, am.map, m.site, m.ptm, p.symbol, p.fullname FROM unimod m, uniprot p, alphamap am, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.source IN ('Ochoa', 'UniProt') AND p.acc=m.acc AND m.acc=am.value AND am.type='uniprot' AND am.version='{uniprot_release}' AND am.best=1 AND am.map=a.acc AND m.site=a.site AND (a.surf='C' OR m.name IN ('" + "', '".join(include_list) + "')) ORDER BY m.name, m.acc, m.site, m.ptm")
        # mainquery = Query(f"SELECT m.name, m.acc, am.map, m.site, m.ptm, p.symbol, p.fullname, COUNT(DISTINCT source) AS sources, GROUP_CONCAT(DISTINCT REPLACE(REPLACE(m.subset, '_small', ' small-scale'), '_large', ' large-scale') ORDER BY subset SEPARATOR ', ') AS source_subset FROM unimod m, uniprot p, alphamap am, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.source IN ('Ochoa', 'UniProt') AND p.acc=m.acc AND m.acc=am.value AND am.type='uniprot' AND am.version='{uniprot_release}' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND (a.surf='C' OR m.name IN ('" + "', '".join(include_list) + "')) GROUP BY m.name, m.acc, am.map, m.site, m.ptm, p.symbol, p.fullname ORDER BY m.name, m.acc, m.site, m.ptm")
        mainquery = Query(f"SELECT m.name, m.acc, am.map, m.site, m.ptm, p.symbol, p.fullname, COUNT(DISTINCT source) AS sources, GROUP_CONCAT(DISTINCT REPLACE(REPLACE(m.subset, '_small', ' small-scale'), '_large', ' large-scale') ORDER BY subset SEPARATOR ', ') AS source_subset FROM unimod m, uniprot p, alphamap am, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND p.acc=m.acc AND m.acc=am.value AND am.type='uniprot' AND am.version='{uniprot_release}' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND (a.surf='C' OR m.name IN ('" + "', '".join(include_list) + "')) GROUP BY m.name, m.acc, am.map, m.site, m.ptm, p.symbol, p.fullname ORDER BY m.name, m.acc, m.site, m.ptm")
        # without ORDER BY (faster)
        # mainquery = Query(f"SELECT m.name, m.acc, m.site, m.ptm, p.symbol, p.fullname FROM unimod m, uniprot p, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.source IN ('Ochoa', 'UniProt') AND p.acc=m.acc AND m.acc=a.acc AND m.site=a.site AND a.surf='C' LIMIT 1")

    for (name, acc, afdb_acc, site, ptm, symbol, fullname, sources, source_subset) in tq(mainquery):

        print(f"\n >> {name} >> {symbol} >> {fullname} >> {acc} >> {afdb_acc} >> {site} >> {ptm} >> {sources} sources >> {source_subset}") if Switch('debug') else None

        # Buried, but only just
        # Using surf10='S':
        # alphasa a WHERE m.acc=a.acc AND m.site=a.site AND a.surf='C' AND a.surf10='S' AND a.dis='.' AND a.dis10='*'
        # alphasa = Panda(f"SELECT * FROM alphasa WHERE acc='{acc}' AND site='{site}' AND surf='C' AND surf10='S'")
        # exception for HEBP2, which is 26% exposed
        # alphasa = Panda(f"SELECT * FROM alphasa WHERE acc='{acc}' AND site='{site}' AND relasa<0.26 AND relasa10>=0.26")
        # Using surf='S' within ±3 aa:
        # surfdist = Panda(f"SELECT MIN(ABS(site - {site})) AS surfdist FROM alphasa WHERE acc='{acc}' AND surf='S' AND site BETWEEN {site - 5} AND {site + 5}")
        surfdist = int(FetchOne(Query(f"SELECT MIN(ABS(site - {site})) AS surfdist FROM alphasa WHERE acc='{afdb_acc}' AND surf='S'")))
        # if len(alphasa) > 1:
        #     Die(f"Error: Multiple matches in table alphasa for acc|site '{acc}|{site}'")
        if surfdist > 5 and name not in include_list:
            print(f" >> too buried >> SKIP") if Switch('debug') else None
            # Log("skipped ptm site because it is not 'buried (but close to a surface region)' for acc|site (skipped)", f"{acc}|{site}")
            Log("skipped ptm site because it is >5 aa away from the nearest surface residue for acc|site (skipped)", f"{acc}|{site}")
            continue

        # Get alphasa data on the PTM site itself
        # alphasa = Panda(f"SELECT * FROM alphasa WHERE acc='{acc}' AND site='{site}' AND surf='C'")
        alphasa = Panda(f"SELECT * FROM alphasa WHERE acc='{afdb_acc}' AND site='{site}'")

        # Filter by pLDDT (AlphaFold confidence)
        if alphasa['plddt'][0] < min_plddt and name not in include_list:
            print(f" >> pLDDT too low at <{min_plddt} >> SKIP")
            Log(f"skipped ptm site because its plddt is <{min_plddt} for acc|site (skipped)", f"{acc}|{site}")
            continue

        # Filter by pLDDT10 (AlphaFold confidence) (smoothed)
        if alphasa['plddt10'][0] < min_plddt10 and name not in include_list:
            print(f" >> pLDDT10 too low >> SKIP")
            Log(f"skipped ptm site because its plddt10 is <{min_plddt10} for acc|site (skipped)", f"{acc}|{site}")
            continue

        # alphasa = Panda(f"SELECT * FROM alphasa WHERE acc='{acc}' AND site='{site}' AND dis='.' AND dis10='*'")
        # if len(alphasa) == 0:
        #     print(f" >> not structured >> SKIP")
        #     Log("skipped ptm site because it is not 'structured (but close to a disordered region)' for acc|site (skipped)", f"{acc}|{site}")
        #     continue

        # Disease mutations (within the protein)
        # snps_clinvar c WHERE m.acc=c.acc AND c.clinsig=1
        snps_clinvar = Panda(f"SELECT site, original, variant, clinsig, genetictest, phenotype, originsimple FROM snps_clinvar WHERE acc='{acc}' AND clinsig=1")
        # # Don't require clinsig=1
        # snps_clinvar = Panda(f"SELECT site, original, variant, clinsig, genetictest, phenotype, originsimple FROM snps_clinvar WHERE acc='{acc}'")
        # if len(snps_clinvar) == 0:
        # if len(snps_clinvar) == 0 and name not in include_list:
        #     print(f" >> not clinsig >> SKIP") if Switch('debug') else None
        #     Log("skipped ptm site because there are no snps_clinvar clinsig=1 disease mutations in this protein for acc|site (skipped)", f"{acc}|{site}")
        #     continue

        # Other NMR and/or X-ray structures of the region would be great (so we know it's tractable)
        unipdb = Panda(f"SELECT DISTINCT pdb, stop - start + 1 AS len FROM unipdb WHERE acc='{acc}' AND {site} BETWEEN start AND stop")
        if len(unipdb) == 0 and name not in include_list:
            print(f" >> no pdb >> SKIP") if Switch('debug') else None
            Log("skipped ptm site because there are no pdb structures of it for acc|site (skipped)", f"{acc}|{site}")
            continue
        
        # Small domain of ≤max_domain_size aa (for a small MD box and for NMR)

        # Check PDB structure sizes
        unipdb = Panda(f"SELECT DISTINCT pdb, stop - start + 1 AS len FROM unipdb WHERE acc='{acc}' AND {site} BETWEEN start AND stop HAVING len <= {max_domain_size}")
        if len(unipdb) == 0 and name not in include_list:
            print(f" >> pdb too large >> SKIP") if Switch('debug') else None
            Log(f"skipped ptm site because all its pdb structures are >{max_domain_size} aa for acc|site (skipped)", f"{acc}|{site}")
            continue
        if name in include_list:
            unipdb = Panda(f"SELECT DISTINCT pdb, stop - start + 1 AS len FROM unipdb WHERE acc='{acc}' AND {site} BETWEEN start AND stop")

        # ≥40 aa (not just a tiny peptide)

        # Check PDB structure sizes
        unipdb = Panda(f"SELECT DISTINCT pdb, stop - start + 1 AS len FROM unipdb WHERE acc='{acc}' AND {site} BETWEEN start AND stop HAVING len >= {min_domain_size}")
        if len(unipdb) == 0 and name not in include_list:
            print(f" >> pdb too small >> SKIP") if Switch('debug') else None
            Log(f"skipped ptm site because all its pdb structures are <{min_domain_size} aa for acc|site (skipped)", f"{acc}|{site}")
            continue
        if name in include_list:
            unipdb = Panda(f"SELECT DISTINCT pdb, stop - start + 1 AS len FROM unipdb WHERE acc='{acc}' AND {site} BETWEEN start AND stop")



        # Check AlphaFold domain size
        query = Query(f"SELECT seq FROM uniseq WHERE acc='{acc}' AND type='AlphaFold'")
        if Numrows(query) == 0:
            print(f" >> no AlphaFold sequence in uniseq for acc (must be a sequence mismatch between UniProt and alphaseq) >> SKIP") if Switch('debug') else None
            Log("skipped ptm site because there is no AlphaFold disorder string in uniseq for acc|site (skipped)", f"{acc}|{site}")
            continue
        (alphaseq) = FetchOne(query)
        # Before the site:
        before = alphaseq[:site-1]
        # After the site:
        after = alphaseq[site:]
        # How far does the (structured) domain extend to the left?
        left = len(before) - len(before.rstrip('.'))
        # How far does the (structured) domain extend to the right?
        right = len(after) - len(after.lstrip('.'))
        # How long is the structured domain?
        domainlen = left + right + 1

        # Check AlphaFold domain size (≤ max_domain_size aa)
        if domainlen > max_domain_size and name not in include_list:
            print(f" >> domain too large >> SKIP") if Switch('debug') else None
            Log(f"skipped ptm site because its AlphaFold alphaseq structured domain is >{max_domain_size} aa for acc|site (skipped)", f"{acc}|{site}")
            continue
        if name in include_list:
            print(f" >> domain size {domainlen} (exception for include_list) (kept)")

        # Check AlphaFold domain size (≥ min_domain_size aa)
        if domainlen < min_domain_size and name not in include_list:
            print(f" >> domain too small >> SKIP") if Switch('debug') else None
            Log(f"skipped ptm site because its AlphaFold alphaseq structured domain is <{min_domain_size} aa for acc|site (skipped)", f"{acc}|{site}")
            continue
        if name in include_list:
            print(f" >> domain size {domainlen} (exception for include_list) (kept)")



        # No natural variants (other than singletons) in gnomAD (i.e. site is functionally important)
        # (snps_gnomad) = Panda(f"SELECT * FROM snps_gnomad WHERE acc='{acc}' AND site='{site}'")
        # Allow some variants:
        (snps_gnomad) = Panda(f"SELECT * FROM snps_gnomad WHERE acc='{acc}' AND site='{site}' AND {max_variant_metric}>{max_variant}")
        if len(snps_gnomad) > 0 and name not in include_list:
            print(f" >> gnomad variants {max_variant_metric}>{max_variant} >> SKIP") if Switch('debug') else None
            Log(f"skipped ptm site because there are snps_gnomad variants with {max_variant_metric}>{max_variant} in this protein for acc|site (skipped)", f"{acc}|{site}")
            continue

        # Conserved residue (i.e. site is functionally important)
        # query = Query(f"SELECT DISTINCT r.rate FROM evorate_lichtarge_einsi_tree r, uniens ue WHERE ue.acc='{acc}' AND ue.ensp=r.ensp AND r.site='{site}'")
        evorate = FetchOne(Query(f"SELECT MAX(r.rate) FROM evorate_rate4site_einsi_tree_1para r, uniens ue WHERE ue.acc='{acc}' AND ue.ensp=r.ensp AND r.site='{site}'"))
        # Tolerate missing evorate data
        if evorate is None:
            print(f" >> no evorate >> SKIP") if Switch('debug') else None
            Log("skipped ptm site because there is no evorate_lichtarge_einsi_tree data for acc|site (skipped)", f"{acc}|{site}")
            continue
            # evorate = -1
        # If evorate data exists, check if evorate is low enough (= if conservation is good enough)
        if evorate >= 0.3 and name not in include_list:
            # print(f" >> evorate_lichtarge_einsi_tree is too high at >= 3 >> SKIP") if Switch('debug') else None
            # Log("skipped ptm site because evorate_lichtarge_einsi_tree was too high (too fast-evolving) for acc|site (skipped)", f"{acc}|{site}")
            print(f" >> evorate_rate4site_einsi_tree_1para is too high at >= 0.3 >> SKIP") if Switch('debug') else None
            Log("skipped ptm site because evorate_rate4site_einsi_tree_1para was too high (too fast-evolving) for acc|site (skipped)", f"{acc}|{site}")
            continue
        
        # Get upstream proline isomerization state (-1)
        npro = ' '
        query = Query(f"SELECT aa, iso FROM alphasa WHERE acc='{afdb_acc}' AND site='{site - 1}'")
        if Numrows(query) == 1:
            npro_aa, npro_iso = FetchOne(query)
            if npro_aa == 'P':
                npro = npro_iso
                # Spell out isomerisation state
                if npro == 'c':
                    npro = 'cis'
                elif npro == 't':
                    npro = 'trans'
        # Get downstream proline isomerization state (+1)
        cpro = ' '
        query = Query(f"SELECT aa, iso FROM alphasa WHERE acc='{afdb_acc}' AND site='{site + 1}'")
        if Numrows(query) == 1:
            cpro_aa, cpro_iso = FetchOne(query)
            if cpro_aa == 'P':
                cpro = cpro_iso
                # Spell out isomerisation state
                if cpro == 'c':
                    cpro = 'cis'
                elif cpro == 't':
                    cpro = 'trans'

        # Successfully passed all filters
        print(f"\n\n\n >> {name} >> {symbol} >> {fullname} >> {acc} >> {site} >> {ptm} >> {evorate}")
        print(f"   >> OK (successfully passed all filters!)")
        Log("successfully passed all filters for ptm site for acc|site (kept)", f"{acc}|{site}")

        print("\n\n")
        print(f"alphasa:\n{alphasa}")
        print(f"unipdb:\n{unipdb}")
        print(f"snps_clinvar:\n{snps_clinvar}")
        print(f"evorate:\n{evorate}")
        print("\n\n")

        # Process for printing

        # ClinVar
        if len(snps_clinvar) > 0:
            clinsig = max(snps_clinvar['clinsig'])
            clingen = max(snps_clinvar['genetictest'])
            clinsites = ', '.join(nsort(set(str(i) for i in snps_clinvar['site'].to_list())))
            clinorigin = '|'.join(nsort(set(snps_clinvar['originsimple'])))
            clinphen = '|'.join(nsort(set(snps_clinvar['phenotype'].dropna())))
        else:
            clinsig = 0
            clingen = 0
            clinsites = ''
            clinorigin = ''
            clinphen = ''

        # PDB
        # pdbs = ', '.join(nsort(set(unipdb['pdb'].to_list())))
        pdbs = ', '.join(unipdb['pdb'].to_list())

        # Write to output file (tab-separated)
        out.write(f"{name}\t{symbol}\t{fullname}\t{acc}\thttps://www.uniprot.org/uniprotkb/{acc}/entry#structure\thttps://alphafold.ebi.ac.uk/entry/{acc}\thttps://alphasync.stjude.org/display/{acc}\t{site}\t{ptm}\t{sources}\t{source_subset}\t{npro}\t{cpro}\t{surfdist}\t{alphasa['relasa'][0]}\t{alphasa['surf'][0]}\t{alphasa['surf10'][0]}\t{alphasa['dis'][0]}\t{alphasa['dis10'][0]}\t{alphasa['plddt'][0]}\t{alphasa['plddt10'][0]}\t{alphasa['sec'][0]}\t{evorate}\t{pdbs}\t{clinsig}\t{clingen}\t{clinsites}\t{clinorigin}\t{clinphen}\n")

Show(lim=0)
Show("successfully passed all filters for ptm site for acc|site (kept)")

print(f"\nWrote to '{outfile}'\n")

Run("cat", f"cat {outfile}")

print("\nDone!")
