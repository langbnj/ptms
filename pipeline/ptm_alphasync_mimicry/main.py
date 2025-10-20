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

# min_mimicspecies = 3
# min_idmimicfrac = 0.9

# Args(0, "", "")
(species, mafftmode, ptm, mimics, mimic_contacts, evorate_type, evorate_percentile, shared_clade, minseqs, min_mimicspecies, min_idmimicfrac, min_orthomimics, min_orthomimicfrac) = Args(13, "[species] [MAFFT mode: einsi_tree_1para, linsi_tree/ginsi_tree or einsi/linsi/ginsi] [PTM of interest: e.g. S-p for serine phosphorylation] [expected alternative residues, e.g. DE for glutamate or aspartate] [expected residues contacting the mimic, e.g. KR for lysine or asparagine] [evolutionary rate type: capra0/lichtarge/rate4site] [evolutionary rate percentile threshold] [shared clade with identical residue] [minimum number of sequences in an aln] [minimum number of mimic residues] [minimum fraction of identical or mimic residues]", "human einsi_tree_1para S-p DE KR rate4site 4 primates 10 1 0.9 1 0.9")

# Define tables
comparafasta = f"comparafasta_{mafftmode}"
evorate_table = f"evorate_{evorate_type}_{mafftmode}"

# Output file
outfile = f"output-{species}-{mafftmode}-{ptm}-{mimics}-{mimic_contacts}-{evorate_type}-{evorate_percentile}-{shared_clade}-{minseqs}-{min_mimicspecies}-{min_idmimicfrac}-{min_orthomimics}-{min_orthomimicfrac}.txt"

# Get expected aa from ptm (e.g. "S" from "S-p")
expected_aa = ptm.split("-")[0]

# Split e.g. DE into a list (['D', 'E'])
mimics = list(mimics)

# Split e.g. KR into a list (['K', 'R'])
mimic_contacts = list(mimic_contacts)

# Convert to lower case
shared_clade = shared_clade.lower()

# Get number of species in Ensembl Compara
(comparaspecies_count) = FetchOne(Query(f"SELECT COUNT(DISTINCT(species)) FROM comparaspecies"))

# Get full species name from unispec name (e.g. "homo_sapiens" from "HUMAN")
(fullspecies) = FetchOne(Query(f"SELECT species FROM comparaspecies WHERE unispec='{species}'"))



# Possible 1-letter codes for contact types:
# a    AromaticContacts
# c    CarbonylContacts
# C    CovalentContacts
# H    HBondContacts
# h    WeakHBondContacts
# p    HydrophobicContacts
# i    IonicContacts
# m    MetalContacts
# P    PolarHBondContacts
# p    WeakPolarHBondContacts
# v    VanDerWaalsContacts




# # Get core ENSPs (that are present in all mafftmode/tree combinations)
# print("Getting core ENSPs (that are present in all mafftmode/tree combinations):")
# core_ensps = set()
# for mafftmode in ("einsi", "ginsi", "linsi"):
#     for tree in ("", "_tree"):
#         # Construct table name
#         table = f"comparafasta_{mafftmode}{tree}"
#         these_ensps = FetchSet(Query(f"SELECT DISTINCT ensp FROM {table} WHERE species='{fullspecies}'"))
#         if len(core_ensps) > 0:
#             core_ensps = core_ensps.intersection(these_ensps)
#         else:
#             core_ensps = these_ensps
#         print(f" >> {table} >> {len(these_ensps)}")
#         print(f"   >> total {len(core_ensps)}")






# Function definitions

# Get most terminal clade that unites a set of species
def get_clade(species_set, reference_species=fullspecies, all=False):
    """Get clades shared by all members of a set of species.

    Args:
        species (set): Set of species in Ensembl "fullspecies" format (e.g. "pan_troglodytes").
        reference_species (string): Reference species in "fullspecies" format (e.g. "homo_sapiens").
        all (boolean): True returns all shared clades, False returns only the most terminal clade ("lowest", most "detailed").

    Returns:
        A list of shared clades (or, with all=False (default), just the most terminal one).
    """

    # # If no reference species given (reference_species=None), use first species in the species_set
    # Not a good approach
    # if reference_species is None:
    #     if len(species_set) > 0:
    #         reference_species = species_set[0]

    # Get parents of the reference species
    # Retrieve taxon ID for this species from table 'comparaspecies'
    tax = FetchOne(Query(f"SELECT tax FROM comparaspecies WHERE species='{reference_species}'"))
    # Retrieve parents from table 'tax' (NCBI Taxonomy)
    reference_parents = FetchOne(Query(f"SELECT DISTINCT parents FROM tax WHERE tax='{tax}'"))
    # Convert to list, starting from the most terminal taxon
    reference_parents = reference_parents.split('|')
    # Initialize occurrences count dictionary
    reference_parents_occurrences = {parent: 0 for parent in reference_parents}

    for species in species_set:
        # Retrieve taxon ID for this species from table 'comparaspecies'
        tax = FetchOne(Query(f"SELECT tax FROM comparaspecies WHERE species='{species}'"))

        # Retrieve parents from table 'tax' (NCBI Taxonomy)
        parents = FetchOne(Query(f"SELECT DISTINCT parents FROM tax WHERE tax='{tax}'"))
        # Convert to list, starting from the most terminal taxon
        parents = parents.split('|')
        
        for parent in parents:
            # Record how many times each parent occured
            if parent in reference_parents_occurrences:
                reference_parents_occurrences[parent] += 1
    
    # Find most terminal ("lowest", listed first) clade taxon ID that is shared by all species in species_set
    cladetax = []
    for parent in reference_parents_occurrences:
        if reference_parents_occurrences[parent] == max(reference_parents_occurrences.values()):
            cladetax.append(parent)
            # By breaking here, we'll get only the most terminal taxon shared among all species
            if all == False:
                break

    # clade should never be empty at this point (there should always be a shared taxon ID "1", at the very least)
    if len(cladetax) == 0:
        Die(f"Error: clade remained empty for species_set '{', '.join(species_set)}' and reference_species '{reference_species}'")

    # Get scientific name (e.g. cladetax=9604 to clade="Hominidae")
    clades = []
    for tax in cladetax:
        clades.append(FetchOne(Query(f"SELECT DISTINCT name FROM tax WHERE tax='{tax}' AND type='scientific name'")))

    # Return
    if len(clades) == 1:
        # Return string
        return clades[0]
    else:
        # Return list
        return clades

# Set up alias
get_clades = functools.partial(get_clade, all=True)

# # Test get_most_terminal_clade function
# print(f"species: {species}")
# Time(1)
# species_set = {'homo_sapiens', 'pongo_abelii', 'pan_paniscus'}
# # species_set = {'homo_sapiens', 'pongo_abelii', 'pan_paniscus', 'saccharomyces_cerevisiae'}
# # species_set = {'homo_sapiens', 'pongo_abelii', 'pan_paniscus', 'caenorhabditis_elegans'}
# # species_set = {'homo_sapiens', 'pongo_abelii', 'pan_paniscus', 'drosophila_melanogaster'}
# # species_set = {'homo_sapiens', 'gorilla_gorilla', 'pongo_abelii', 'pan_paniscus', 'drosophila_melanogaster', 'caenorhabditis_elegans', 'saccharomyces_cerevisiae'}
# # species_set = {'homo_sapiens', 'gorilla_gorilla', 'pongo_abelii', 'pan_paniscus', 'mus_musculus', 'canis_lupus_familiaris'}
# print(species_set)
# # clade = get_clade(species_set)
# # clade = get_clade(species_set, all=True)
# clades = get_clades(species_set)
# print(clades)
# print(clades[0])
# Time(1)
# print(f"species: {species}")



# Utility function to get nth amino acid character in a sequence that contains gaps (i.e. to get the alignment position corresponding to the nth amino acid in the reference species)
def site_to_alnsite(seq, site, pattern=r"[ACDEFGHIKLMNPQRSTVWY]"):
    """Find the index of the nth occurrence of a regular expression pattern in a string.

    Args:
        seq (str): The input string to search within.
        site (int): The desired occurrence number (1-based).
        pattern (str, optional): The regular expression pattern to match. Defaults to amino acid characters [ACDEFGHIKLMNPQRSTVWY].

    Returns:
        int: The index of the nth occurrence of the pattern in the string (1-based).
             If the pattern is not found or site is invalid (0 or negative), -1 is returned.
    """
    occ = [m.start() + 1 for m in re.finditer(pattern, seq, re.IGNORECASE)]
    return occ[site-1] if len(occ) >= site and site > 0 else Die("Error: site_to_alnsite(): site out of range")



# Utility function to get the occurrence number of a regular expression pattern in a string at a specific index (i.e. to get the species position corresponding to the nth amino acid in the alignment)
def alnsite_to_site(seq, alnsite, pattern=r"[ACDEFGHIKLMNPQRSTVWY]"):
    """Find the occurrence number of a regular expression pattern in a string at a specific index.

    Args:
        seq (str): The input string to search within.
        alnsite (int): The index of the pattern in the string (1-based).
        pattern (str, optional): The regular expression pattern to match. Defaults to amino acid characters [ACDEFGHIKLMNPQRSTVWY].

    Returns:
        int: The occurrence number of the pattern in the string (1-based).
             If the pattern is not found or site is invalid (0 or negative), -1 is returned.
    """
    occ = [m.start() + 1 for m in re.finditer(pattern, seq, re.IGNORECASE)]
    # return occ.index(alnsite) + 1 if alnsite in occ else Die("Error: alnsite_to_site(): alnsite out of range")
    return occ.index(alnsite) + 1 if alnsite in occ else None



# Get alignment column at a specific position
def get_alncol(ensp, site, mafftmode=mafftmode, reference_species=fullspecies):

    # Construct table name
    table = f"comparafasta_{mafftmode}"

    # Get reference sequence (human) and alignment (aln)
    (aln, seq) = FetchOne(Query(f"SELECT aln, seq FROM {table} WHERE ensp='{ensp}' AND species='{reference_species}'"))
    
    # Convert site in reference sequence to alignment position
    alnsite = site_to_alnsite(seq, site)
    
    # Check column for species with an identical amino acid to the reference species (human)
    # Get reference amino acid (human)
    aa = seq[alnsite - 1]
    if aa != expected_aa:
        Die(f"Error: site is not the expected amino acid ('{aa}' instead of '{expected_aa}')")

    # Skip gaps
    if aa == "-":
        Die("Error: site is a gap")
    else:
        # Order these by comparaspecies table and return as a concatenated string (always length 200, i.e. all species in comparaspecies)
        # e.g. SELECT GROUP_CONCAT(IFNULL(SUBSTRING(seq, 983, 1), ' ') ORDER BY s.id SEPARATOR '') FROM comparaspecies s LEFT OUTER JOIN comparafasta_einsi_tree_1para f on f.species=s.species AND f.aln=98;
        # >> "VVV V VVVVVVVVVAAAA   VAVVVVV GGVVV EVVVVV VVV  VV      GVV   VV   V -VVG  G GGG   G  -                                                                                                                 "
        # Spaces are NULL values (species that don't occur in the aln)
        alncol = FetchOne(Query(f"SELECT GROUP_CONCAT(IFNULL(SUBSTRING(seq, {alnsite}, 1), ' ') ORDER BY s.id SEPARATOR '') FROM comparaspecies s LEFT OUTER JOIN {table} f on f.species=s.species AND f.aln={aln}"))
        if len(alncol) != comparaspecies_count:
            Die(f"Error: Expected length of alncol to be equal to the number of species in Ensembl Compara, but it is {len(alncol)} instead:\n\n['{alncol}']\n\n")
        return alncol




# Get species with identical residues at a specific position
def get_species_with_identical_residues(ensp, site, mafftmode=mafftmode, reference_species=fullspecies):
    """Get species with identical residues at a specific position

    Args:
        ensp (str): ENSP ID of the reference protein (i.e. human)
        site (int): Position within the reference protein (not counting gaps, i.e. this isn't identical to the alignment position)
        mafftmode (str, optional): MAFFT mode and whether a species tree is used. Defaults to "einsi_tree_1para", which performed best in benchmarks.
        reference_species (str, optional): Reference species to use, e.g. "homo_sapiens". Defaults to argument value.

    Returns:
        list: List of species that have identical residues at the specified position, ordered as in table 'comparaspecies'.
    """

    # Construct table name
    table = f"comparafasta_{mafftmode}"

    # Get reference sequence (human) and alignment (aln)
    (aln, seq) = FetchOne(Query(f"SELECT aln, seq FROM {table} WHERE ensp='{ensp}' AND species='{reference_species}'"))
    
    # Convert site in reference sequence to alignment position
    alnsite = site_to_alnsite(seq, site)
    
    # Check column for species with an identical amino acid to the reference species (human)
    # Get reference amino acid (human)
    aa = seq[alnsite - 1]
    if aa != expected_aa:
        Die(f"Error: Site is not the expected amino acid ('{aa}' instead of '{expected_aa}')")

    # Skip gaps
    if aa == "-":
        Die("Error: Site is a gap")

    # # No point in ordering these since they're not ordered in the comparafasta tables anyway (human comes first, then alphabetic, but they're not ordered by phylogeny)
    # return FetchSet(Query(f"SELECT DISTINCT species FROM {table} WHERE aln='{aln}' AND species!='{reference_species}' AND SUBSTRING(seq, {alnsite}, 1)='{aa}'"))
    # Order species by table 'comparaspecies'
    return FetchList(Query(f"SELECT s.species FROM {table} f, comparaspecies s WHERE f.species=s.species AND f.aln='{aln}' AND f.species!='{reference_species}' AND SUBSTRING(f.seq, {alnsite}, 1)='{aa}' GROUP BY s.id ORDER BY s.id"))



# Get species with mimic residues at a specific position
def get_species_with_mimics(ensp, site, mimics=mimics, mafftmode=mafftmode, reference_species=fullspecies):
    """Get species with mimic residues at a specific position

    Args:
        ensp (str): ENSP ID of the reference protein (i.e. human)
        site (int): Position within the reference protein (not counting gaps, i.e. this isn't identical to the alignment position)
        mimics (list, optional): List of mimic residues, e.g. ['D', 'E'] for negatively charged residues. Defaults to argument value.
        mafftmode (str, optional): MAFFT mode and whether a species tree is used. Defaults to "einsi_tree_1para", which performed best in benchmarks.
        reference_species (str, optional): Reference species to use, e.g. "homo_sapiens". Defaults to argument value.

    Returns:
        list: List of species that have mimic residues at the specified position, ordered as in table 'comparaspecies'.
    """
    
    # Construct table name
    table = f"comparafasta_{mafftmode}"

    # Get reference sequence (human) and alignment (aln)
    (aln, seq) = FetchOne(Query(f"SELECT aln, seq FROM {table} WHERE ensp='{ensp}' AND species='{reference_species}'"))

    # Convert site in reference sequence to alignment position
    alnsite = site_to_alnsite(seq, site)
    
    # Check column for species with an identical amino acid to the reference species (human)
    # Get reference amino acid (human)
    aa = seq[alnsite - 1]
    if aa != expected_aa:
        Die(f"Error: Site is not the expected amino acid ('{aa}' instead of '{expected_aa}')")

    # Skip gaps
    if aa == "-":
        Die("Error: Site is a gap")

    # # No point in ordering these since they're not ordered in the comparafasta tables anyway (human comes first, then alphabetic, but they're not ordered by phylogeny)
    # return FetchSet(Query(f"SELECT DISTINCT species FROM {table} WHERE aln='{aln}' AND species!='{reference_species}' AND SUBSTRING(seq, {alnsite}, 1) IN ('" + "', '".join(mimics) + "')"))
    # Order species by table 'comparaspecies'
    return FetchList(Query(f"SELECT s.species FROM {table} f, comparaspecies s WHERE f.species=s.species AND f.aln='{aln}' AND f.species!='{reference_species}' AND SUBSTRING(f.seq, {alnsite}, 1) IN ('" + "', '".join(mimics) + "') GROUP BY s.id ORDER BY s.id"))



# Get species with mimic residues at a specific position
def get_species_with_identical_residues_or_mimics(ensp, site, mimics=mimics, mafftmode=mafftmode, reference_species=fullspecies):
    """Get species with identical or mimic residues at a specific position

    Args:
        ensp (str): ENSP ID of the reference protein (i.e. human)
        site (int): Position within the reference protein (not counting gaps, i.e. this isn't identical to the alignment position)
        mimics (list, optional): List of mimic residues, e.g. ['D', 'E'] for negatively charged residues. Defaults to argument value.
        mafftmode (str, optional): MAFFT mode and whether a species tree is used. Defaults to "einsi_tree_1para", which performed best in benchmarks.
        reference_species (str, optional): Reference species to use, e.g. "homo_sapiens". Defaults to argument value.

    Returns:
        list: List of species that have identical or mimic residues at the specified position, ordered as in table 'comparaspecies'.
    """
    
    # Construct table name
    table = f"comparafasta_{mafftmode}"

    # Get reference sequence (human) and alignment (aln)
    (aln, seq) = FetchOne(Query(f"SELECT aln, seq FROM {table} WHERE ensp='{ensp}' AND species='{reference_species}'"))

    # Convert site in reference sequence to alignment position
    alnsite = site_to_alnsite(seq, site)
    
    # Check column for species with an identical amino acid to the reference species (human)
    # Get reference amino acid (human)
    aa = seq[alnsite - 1]
    if aa != expected_aa:
        Die(f"Error: Site is not the expected amino acid ('{aa}' instead of '{expected_aa}')")

    # Skip gaps
    if aa == "-":
        Die("Error: Site is a gap")

    # # No point in ordering these since they're not ordered in the comparafasta tables anyway (human comes first, then alphabetic, but they're not ordered by phylogeny)
    # return FetchSet(Query(f"SELECT DISTINCT species FROM {table} WHERE aln='{aln}' AND species!='{reference_species}' AND SUBSTRING(seq, {alnsite}, 1) IN ('" + "', '".join(mimics) + "')"))
    # Order species by table 'comparaspecies'
    # Add aa to mimics (to also get identical residues)
    return FetchList(Query(f"SELECT s.species FROM {table} f, comparaspecies s WHERE f.species=s.species AND f.aln='{aln}' AND f.species!='{reference_species}' AND SUBSTRING(f.seq, {alnsite}, 1) IN ('{aa}', '" + "', '".join(mimics) + "') GROUP BY s.id ORDER BY s.id"))



# Convert contact type to 1-letter code
def contact_type_to_one_letter(type):
    """Convert contact type to 1-letter code

    Args:
        type (str): Contact type (e.g. "AromaticContacts")

    Returns:
        code: 1-letter code (e.g. "a")
    """
    # a    AromaticContacts
    # c    CarbonylContacts
    # C    CovalentContacts
    # H    HBondContacts
    # h    WeakHBondContacts
    # p    HydrophobicContacts
    # i    IonicContacts
    # m    MetalContacts
    # P    PolarHBondContacts
    # p    WeakPolarHBondContacts
    # v    VanDerWaalsContacts
    if type == "AromaticContacts":
        return "a"
    elif type == "CarbonylContacts":
        return "c"
    elif type == "CovalentContacts":
        return "C"
    elif type == "HBondContacts":
        return "H"
    elif type == "WeakHBondContacts":
        return "h"
    elif type == "HydrophobicContacts":
        return "p"
    elif type == "IonicContacts":
        return "i"
    elif type == "MetalContacts":
        return "m"
    elif type == "PolarHBondContacts":
        return "P"
    elif type == "WeakPolarHBondContacts":
        return "p"
    elif type == "VanDerWaalsContacts":
        return "v"
    else:
        Die(f"Error: contact_type_to_one_letter(): Unknown contact type '{type}'")



# Convert a list of contact types to 1-letter codes
def contact_types_to_one_letter(types):
    """Convert a list of contact types to 1-letter codes

    Args:
        types (list): List of contact types (e.g. ["AromaticContacts", "HBondContacts"])

    Returns:
        codes: List of 1-letter codes (e.g. ["a", "H"])
    """
    codes = []
    for type in sorted(types):
        codes.append(contact_type_to_one_letter(type))
    return codes





# Start

# Calculate evolutionary rate quantiles

# # Get all evolutionary rate values within the set of modified proteins
# Get all evolutionary rate values
# query = Query(f"SELECT rate FROM {evorate_table}")

# Find interesting PTM sites (mimicry/regulated unfolding)

# Get PTM sites
mainquery = Query(f"SELECT DISTINCT acc, name, site FROM unimod WHERE species='{species}' AND ptm='{ptm}' ORDER BY name, acc, site")

with open(outfile, 'w') as out:
    
    # Write header
    out.write(f"acc\tsite\tfullname\taln\talnseqs\taa\tidfrac\tmimicfrac\tidmimicfrac\trelasa\tplddt\tconsites\torthoconfrac\tidclade\tmimicclade\tidmimicclade\torthoconspecies_all_clade\torthocontypes\torthoconfracs\torthomimicfracs\n")

    for acc, name, site in tq(mainquery):

        if Switch('debug'):
            print(f"\n >> acc {acc}") 

        # Get full protein name
        fullname = FetchOne(Query(f"SELECT fullname FROM uniprot WHERE name='{name}' AND species='{species}'"))

        if Switch('debug'):
            print(f"   >> name {name}")
            print(f"   >> fullname {fullname}")
            print(f"   >> site {site}")

        # Need to skip isoform accessions since alphaseq doesn't contain them
        # Skip any accessions that contain a hyphen (i.e. isoforms)
        if "-" in acc:
            Log("reference species isoform acc skipped since it wouldn't be in alphaseq for acc (skipped)", acc)
            continue

        # Skip if uniseq's sequence doesn't match alphaseq's
        # Try to get alphaseq sequence
        query = Query(f"SELECT DISTINCT seq FROM alphaseq WHERE acc='{acc}'")
        # Skip if this acc isn't in alphaseq
        if Numrows(query) == 0:
            Log("no seq in alphaseq for acc (skipped)", acc)
            continue
        alphaseq = FetchOne(query)

        # Get UniProt sequence
        uniseq = FetchOne(Query(f"SELECT DISTINCT seq FROM uniseq WHERE type IN ('UniProt', 'UniIso') AND acc='{acc}'"))
        if alphaseq != uniseq:
            Log("reference species's alphaseq sequence doesn't match uniseq's for acc (skipped)", acc)
            continue
        
        # Get ENSP
        enspquery = Query(f"SELECT f.ensp FROM {comparafasta} f, uniens ue WHERE ue.acc='{acc}' AND ue.ensp=f.ensp AND f.species='{fullspecies}'")
        # Skip if no ENSP for this acc
        if Numrows(enspquery) == 0:
            Log(f"no ensp found in either {comparafasta} or uniens for acc (skipped)", acc)
            continue
        for (ensp,) in enspquery:
            # ensp = FetchOne(query)
            if Switch('debug'):
                print(f"   >> ensp {ensp}")

            # Get alignment (aln)
            (aln, seq) = FetchOne(Query(f"SELECT aln, seq FROM {comparafasta} WHERE ensp='{ensp}' AND species='{fullspecies}'"))
            if Switch('debug'):
                print(f"   >> aln {aln}")
                print(f"   >> seq {seq}")

            # Get alignment position
            alnsite = site_to_alnsite(seq, site)
            if Switch('debug'):
                print(f"   >> alnsite {alnsite}")
        
            # Get amino acid
            aa = seq[alnsite - 1]
            if Switch('debug'):
                print(f"   >> aa {aa}")
                
            # Get alignment column
            alncol = get_alncol(ensp, site)
            if Switch('debug'):
                print(f"   >> alncol '{alncol}'")
            
            # Get total number of species in this aln
            alnseqs = FetchOne(Query(f"SELECT COUNT(*) FROM {comparafasta} WHERE aln='{aln}' AND species!='{fullspecies}'"))
            if Switch('debug'):
                print(f"   >> alnseqs {alnseqs}")

            
            # Identical to reference species
            
            # Get species with identical residues
            idspecies = get_species_with_identical_residues(ensp, site)
            if Switch('debug'):
                print(f"   >> idspecies")
                print(f"     >> {idspecies}")

            # Get fraction of species with identical residues
            idfrac = len(idspecies) / alnseqs
            if Switch('debug'):
                print(f"   >> idfrac {idfrac}")

            # Get shared clades between species with identical residues
            idclades = get_clades(idspecies)
            if Switch('debug'):
                print(f"   >> idclades")
                print(f"     >> {idclades}")
            
            
            # Mimics

            # Get species with mimic residues
            mimicspecies = get_species_with_mimics(ensp, site, mimics)
            if Switch('debug'):
                print(f"   >> mimicspecies")
                print(f"     >> {mimicspecies}")
            # Apply min_mimicspecies threshold
            if len(mimicspecies) < min_mimicspecies:
                Log(f"site skipped because there were fewer than {min_mimicspecies} species with mimic residues in the aln for acc|site (skipped)", f"{acc}|{site}")
                continue

            # Get fraction of species with mimic residues
            mimicfrac = len(mimicspecies) / alnseqs
            if Switch('debug'):
                print(f"   >> mimicfrac {mimicfrac}")

            # Get shared clades between species with mimic residues
            mimicclades = get_clades(mimicspecies)
            if Switch('debug'):
                print(f"   >> mimicclades")
                print(f"     >> {mimicclades}")

            
            # Identical or mimics
            
            # Get species with mimic residues
            idmimicspecies = get_species_with_identical_residues_or_mimics(ensp, site, mimics)
            if Switch('debug'):
                print(f"   >> idmimicspecies")
                print(f"     >> {idmimicspecies}")

            # Get fraction of species with identical or mimic residues
            idmimicfrac = len(idmimicspecies) / alnseqs
            if Switch('debug'):
                print(f"   >> idmimicfrac {idmimicfrac}")
            # Apply idmimicfrac threshold
            if idmimicfrac < min_idmimicfrac:
                Log(f"site skipped because the fraction of species with identical or mimic residues was less than {min_idmimicfrac} in the aln for acc|site (skipped)", f"{acc}|{site}")
                continue

            # Get shared clades between species with identical or mimic residues
            idmimicclades = get_clades(idmimicspecies)
            if Switch('debug'):
                print(f"   >> idmimicclades")
                print(f"     >> {idmimicclades}")

            # # Verify that idmimicfrac is the sum of idfrac and mimicfrac (should always be the case) (i.e. no overlap)
            # if idmimicfrac != idfrac + mimicfrac:
            #     Die(f"Error: idmimicfrac ({idmimicfrac}) != idfrac + mimicfrac ({idfrac + mimicfrac})")
            # Verify that there is no overlap between idspecies and mimicspecies (should always be the case) (i.e. no overlap)
            if len(set(idspecies).intersection(set(mimicspecies))) > 0:
                Die(f"Error: There is overlap between the idspecies and mimicspecies sets ({set(idspecies).intersection(set(mimicspecies))})")
            if set(idspecies).union(set(mimicspecies)) != set(idmimicspecies):
                Die(f"Error: The union of idspecies and mimicspecies ({set(idspecies).union(set(mimicspecies))}) is not equal to idmimicspecies ({set(idmimicspecies)})")



            # Get relative accessible surface area and pLDDT (AlphaFold confidence score)
            (relasa, plddt) = FetchOne(Query(f"SELECT relasa, plddt FROM alphasa WHERE acc='{acc}' AND site='{site}'"))
            if Switch('debug'):
                print(f"   >> relasa {relasa}")
                print(f"   >> plddt {plddt}")

            # Get contacts
            # Add "...1" variables for clarity
            site1 = site
            alnsite1 = alnsite
            aa1 = aa
            # Contact and ortholog counters
            # orthocontacts = 0
            # orthologs = 0
            orthoconfracs = []
            orthomimiccounts = []
            orthomimicfracs = []
            orthoconspecies_all = set()
            consites = []
            # First, get human contacts
            for (site2,) in Query(f"SELECT site1 FROM alphacon WHERE acc='{acc}' AND site2='{site}' UNION SELECT site2 FROM alphacon WHERE acc='{acc}' AND site1='{site}'"):

                # Count orthologs with contacts and total orthologs for this contact (site2)
                orthoconspecies = []
                orthocontacts = 0
                orthologs = 0

                # Get alignment position of contacted residue
                alnsite2 = site_to_alnsite(seq, site2)

                # Get amino acid at contacted residue
                aa2 = seq[alnsite2 - 1]

                if Switch('debug'):
                    print(f"     >> Contact {site1}-{site2}")
                    print(f"       >> alnsites {alnsite1}-{alnsite2}")
                    print(f"       >> aas {aa1}-{aa2}")
                    

                # Get fraction of orthologs that also have a contact between the orthologous residues
                # for (orthoensp, orthoseq) in Query(f"SELECT ensp, seq FROM {comparafasta} WHERE aln='{aln}' AND species!='{fullspecies}'"):
                # Order by comparaspecies, and include species that don't have an ortholog
                # f"SELECT GROUP_CONCAT(IFNULL(SUBSTRING(seq, {alnsite}, 1), ' ') ORDER BY s.id SEPARATOR '') FROM comparaspecies s LEFT OUTER JOIN {table} f on f.species=s.species AND f.aln={aln}"
                orthocol1 = ""
                orthocol2 = ""
                orthoconcols = []
                orthocontypes = set()
                for (orthoensp, orthoseq, orthospecies) in Query(f"SELECT f.ensp, f.seq, f.species FROM comparaspecies s LEFT OUTER JOIN {comparafasta} f ON f.species=s.species AND aln='{aln}' AND f.species!='{fullspecies}' ORDER BY s.id"):
                    
                    # For species where there is no ortholog: append space to orthocol
                    if orthoensp is None:
                        orthocol1 += " "
                        orthocol2 += " "
                        continue
                    
                    # Get orthologous contact positions within the ortholog's sequence
                    orthosite1 = alnsite_to_site(orthoseq, alnsite1)
                    orthosite2 = alnsite_to_site(orthoseq, alnsite2)

                    # Skip if one or both of these sites are at gaps or outside the ortholog's sequence (since it obviously doesn't have orthologous residues then)
                    if orthosite1 is None or orthosite2 is None:
                        continue

                    # Strip gaps from orthoseq for alphaseq sequence query
                    orthoseq_nogaps = orthoseq.replace("-", "")

                    # # Get accession for the ortholog via uniens
                    # query = Query(f"SELECT acc FROM uniens WHERE ensp='{orthoensp}'")
                    # # Get accession for the ortholog via sequence
                    # query = Query(f"SELECT acc FROM alphaseq WHERE seq='{orthoseq_nogaps}'")
                    # Get accession for the ortholog via sequence, and order by descending pLDDT score at the PTM site and its contact (get the top structure prediction for this residue)
                    # # Sort by the average pLDDT score across the PTM site and the contacted residue
                    # query = Query(f"SELECT s.acc FROM alphaseq s, alphasa a1, alphasa a2 WHERE a1.site='{orthosite1}' AND a2.site='{orthosite2}' AND a1.acc=s.acc AND a2.acc=s.acc AND s.seq='{orthoseq_nogaps}' ORDER BY ((a1.plddt + a2.plddt) / 2) DESC LIMIT 1")
                    # First sort by the PTM site's pLDDT score (descending), then by the contacted residue's (descending), and get only the "best" structure ("orthoacc", best by pLDDT score)
                    query = Query(f"SELECT s.acc FROM alphaseq s, alphasa a1, alphasa a2 WHERE a1.site='{orthosite1}' AND a2.site='{orthosite2}' AND a1.acc=s.acc AND a2.acc=s.acc AND s.seq='{orthoseq_nogaps}' ORDER BY a1.plddt DESC, a2.plddt DESC LIMIT 1")
                    # Skip if no accession for this ortholog
                    if Numrows(query) == 0:
                        Log("no orthoacc found via sequence in alphaseq for orthoensp (skipped)", orthoensp)
                        Log("no orthoacc found via sequence in alphaseq for reference species acc (skipped)", acc)
                        # if Switch('debug'):
                        if Switch('debug2'):
                            print(f"         >> no acc found in alphaseq via sequence for orthoensp '{orthoensp}' (no AlphaFold prediction exists for this sequence) (skipped)")
                        continue
                    Log("orthoacc successfully found via sequence in alphaseq for orthoensp (kept)", orthoensp)
                    Log("orthoacc successfully found via sequence in alphaseq for reference species acc (kept)", acc)
                    orthoacc = FetchOne(query)

                    # Get amino acids at orthologous residue positions (before swapping, so they're equivalent to aa1 and aa2 above)
                    orthoaa1 = orthoseq[alnsite1 - 1]
                    orthoaa2 = orthoseq[alnsite2 - 1]
                    # Verify
                    if orthoaa1 == '-':
                        Die("Error: orthoaa1 is a gap ('-')")
                    if orthoaa2 == '-':
                        Die("Error: orthoaa2 is a gap ('-')")



                    # Swap sites if orthosite1 > orthosite2
                    if orthosite1 > orthosite2:
                        (orthosite1, orthosite2) = (orthosite2, orthosite1)

                    # Check if orthologous contact exists
                    # query = Query(f"SELECT * FROM alphacon WHERE acc='{orthoacc}' AND site1='{orthosite1}' AND site2='{orthosite2}'")
                    query = Query(f"SELECT type FROM alphacon WHERE acc='{orthoacc}' AND site1='{orthosite1}' AND site2='{orthosite2}'")

                    # Count cases where there is a contact between the orthologous residues
                    if Numrows(query) == 0:
                        
                        # No orthologous contact exists: append slash to orthocols
                        # Build alignment columns across orthologs (including ' ' for species without an ortholog, and '/' for species with an ortholog, but without a contact)
                        orthocol1 += '/'
                        orthocol2 += '/'
                        orthoconcols.append(['/'])
                        
                    else:
                        
                        # Get contact type
                        contact_types = nsort(FetchList(query))
                        contypes = contact_types_to_one_letter(contact_types)

                        # # Filter: Require orthoaa2 to be identical to aa2 (i.e. the orthologous contacted residue must be identical to the original contacted residue)
                        # if orthoaa2 == aa2:
                        orthocontacts += 1

                        orthoconspecies.append(orthospecies)

                        # if Switch('debug'):
                        if Switch('debug2'):
                            print(f"         >> Ortholog {orthoensp} ({orthoacc})")
                            print(f"           >> orthoseq {orthoseq}")
                            print(f"           >> orthosites {orthosite1}-{orthosite2}")
                            print(f"           >> orthoaas {orthoaa1}-{orthoaa2}")

                        # Build alignment columns across orthologs (including ' ' for species without an ortholog, and '/' for species with an ortholog, but without a contact)
                        orthocol1 += orthoaa1
                        orthocol2 += orthoaa2

                        # Add contact type to list of contact types for this contact
                        orthoconcols.append(contypes)

                        # Add these contact types to the set
                        orthocontypes.update(contact_types)

                    # Count the total number of orthologs (counting these here since alnseqs might include ENSPs that don't have a corresponding acc)
                    orthologs += 1

                # Get fraction of orthologs that also have a contact between the orthologous residues
                if orthologs > 0:
                    orthoconfrac = orthocontacts / orthologs
                else:
                    orthoconfrac = 0

                # Add to list of orthoconfracs (so we can get the maximum, i.e. the best "orthologous contact conservation fraction" across all contacts the reference species residue makes, below)
                orthoconfracs.append(orthoconfrac)
                # orthoconspecies_sets.append(orthoconspecies)
                # Add species where the current contact is conserved to set of all species with conserved contacts
                orthoconspecies_all.update(set(orthoconspecies))
                # Get clade shared by all species where the present contact is conserved
                orthoconspecies_clade = get_clade(orthoconspecies)

                # Get fraction of contacted residues that are expected mimic_contacts (e.g. KR for S-p)
                orthomimics = 0
                orthototal = 0
                for orthoaa in orthocol2:
                    if orthoaa != ' ' and orthoaa != '/':
                        orthototal += 1
                        if orthoaa in mimic_contacts:
                            orthomimics += 1
                if orthototal != orthocontacts:
                    Die(f"Error: orthototal {orthototal} != orthocontacts {orthocontacts}")
                # Calculate fraction
                orthomimicfrac = 0
                if orthototal > 0:
                    orthomimicfrac = orthomimics / orthototal
                # Add to list of orthomimicfracs (so we can get the maximum, i.e. the best "orthologous contact mimic fraction" across all contacts the reference species residue makes, below)
                orthomimiccounts.append(orthomimics)
                orthomimicfracs.append(orthomimicfrac)

                # Add to list of contacted (human) residues
                consites.append(site2)

                if Switch('debug'):
                    print(f"       >> orthocol1 '{orthocol1}'")
                    print(f"       >> orthocol2 '{orthocol2}'")

                    # Unfinished: was supposed to display 1-letter code for contact types below each residue 
                    # print(f"       >> contypes  '", end="")
                    # # for typenum in range(max(len(l) for l in orthoconcols)):
                    # #     for pos in range(len(orthoconcols)):
                    # #         print(orthoconcols[pos][typenum], end="")
                    # print(f"'")

                    print(f"       >> orthocontypes {nsort(orthocontypes)}")
                    print(f"       >> orthoconfrac {orthoconfrac} ({orthocontacts} / {orthologs})")
                    print(f"       >> orthoconspecies {orthoconspecies}")
                    print(f"       >> orthoconspecies_clade {orthoconspecies_clade}")
                    print(f"       >> orthomimicfrac {orthomimicfrac} ({orthomimics} / {orthototal})")


            # Get clade shared by all species where there is some type of conserved contact
            orthoconspecies_all_clade = get_clade(orthoconspecies_all)
            if Switch('debug'):
                print(f"     >> orthoconspecies_all_clade {orthoconspecies_all_clade}")



            # Filtering
            print("\n\n")

            # Filter by number of species
            if alnseqs < minseqs:
                Log(f"site skipped because there were fewer than {minseqs} species in the aln for acc|site (skipped)", f"{acc}|{site}")
                continue

            # # Filter by number of species with identical residues (at least 1 species must have an identical residue, other than human)
            # if len(idspecies) < 1:
            #     Log(f"site skipped because there were fewer than 1 species with identical residues in the aln for acc|site (skipped)", f"{acc}|{site}")
            #     continue

            # Filter by number of species with mimic residues (at least 3 species must have a mimic residue, other than human)
            if len(mimicspecies) < min_mimicspecies:
                Log(f"site skipped because there were fewer than {min_mimicspecies} species with mimic residues in the aln for acc|site (skipped)", f"{acc}|{site}")
                continue

            # Filter by mimic fraction (require at least one species with a mimic residue)
            if mimicfrac == 0:
                Log(f"site skipped because the fraction of species with mimic residues was 0 in the aln for acc|site (skipped)", f"{acc}|{site}")
                continue
            
            # Filter by idmimic fraction (require all species to have either an identical or a mimic residue) (very stringent)
            if idmimicfrac < min_idmimicfrac:
                Log(f"site skipped because the fraction of species with identical or mimic residues was less than {min_idmimicfrac} in the aln for acc|site (skipped)", f"{acc}|{site}")
                continue
            
            # # Filter by evolutionary rate quantile (use the top quantile)
            # evorate = FetchOne(Query(f"SELECT r.rate FROM {evorate_table} r, uniens ue WHERE ue.acc='{acc}' AND ue.ensp=r.ensp AND r.site='{site}'"))
            # d()

            # Filter by shared species clade
            # Note that this ensures that all species with identical residues share the clade of interest (shared_clade). However, it does not ensure that all members of the clade of interest have identical residues. It only ensures that there are no identical residues outside of the clade.
            # It provides the "boundary" around the presence of identical residues. All identical residues are contained within this clade.
            # Convert all to lower case
            idclades = [clade.lower() for clade in idclades]
            # Skip if desired clade is not shared among these species
            if shared_clade not in idclades:
                Log(f"site skipped because the clade in which identical residues were shared was not contained within '{shared_clade}' for acc|site (skipped)", f"{acc}|{site}")
                continue
            
            # Filter by relative accessible surface area
            # Filter by pLDDT (AlphaFold confidence score)
            # Filter by contacts
            # Filter by contact conservation
            if len(orthoconfracs) == 0 or max(orthoconfracs) == 0:
                Log(f"site skipped because the fraction of species with conserved contacts between orthogonal residues was 0 for all reference species contacts in the aln for acc|site (skipped)", f"{acc}|{site}")
                continue
            # Filter by mimic_contacts (e.g. KR for S-p)
            if len(orthomimiccounts) == 0 or max(orthomimiccounts) < min_orthomimics:
                Log(f"site skipped because the number of contacted amino acids that were expected mimic_contacts was less than {min_orthomimics} in all orthocol2s for acc|site (skipped)", f"{acc}|{site}")
                continue
            # Filter by mimic_contact fraction (e.g. KR for S-p)
            if len(orthomimicfracs) == 0 or max(orthomimicfracs) < min_orthomimicfrac:
                Log(f"site skipped because the number of contacted amino acids that were expected mimic_contacts was less than {min_orthomimicfrac} in all orthocol2s for acc|site (skipped)", f"{acc}|{site}")
                continue
            
            # Filter by PDB structure availability for the residue and contacted residues (via unipdb)
            if (Numrows(Query(f"SELECT * FROM unipdb WHERE acc='{acc}' AND {site} >= start AND {site} <= stop")) == 0):
                Log(f"site skipped because there were no PDB structures for acc|site (skipped)", f"{acc}|{site}")
                continue
            
            # Show adjacent proline isomerization states? Not sure if relevant here. Implemented it in julian_ptm_site_candidates_md_nmr script instead.
                        
            # Filter by contact with identical as well as mimic residues

            # Filter by alignment quality?

            # Add full protein name etc.

            
            # Site passed all filters
            if Switch('debug'):
                print("SUCCESS\n\n")
            Log("success: site passed all filters for acc|site", f"{acc}|{site}")
            # d()

            # Write to output file
            print(f"{acc}\t{site}\t{ptm}\t{aln}\t{alnseqs}\t{aa}\t{idfrac}\t{mimicfrac}\t{idmimicfrac}\t{relasa}\t{plddt}\t{'|'.join([str(x) for x in consites])}\t{max(orthoconfracs)}\t{idclades[0]}\t{mimicclades[0]}\t{idmimicclades[0]}\t{orthoconspecies_all_clade}\t{'|'.join(nsort(contact_types_to_one_letter(orthocontypes)))}\t{'|'.join([str(x) for x in orthoconfracs])}\t{'|'.join([str(x) for x in orthomimicfracs])}", file=out)
            



# Find interesting PTM sites (buried sites, i.e. regulated unfolding)

Show(lim=0)
Show("success: site passed all filters for acc|site")

print("\nDone!")
