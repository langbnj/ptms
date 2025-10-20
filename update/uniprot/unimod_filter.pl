#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$minsites = 1000;   # Minimum sites to show still-unlabelled, but common PTM types

our $usage = "$0 [-hc] [-other] [-debug]\n\n".
" -hc:    Only use experimentally determined sites within a given species (no 'by similarity', 'probable', 'potential', automatic sequence motifs etc.)\n".
# # " -lcng:  Use low-confidence sites for N-glycosylation (i.e. 'by similarity', 'probable' etc.) (since they are so stereotypical, tripeptide)\n".              # Not a good idea - I need traceable PTM evidence
# " -other: Fill unimod ptm field with 'other' for PTMs other than K-ac etc., instead of leaving it empty\n".                                                     # It doesn't make sense to study "other" as one group. Leaving this out.
" -debug: Don't do any actual updates, just simulate\n";
# # "Note: There used to be a -nohistones, but not anymore! It'd be strange to leave those out, and it's very few proteins.";
args(0);


# # Useful query to get an overview of type|description combinations (all PTMs)
# # SELECT type, description, COUNT(DISTINCT name, site) AS c, GROUP_CONCAT(DISTINCT ptm) AS ptms FROM unimod WHERE (evid LIKE '%ECO:0000006%' OR evid LIKE '%ECO:0000269%') GROUP BY type, description ORDER BY c DESC;


# Set PTM aliases
# Useful query: e.g. SELECT aa, type, description, COUNT(DISTINCT name, site) AS sites FROM unimod WHERE description LIKE '%acetyllysine%' AND species='human' AND evid IN (evid LIKE '\%ECO:0007744%' OR evid LIKE '\%ECO:0000269%' OR evid LIKE '\%ECO:0000305%' OR evid LIKE '\%ECO:0000312%') GROUP BY aa, type, description ORDER BY sites DESC, aa, type, description;
%ptm = ();


# Summary list:
# - Lysine acylations: acetyllysine, crotonyllysine, formyllysine, glutaryllysine, lactoylation, malonyllysine, succinyllysine
# - Lysine methylation: methyllysine
# - Lysine ubiquitination
# - Lysine sumoylation
# - Glycosylation: GlcNAc and GalNAc (N-linked and O-linked) on asparagine, serine and threonine
# - Arginine methylation: methylarginine
# - Phosphorylation: phosphoserine, phosphothreonine, phosphotyrosine
# - S-glutathionyl cysteine
# - Hydroxyproline
# - Serine acetylation
# - Methionine sulfoxidation
# - Cysteine palmitoylation
# - Cysteine nitrosylation
# - Lysine neddylation
# - Glutamine to pyrrolidone carboxylic acid
# - Glutamate to gamma-carboxyglutamic acid
# - Arginine citrullination
# - ADP-ribosylation (on S, C, D, E, K and R)
# - Lysine butyrylation: hydroxyisobutyryllysine, hydroxybutyryllysine, butyryllysine
# End summary

# The PTM types where more sites would be most beneficial are:
# - Lysine methylation (5,163 human sites)
# - Lysine succinylation (1,608)
# - Lysine butyrilation (418)
# - Lysine glutarylation (227)
# - Lysine lactoylation (222)
# - Lysine crotonylation (208)
# - Arginine citrullination (57)


# Lysine acylations
$ptm{'K|modified residue|%acetyllysine%'} = 'K-ac';
$ptm{'K|modified residue|%crotonyllysine%'} = 'K-cro';
$ptm{'K|modified residue|%glutaryllysine%'} = 'K-glu';  # negatively charged
$ptm{'K|modified residue|%malonyllysine%'} = 'K-mal';   # negatively charged
$ptm{'K|modified residue|%succinyllysine%'} = 'K-suc';  # negatively charged

# Lysine methylation
$ptm{'K|modified residue|%methyllysine%'} = 'K-me';
$ptm{'K|modified residue|%methylated lysine%'} = 'K-me';

# Lysine ubiquitination and sumoylation
$ptm{'K|cross-link|%Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter%ubiquitin%'} = 'K-ub';
$ptm{'K|cross-link|%Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter%SUMO%'} = 'K-sum';

# Glycosylation: GlcNAc and GalNAc (N-linked and O-linked). O-linked is catalysed by OGT, N-linked is canonical glycosylation, almost always starting with a GlcNAc.
# 13,110 untraceable human N-linked (GlcNAc...) asparagine sites will be lost without -lcng, but there'll still be 3,060 experimental sites, which is still plenty!
# # Distinguishing between GlcNAc and GalNAc
# $ptm{'N|glycosylation site|%N%link%GlcNAc%'} = 'N-glc';
# $ptm{'N|glycosylation site|%N%link%GalNAc%'} = 'N-gal';
# $ptm{'S|glycosylation site|%O%link%GlcNAc%'} = 'S-glc';
# $ptm{'T|glycosylation site|%O%link%GlcNAc%'} = 'T-glc';
# $ptm{'S|glycosylation site|%O%link%GalNAc%'} = 'S-gal';
# $ptm{'T|glycosylation site|%O%link%GalNAc%'} = 'T-gal';
# Not distinguishing between GlcNAc and GalNac (e.g. dbPTM doesn't have any distinction)
$ptm{'N|glycosylation site|%N%link%GlcNAc%'} = 'N-gly';
$ptm{'N|glycosylation site|%N%link%GalNAc%'} = 'N-gly';
$ptm{'S|glycosylation site|%O%link%GlcNAc%'} = 'S-gly';
$ptm{'T|glycosylation site|%O%link%GlcNAc%'} = 'T-gly';
$ptm{'S|glycosylation site|%O%link%GalNAc%'} = 'S-gly';
$ptm{'T|glycosylation site|%O%link%GalNAc%'} = 'T-gly';

# Arginine methylation
$ptm{'R|modified residue|%methylarginine%'} = 'R-me';
$ptm{'R|modified residue|%methylated arginine%'} = 'R-me';

# Phosphorylation (serine, threonine, tyrosine)
$ptm{'S|modified residue|%phosphoserine%'} = 'S-p';
$ptm{'T|modified residue|%phosphothreonine%'} = 'T-p';
$ptm{'Y|modified residue|%phosphotyrosine%'} = 'Y-p';



# Other PTMs with ≥1000 sites:

# N-terminal acetylation
# # N-acetylmethionine (almost always N-terminal on residue 1)
#  Methionine 1 obviously will be highly conserved, so point looking at that.
# $ptm{'P|modified residue|%acetylmethionine%'} = 'M-ac';
# Alanine acetylation (almost always N-terminal on residue 2)
# I think it's interesting to study whether A-ac residues are highly conserved.
$ptm{'A|modified residue|%acetylalanine%'} = 'A-ac';

# S-glutathionyl cysteine
# UniProt also ~100 sites
$ptm{'C|modified residue|%S-glutathionyl cysteine%'} = 'C-glt';

# Hydroxyproline (mostly 4-hydroxyproline, some unspecified, and some 3-hydroxyproline)
# I think this is interesting as, I think, a hypoxia sensor.
$ptm{'P|modified residue|%hydroxyproline%'} = 'P-hyd';

# Serine acetylation is not a rarity! There are actually many sites! This is the Yersinia pestis PTM, where it acetylates host proteins (2006 paper).
$ptm{'S|modified residue|%acetylserine%'} = 'S-ac';

# Methionine acetylation (nearly always at the N -terminus). Might stablilise/destabilise proteins (I remember reading it can destabilise them - somehow related to the N-end rule, perhaps?).
$ptm{'M|modified residue|%acetylmethionine%'} = 'M-ac';



# Others (from dbPTM)

# Methionine sulfoxidation: interesting since it happens on internal methionines, not terminal ones. >1000 sites in dbPTM.
# "Oxidation" is also a dbPTM term for methionine sulfoxidation (PMID 23648414). Pretty ridiculous that dbPTM has it separately.
$ptm{'M|modified residue|%sulfoxid%'} = 'M-ox';

# Lipid: cysteine palmitoylation (internal, not terminal, hence interesting)
$ptm{'C|lipid moiety-binding region|%S-palmitoyl cysteine%'} = 'C-pal';

# Cysteine nitrosylation
$ptm{'C|modified residue|%nitroso%'} = 'C-nit';

# Amidation: mostly terminal or on small peptides, not interesting (skip)

# Lysine neddylation
$ptm{'K|cross-link|%Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter%NEDD%'} = 'K-ned';

# Glutamine to Pyrrolidone carboxylic acid
$ptm{'Q|modified residue|%pyrrolido%'} = 'Q-pyr';

# Glutamate to Gamma-carboxyglutamic acid
$ptm{'E|modified residue|%carboxyglutam%'} = 'E-car';

# Lysine lactoylation (only one PMID in dbPTM and in UniProt: 31645732)
$ptm{'K|modified residue|%lactoyllysine%'} = 'K-lac';

# Lysine lactylation (only 2 PMIDs in dbPTM: 31645732 (same as for lactoylation), 33193272). Both PMIDs actually call it "lactylation" in the text. However, UniProt has it as "lactoyllysine", as does PMID 31767537. I think these latter guys were trying to correct it from lactylation to lactoylation. I should go with lactoylation.
# It's pretty ludicrous that dbPTM has both spellings as separate files, but oh well.
# Not necessary since UniProt never has "lactylation"
# $ptm{'K|modified residue|%lacty%'} = 'K-lac';

# Myristoylation: N-terminal on glycine, not interesting (skip)

# Lysine formylation
$ptm{'K|modified residue|%formyllysine%'} = 'K-for';

# Arginine citrullination (replaces the residue, irreversible)
$ptm{'R|modified residue|%citrulline%'} = 'R-cit';

# ADP-ribosylation (occurs on lots of different AAs) (just a few hundred sites though)
$ptm{'S|modified residue|%adp%ribosyl%'} = 'S-adp';
$ptm{'C|modified residue|%adp%ribosyl%'} = 'C-adp';
$ptm{'D|modified residue|%adp%ribosyl%'} = 'D-adp';
$ptm{'E|modified residue|%adp%ribosyl%'} = 'E-adp';
$ptm{'K|modified residue|%adp%ribosyl%'} = 'K-adp';
$ptm{'R|modified residue|%adp%ribosyl%'} = 'R-adp';

# Lysine butyrylation (includes hydroxyisobutyryllysine, hydroxybutyryllysine, and butyryllysine)
$ptm{'K|modified residue|%buty%lysine%'} = 'K-but';




# start

if (!switch('debug'))
{
    starttime();
    Query("UPDATE unimod SET ptm=NULL");
    # Query("UPDATE unimod SET ptm=''");		# For historical reasons (it'd break many scripts), I'll keep using '' instead of NULL here
    state("Cleared 'ptm' field in table 'unimod'");
    stoptime();
}



starttime();
$evid = '';
if (switch('hc'))
{
	# Experimental only

    # New (2022_04)
    # https://www.uniprot.org/help/evidences
    # 
    # Evidence codes in UniProt 2022_04:
    # 
    # All species
    # 144933 (15261 human)   ECO:0000250     X SKIP "By similarity": "We use the ECO code ECO:0000250 for manually curated information which has been propagated from a related experimentally characterized protein. The accession number of the experimental source is provided (except for some historic UniProtKB/Swiss-Prot annotations, see Why don't all UniProtKB/Swiss-Prot annotations have evidence? )."
    # 114979 (15102 human)   [NULL]          X SKIP In human, this only happens for N-GlcNAc and O-GalNAc. In other species, it mostly seems to be N-glycosylation and C-palmitoylation, and T-GalNAc, but also N6-(pyridoxal phosphate)lysine (with no PMID and untraceable). Sites with this ECO don't have PMIDs (very bad), nor are they "By similarity", "Probable" or "Potential" (which would be good): they could not be automatically evidence-annotated by looking for "By similarity, Probable and Potential qualifiers that marked information that was not experimental as well as the Cited for section of literature citations which indicates the type of information retrieved from a given reference": "The annotations which are missing evidence were created before we started to manually curate information with evidence attribution in UniProtKB/Swiss-Prot. The manual attribution of evidence to these existing annotations was not possible due to the huge amount of existing data. Therefore we wrote a program to automatically add the evidence based on how provenance information was previously curated in Swiss-Prot including the By similarity, Probable and Potential qualifiers that marked information that was not experimental as well as the Cited for section of literature citations which indicates the type of information retrieved from a given reference. This allowed the attribution of evidence to a substantial fraction of annotations in Swiss-Prot but not to all annotations. The remaining annotations without evidence must be updated and retrofitted manually, which is an ongoing process. (https://www.uniprot.org/help/evidence_in_swissprot)"
    # 71864  (37909 human)   ECO:0007744     KEEP Experimental large-scale mass spec or PDB structure evidence. "published large-scale proteomics data and a subset of 3D structural data for which UniProt makes use of computational procedures to generate the data": "Combinatorial evidence. We use the ECO code ECO:0007744 in manual assertions for information inferred from a combination of experimental and computational evidence. It is currently used in UniProtKB for published large-scale proteomics data and a subset of 3D structural data for which UniProt makes use of computational procedures to generate the data. If the experimental evidence is derived from a paper, the PubMed identifier of the paper is provided. For experimental data derived from the Protein Data Bank, the PDB code is provided."
    # 48875    (532 human)   ECO:0000255     X SKIP Sequence motifs from HAMAP and UniRule. Not experimental. "Sequence model evidence. We use the ECO code ECO:0000255 in manual assertions and ECO:0000256 (and its descendant code ECO:0000259 ) in automatic assertions, respectively, for information which has been generated by the UniProtKB automatic annotation system. The database name and identifier of the sequence model rule that has been used by the system are provided. The ECO code ECO:0000255 is also used for information which has been generated by various sequence analysis programs that are used during the manual curation process and which has been verified by a curator."
    # 40888  (13498 human)   ECO:0000269     KEEP Experimental small-scale evidence. "Evidence that is used only in manual assertions. Experimental evidence. We use the ECO code ECO:0000269 for manually curated information for which there is published experimental evidence. The PubMed identifier of the publication(s) which is the original source of the information is provided (for publications that are not in PubMed we indicate instead the corresponding UniProtKB reference number(s))."
    # 2458     (609 human)   ECO:0000305     KEEP Curator-added. Seems to be based on very reasonable experimental evidence (e.g. S191 in ENAM_HUMAN from PMID 25789606, which explicitly states S191 is phosphorylated). "Curator inference evidence. We use the ECO code ECO:0000305 for manually curated information which has been inferred by a curator based on his/her scientific knowledge or on the scientific content of an article. If the information is inferred from the scientific content of an article, the PubMed identifier of the supporting paper is provided (for publications that are not in PubMed we indicate instead the corresponding UniProtKB reference number(s))."
    # 92         (5 human)   ECO:0000303     X SKIP No evidence. "Non-traceable author statement evidence. We use the ECO code ECO:0000303 for manually curated information that is based on statements in scientific articles for which there is no experimental support. The PubMed identifier of the paper(s) which is the original source of the information is provided (for publications that are not in PubMed we indicate instead the corresponding UniProtKB reference number(s))."
    # 48        (18 human)   ECO:0000312     KEEP From another database. Using this myself for dbPTM etc. Always coincides with ECO:0000269 (small-scale experimental evidence) in UniProt, except for one case, https://www.uniprot.org/uniprotkb/K7N5M5/entry#ptm_processing, where the other database is PDB (perfect). "Imported information evidence. We use the ECO code ECO:0000312 in manual assertions for information which has been imported from another database. The database name and identifier of the entry from which the information has been imported are provided."
    # 
    # # Grouped by 'aa', 'type', 'description' fields
    # SELECT * FROM (
    # SELECT NULL          AS evid, 'skip' AS `skip`, aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE evid IS NULL              GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000250' AS evid, 'skip' AS `skip`, aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE evid LIKE '%ECO:0000250%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000255' AS evid, 'skip' AS `skip`, aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE evid LIKE '%ECO:0000255%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000269' AS evid, 'ok' AS `skip`,   aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE evid LIKE '%ECO:0000269%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000303' AS evid, 'skip' AS `skip`, aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE evid LIKE '%ECO:0000303%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000305' AS evid, 'ok' AS `skip`,   aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE evid LIKE '%ECO:0000305%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000312' AS evid, 'ok' AS `skip`,   aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE evid LIKE '%ECO:0000312%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0007744' AS evid, 'ok' AS `skip`,   aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE evid LIKE '%ECO:0007744%' GROUP BY aa, type, description
    # ) AS t ORDER BY t.namesites DESC, t.names DESC, t.evid, t.aa, t.type, t.description;
    # 
    # # Grouped by 'aa', 'type', 'description' fields (human only)
    # SELECT * FROM (
    # SELECT NULL          AS evid, 'skip' AS `skip`, aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE species='human' AND evid IS NULL              GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000250' AS evid, 'skip' AS `skip`, aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE species='human' AND evid LIKE '%ECO:0000250%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000255' AS evid, 'skip' AS `skip`, aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE species='human' AND evid LIKE '%ECO:0000255%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000269' AS evid, 'ok' AS `skip`,   aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE species='human' AND evid LIKE '%ECO:0000269%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000303' AS evid, 'skip' AS `skip`, aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE species='human' AND evid LIKE '%ECO:0000303%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000305' AS evid, 'ok' AS `skip`,   aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE species='human' AND evid LIKE '%ECO:0000305%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0000312' AS evid, 'ok' AS `skip`,   aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE species='human' AND evid LIKE '%ECO:0000312%' GROUP BY aa, type, description UNION
    # SELECT 'ECO:0007744' AS evid, 'ok' AS `skip`,   aa, type, description, COUNT(DISTINCT name) AS names, COUNT(DISTINCT name, site) AS namesites, COUNT(DISTINCT pmids), GROUP_CONCAT(DISTINCT pmids ORDER BY pmids SEPARATOR '|'), COUNT(DISTINCT species), GROUP_CONCAT(DISTINCT species ORDER BY species SEPARATOR ',') FROM unimod WHERE species='human' AND evid LIKE '%ECO:0007744%' GROUP BY aa, type, description
    # ) AS t ORDER BY t.namesites DESC, t.names DESC, t.evid, t.aa, t.type, t.description;
    # 
    # >> To keep: ECO:0007744, ECO:0000269, ECO:0000305, ECO:0000312

    # -hc (experimental only)
    # $mainquery = Query("SELECT type, description, COUNT(DISTINCT name, site) AS c FROM unimod WHERE (evid LIKE '%ECO:0000006%' OR evid LIKE '%ECO:0000269%' OR evid LIKE '%ECO:0007744%' OR evid LIKE '%ECO:0007829%') GROUP BY type, description ORDER BY c DESC");
    # $mainquery = Query("SELECT aa, type, description, COUNT(DISTINCT name, site) AS c FROM unimod WHERE (evid LIKE '\%ECO:0007744%' OR evid LIKE '\%ECO:0000269%' OR evid LIKE '\%ECO:0000305%' OR evid LIKE '\%ECO:0000312%') GROUP BY type, description ORDER BY c DESC");
    $evid = " AND (evid LIKE '\%ECO:0007744%' OR evid LIKE '\%ECO:0000269%' OR evid LIKE '\%ECO:0000305%' OR evid LIKE '\%ECO:0000312%')";
}
else
{
    # Including potential/probable/by similarity
    # $mainquery = Query("SELECT aa, type, description, COUNT(DISTINCT name, site) AS c FROM unimod WHERE GROUP BY type, description ORDER BY c DESC");
    $evid = '';
}




# Set ptm column (apply PTM labels). Also set e.g. ptm='K-ac' if the description is 'K-ac'.
$affected = 0;
foreach $key (nsort(keys(%ptm)))
{
    ($aa, $type, $desc) = split(/\|/, $key);
    $ptm = $ptm{$key};

    print " >> $aa >> $type >> $desc";

    if (switch('hc'))
    {
        if (!switch('debug'))
        {
            # Update
            $query = Query("UPDATE unimod SET ptm='$ptm' WHERE aa='$aa' AND ((type='$type' AND description LIKE '".esc($desc)."') OR description='$ptm') AND (evid LIKE '\%ECO:0007744%' OR evid LIKE '\%ECO:0000269%' OR evid LIKE '\%ECO:0000305%' OR evid LIKE '\%ECO:0000312%')");
        }
        else
        {
            # -debug: print only 
            $query = Query("SELECT id FROM unimod WHERE aa='$aa' AND ((type='$type' AND description LIKE '".esc($desc)."') OR description='$ptm') AND (evid LIKE '\%ECO:0007744%' OR evid LIKE '\%ECO:0000269%' OR evid LIKE '\%ECO:0000305%' OR evid LIKE '\%ECO:0000312%')");
        }
    }
    else
    {
        if (!switch('debug'))
        {
            # Update
            $query = Query("UPDATE unimod SET ptm='$ptm' WHERE aa='$aa' AND ((type='$type' AND description LIKE '".esc($desc)."') OR description='$ptm')");
        }
        else
        {
            # -debug: print only 
            $query = Query("SELECT id FROM unimod WHERE aa='$aa' AND ((type='$type' AND description LIKE '".esc($desc)."') OR description='$ptm')");
        }
    }

    print " >> ".commify(Numrows($query))." sites ($ptm)\n";

    addme("total aa|type|descriptions", "$aa|$type|$desc");
    if (Numrows($query) > 0)
    {
        addme("ptm label successfully added for ptm label", $ptm);
        addme("ptm label successfully added for aa|type|description", "$aa|$type|$desc");
        
        $affected += Numrows($query);
    }
}
state("Rows affected: ".commify($affected));





# Set subsets
state("Setting subset '_small'...");
$query = Query("UPDATE unimod SET subset=CONCAT(source, '_small') WHERE scale IN ('small', 'both')");
state("Rows affected: ".Numrows($query));

state("Setting subset '_large'...");
$query = Query("UPDATE unimod SET subset=CONCAT(source, '_large') WHERE scale IN ('large')");
state("Rows affected: ".Numrows($query));

state("Setting subset to source where scale is NULL...");
$query = Query("UPDATE unimod SET subset=source WHERE scale IS NULL");
state("Rows affected: ".Numrows($query));






# Show very common, but still unlabelled PTMs (≥$minsites, i.e. ≥1000 sites)
# $mainquery = Query("SELECT aa, type, description, COUNT(DISTINCT name, site) AS c FROM unimod WHERE ptm=''$evid GROUP BY type, description HAVING c >= $minsites ORDER BY c DESC");
$mainquery = Query("SELECT aa, type, description, COUNT(DISTINCT name, site) AS c FROM unimod WHERE ptm IS NULL$evid GROUP BY type, description HAVING c >= $minsites ORDER BY c DESC");
# $affectedcomplete = 0;
# $affectedpartial = 0;
if (Numrows($mainquery) > 0)
{
    warn("\n\nWarning: There are still unlabelled PTMs with ≥$minsites sites, perhaps add labels for these!\n\n");
}
while (($aa, $type, $desc, $sites) = Fetch($mainquery))
{
    # last if ($sites <= $minsites);
    print " >> $aa >> $type >> $desc >> $sites sites >> Still unlabelled!\n";
    # warn("Warning: >> $aa >> $type >> $desc >> $sites sites >> Still unlabelled!\n");
    
    # # Figure out which PTM label to use
    # $ptm = '';
    # foreach $desc (split(/; /, $desc))
    # {
    #     # print "   >> $desc\n";
    #     $key = $type.'|'.$desc;
    #     
    #     if (exists($ptm{$key}))
    #     {
    #         die("Error: Multiple PTM labels found for type '$type' description '$desc'") if ($ptm ne '');
    #         $ptm = $ptm{$key};
    #     }
    # }
    # if ($ptm eq '')
    # {
    #     if (switch('other'))
    #     {
    #         # Note - Sometimes both type and description are ''. This is from PhosphoSitePlus and is OK, I'm just not parsing the description from it.
    #         $ptm = 'other';
    #     }
    #     else
    #     {
    #         next;
    #     }
    # }
    # 
    #  
    #     # Update unimod table
    #     print " >> $type >> $desc >> $sites sites >> $ptm \n";
    # 
    #     if (switch('hc'))
    #     {
    #         # print "       >> Updating '$type' '$desc' to '$ptm' (hc only)...\n";
    #         if (!switch('debug'))
    #         {
    #     		# $query = Query("UPDATE unimod SET ptm='$ptm' WHERE type='$type' AND description='".esc($desc)."' AND (evid LIKE '\%ECO:0000006\%' OR evid LIKE '\%ECO:0000269\%');");
    #     		$query = Query("UPDATE unimod SET ptm='$ptm' WHERE type='$type' AND description='".esc($desc)."' AND (evid LIKE '\%ECO:0007744%' OR evid LIKE '\%ECO:0000269%' OR evid LIKE '\%ECO:0000305%' OR evid LIKE '\%ECO:0000312%');");
    #         }
    #         else
    #         {
    #             # state("UPDATE unimod SET ptm='$ptm' WHERE type='$type' AND description='".esc($desc)."' AND (evid LIKE '\%ECO:0000006\%' OR evid LIKE '\%ECO:0000269\%');");
    #             # $query = Query("SELECT id FROM unimod WHERE type='$type' AND description='".esc($desc)."' AND (evid LIKE '\%ECO:0000006\%' OR evid LIKE '\%ECO:0000269\%');");
    #             $query = Query("SELECT id FROM unimod WHERE type='$type' AND description='".esc($desc)."' AND (evid LIKE '\%ECO:0007744%' OR evid LIKE '\%ECO:0000269%' OR evid LIKE '\%ECO:0000305%' OR evid LIKE '\%ECO:0000312%');");
    #         }
    #         # print "         >> ".Numrows($query)."\n";
    # 
    # 		# This can happen when there are multiple records of the same modification for the same residue from the same source, e.g. for "SELECT * FROM unimod WHERE name='MBP_HUMAN' AND description='N-acetylalanine';".
    # 		# This is why I always do SELECT DISTINCT name, site.
    #         # die("Error: Expected to change $sites sites, but numrows was only ".Numrows($query)) if ($sites != Numrows($query));
    #     }
    #     else
    #     {
    #         # print "       >> Updating '$type' '$desc' to '$ptm' (including potential/probable/by similarity)...\n";
    #         if (!switch('debug'))
    #         {
    # 		    $query = Query("UPDATE unimod SET ptm='$ptm' WHERE type='$type' AND description='".esc($desc)."'");
    # 	    }
    # 	    else
    # 	    {
    #             # state("UPDATE unimod SET ptm='$ptm' WHERE type='$type' AND description='".esc($desc)."'");
    # 		    $query = Query("SELECT id FROM unimod WHERE type='$type' AND description='".esc($desc)."'");
    # 	    }
    #         # print "         >> ".Numrows($query)."\n";
    # 
    # 		# This can happen when there are multiple records of the same modification for the same residue from the same source, e.g. for "SELECT * FROM unimod WHERE name='MBP_HUMAN' AND description='N-acetylalanine';".
    # 		# This is why I always do SELECT DISTINCT name, site.
    #         # die("Error: Expected to change $sites sites, but numrows was only ".Numrows($query)) if ($sites != Numrows($query));
    #     }
    #     
    #     addme("total type|descriptions", "$type|$desc");
    #     if (Numrows($query) > 0)
    #     {
    #         addme("ptm label successfully added for ptm label", $ptm);
    #         addme("ptm label successfully added for type|description", "$type|$desc");
    #         $affected += Numrows($query);
    #         
    #         # if ($desc eq $desc)
    #         # {
    #         #     addme("ptm label successfully added due to a standard, complete match for type|description", "$type|$desc");
    #         #     $affected += Numrows($query);
    #         #     $affectedcomplete += Numrows($query);
    #         # }
    #         # else
    #         # {
    #         #     addme("ptm label successfully added due to a partial match for type|description", "$type|$desc");
    #         #     $affected += Numrows($query);
    #         #     $affectedpartial += Numrows($query);
    #         # }
    #     }
}
stoptime();

showmeall(1);
# showmesome(50);
nl();

showme("ptm label successfully added for ptm label");
# showme("ptm label successfully added for type|desc");

if (!switch('debug'))
{
    Optimize('unimod');
}

done();
