#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# # Download
# run("Download UniProt XML", "download.pl");
# # # run("Download PHOSIDA data", "download_phosida.pl");	# PHOSIDA was last updated 2011-09, too old
# # # run("Download Phospho.ELM data", "download_phospho_elm.pl");	# PhosphoELM was last updated 2011-11, too old
# run("Download PhosphoSitePlus data", "download_phosphositeplus.pl");
# # # run("Split UniProt XML into chunks", "chunks.pl");



# Clear tables
Clear('uniacc');
Clear('uniens');
Clear('uniens_uniprot');
Clear('unifeat');
Clear('unifunc');
Clear('unigo');
Clear('unihost');
Clear('uniid');
Clear('uniid_counts');
Clear('uniinterpro');
Clear('uniiso');
Clear('unikey');
Clear('uniloc');
Clear('unimod');
Clear('unimod_control');
Clear('unipdb');
Clear('unipfam');
Clear('unipfam_pfam');
Clear('uniprot');
Clear('uniproteome');
# Clear('unirefprot');	# Note - unirefprot is updated by its own pipeline at update/unirefprot. The reference proteomes aren't particularly useful, though (they don't contain all of Swiss-Prot for e.g. human).
Clear('uniseq');
Clear('unisnp');
Clear('unisupfam');
Clear('unitax');

# warn("WARNING: Didn't clear 'uniseq', only deleted the type='UniProt' entries (to preserve disorder predictions)!");
# Query("DELETE FROM uniseq WHERE type='UniProt'");




# Read UniProt XMLs
# All species
run("Fill tables from UniProtKB/Swiss-Prot XML file (all species)", "main.pl input/uniprot_sprot.xml.gz");
# exit;

# # Taxonomic divisions
#
# # Human
# run("Fill tables from UniProtKB/Swiss-Prot XML file (human)", "main.pl input/uniprot_sprot_human.xml.gz -humanonly");
# # run("Fill tables from UniProtKB/TrEMBL XML file (human)", "main.pl input/uniprot_trembl_human.xml.gz -humanonly");
#
# # Viruses
# run("Fill tables from UniProtKB/Swiss-Prot XML file (viruses, human host only)", "main.pl input/uniprot_sprot_viruses.xml.gz -humanonly");
# # run("Fill tables from UniProtKB/TrEMBL XML file (viruses, human host only)", "main.pl input/uniprot_trembl_viruses.xml.gz -humanonly");
#
# # Bacteria
# # run("Fill tables from UniProtKB/Swiss-Prot XML file (bacteria, human host only)", "main.pl input/uniprot_sprot_bacteria.xml.gz -humanonly");	# Swiss-Prot bacteria: no hosts annotated, skipping
# # run("Fill tables from UniProtKB/TrEMBL XML file (bacteria, human host only)", "main.pl input/uniprot_trembl_bacteria.xml.gz -humanonly");	# TrEMBL Bacteria: no hosts annotated, and also too huge (~100 GB as .xml.gz), skipping




# # Add varsplic isoform sequences into table uniiso
# No longer needed since this is integrated in main.pl now (it parses the varsplic fasta file)
# # run("Add sequences to table 'uniiso' from varsplic isoform FASTA file", "uniiso_seq.pl -debug");
# run("Add sequences to table 'uniiso' from varsplic isoform FASTA file", "uniiso_seq.pl");
# # run("Add sequences to table 'uniiso' from varsplic isoform FASTA file", "uniiso_seq.pl -humanonly");

# # Read ID Mapping
# # run("Fill table 'uniid' from ID Mapping file", "uniid.pl");
# run("Fill table 'uniid' from ID Mapping file", "uniid.pl -humanonly");
# run("Fill table 'uniid_counts' (helps to select ID types)", "uniid_counts.pl");
#
# # Read Pfam information on UniProt proteins (may always be somewhat dated)
# Superseded by InterPro, which is up-to-date with UniProt, unlike Pfam's own irregular releases
# run("Fill table 'unipfam_pfam' from Pfam file", "unipfam_pfam.pl");


# Extract species sequences from 'uniseq' to FASTA file for disorder prediction (DISOPRED2, IUPred, MULTICOM)
run("Extract FASTA for disorder predictors", "extract_fasta.pl human");
# run("Extract FASTA for disorder predictors", "extract_fasta.pl arath");
# run("Extract FASTA for disorder predictors", "extract_fasta.pl drome");
# run("Extract FASTA for disorder predictors", "extract_fasta.pl ecoli");
# run("Extract FASTA for disorder predictors", "extract_fasta.pl mouse");
# run("Extract FASTA for disorder predictors", "extract_fasta.pl salty");
# run("Extract FASTA for disorder predictors", "extract_fasta.pl yeast");



# Fill unihost host species mnemonic field
# state("Filling host species mnemonic field in table 'unihost' from same table via taxon ID...\n$q\n");
# $q = "UPDATE unihost h1 JOIN unihost h2 ON h1.hosttax=h2.hosttax SET h1.host=h2.host WHERE h1.host IS NULL AND h2.host IS NOT NULL";
nl();
state("Filling host species mnemonic field in table 'unihost' from table 'unitax' via taxon ID...\n$q\n", 1);
$q = "UPDATE unihost h JOIN unitax t ON h.hosttax=t.tax SET h.host=t.species WHERE h.host IS NULL";
$query = Query($q);
state(commify(Numrows($query))." rows affected", 1);

$query = Query("SELECT id FROM unihost WHERE host IS NULL");
if (Numrows($query) > 0)
{
	# warn("Warning: There are still ".commify(Numrows($query))." rows in table 'unihost' with host=NULL (shouldn't happen, may need to do the update query via table 'unitax' instead)");
	state("Note: There are still ".commify(Numrows($query))." rows in table 'unihost' with host=NULL. These are hosts that don't have any Swiss-Prot proteins. Sometimes they're higher-level taxonomic terms such as 'genus' (e.g. tax 6937 are the Ornithodoros ticks). This is safe to ignore.", 1);
}
nl();


# Add more PTMs - Large databases
# Test run without actually inserting anything (and -quiet so it doesn't print out every insert query)



# PhosphoSitePlus
# run("PhosphoSitePlus (without inserting anything)", "unimod_phosphositeplus.pl -debug -quiet");
run("Add PTM sites from PhosphoSitePlus", "unimod_phosphositeplus.pl");



# dbPTM
run("Add PTM sites from dbPTM", "unimod_dbptm.pl");


# Ochoa et al.
run("Import Ochoa et al. into unimod_ochoa_full... tables", "~/scripts/import.pl input/unimod_ochoa/unimod_ochoa_table_s2_annotated_phosphoproteome_features.csv unimod_ochoa_full -overwrite");
run("Import Ochoa et al. into unimod_ochoa_full... tables", "~/scripts/import.pl input/unimod_ochoa/unimod_ochoa_table_s3_phosphosite_functional_scores.csv unimod_ochoa_full_score -overwrite");
run("Import Ochoa et al. into unimod_ochoa_full... tables", "~/scripts/import.pl input/unimod_ochoa/unimod_ochoa_table_s5_clinvar_variants.csv unimod_ochoa_full_clinvar -overwrite");
run("Import Ochoa et al. into unimod_ochoa_full... tables", "~/scripts/import.pl input/unimod_ochoa/unimod_ochoa_table_s6_thermal_proteome_profiling.csv unimod_ochoa_full_yeast_thermal -overwrite");
run("Add PTM sites from table 'unimod_ochoa_full' to table 'unimod'", "unimod_ochoa.pl");


# Phospho.ELM
run("Phospho.ELM", "unimod_phospho_elm.pl -debug -quiet");          # 16,868 known,  11,329 new
# PHOSIDA
run("PHOSIDA", "unimod_phosida.pl -debug -quiet");                  # 2,446 known,   4,885 new

# # # Add more PTMs - Individual studies
# # These are the only datasets of their kind in these species:
# run("Pang et al.", "unimod_pang.pl -debug -quiet");                 # 0 known,      76 new (yeast K,R-me)
# run("Wang et al.", "unimod_wang.pl -debug -quiet");                 # 0 known,      250 new (salmonella K-ac)
# run("Weinert et al.", "unimod_weinert.pl -debug -quiet");               # 2 known,      799 new (drosophila K-ac)
# run("Zhang et al.", "unimod_zhang.pl -debug -quiet");                   # 0 known,      0 new (e. coli K-suc)
# # Skipping these because they are very minor additions compared to what I've already got (probably just introducing false matches)
# # run("Danielsen et al.", "unimod_danielsen.pl -debug -quiet");           # 37 known,     582 new (human K-ub) (best of these four, but apparently already in one of the big databases)
# # # run("Xu et al.", "unimod_xu.pl -debug -quiet");                     # 44 known,     338 new (human K-ub)
# # # run("Zhao et al.", "unimod_zhao.pl -debug -quiet");                     # 98 known,     612 new (human K-ac)
# # # run("Zielinska et al.", "unimod_zielinska.pl -debug -quiet");           # 3797 known,   936 new (mouse N-gly) (definitely already in PHOSIDA)




# # Fill unimod 'aa' field
# nl();
# $q = "UPDATE unimod SET aa=NULL";
# state("Clearing 'aa' column in table 'unimod'...\n$q\n", 1);
# $query = Query($q);
# state(commify(Numrows($query))." rows affected", 1);

nl();
$q = "UPDATE unimod m, uniseq s SET m.aa=SUBSTRING(s.seq, m.site, 1) WHERE m.acc=s.acc AND s.type IN ('UniProt', 'UniIso')";
state("Filling aa column in table 'unimod' from table 'uniseq'...\n$q\n", 1);
$query = Query($q);
state(commify(Numrows($query))." rows affected", );



# Fill 'ptm' field + filter out unwanted PTMs
run("Fill 'ptm' field, leave out non-experimental", "unimod_filter.pl -hc");
# run("Fill 'ptm' field, leave out non-experimental, and use 'other' rather than blank for PTMs other than K-ac etc.", "unimod_filter.pl -hc -other");

# # Fill 'unimod_control' with non-PTM control residues
# run("Fill unimod_control table", "unimod_control.pl human hc all 0");

# uniens: Now moved to ~/update/ensembl

# Optimize tables
run("Optimize tables", "optimize.pl");

done();
