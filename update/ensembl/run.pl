#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
# $rel = 106;
# $rel = 107;
$rel = 108;

# run("Download", "download.pl");

# Import species table from Ensembl FTP server to map Ensembl species names to taxon IDs, and then to UniProt "unispec" species mnemonics (5 letters)
run("Import ensembl_species table", "~/scripts/import.pl input-uniprot_report_EnsemblVertebrates-$rel.txt ensembl_species -keepcomments -overwrite -quiet");

# Add unispec column to table structure, and rename taxonomy_id to tax
Query("ALTER TABLE ensembl_species ADD COLUMN `unispec` VARCHAR(25) NULL DEFAULT NULL AFTER `species`, ADD INDEX `Unispec` (`unispec` ASC) VISIBLE;");
Query("ALTER TABLE ensembl_species CHANGE COLUMN `taxonomy_id` `tax` INT NULL DEFAULT NULL AFTER `species`, CHANGE COLUMN `species` `fullspecies` VARCHAR(250) NULL DEFAULT NULL, CHANGE COLUMN `unispec` `species` VARCHAR(25) NULL DEFAULT NULL;");

# Manual fix for an Ensembl 106 and 107 (but not 105) problem: Cyprinus_carpio clashes with the newer Cyprinus_carpio_carpio (and only the latter is included on the website, for a total of 314 species, one less than is on the FTP server and in ensembl_species)
# Still not fixed in 108 even though I contacted Ensembl's support about it and they said it was an error. Should be fixed at some point (perhaps in 109)
if (($rel == 106) or ($rel == 107) or ($rel == 108))
{
	$query = Query("DELETE FROM ensembl_species WHERE fullspecies='cyprinus_carpio'");
	if (Numrows($query) > 0)
	{
		warn("Warning: DEBUG: Manually deleted 'cyprinus_carpio' from table 'ensembl_species' because: manually skipped species 'Cyprinus_carpio' since it wasn't on the species list at https://useast.ensembl.org/info/about/species.html, and it clashes with 'Cyprinus_carpio_carpio' (same ENSP identifiers etc.)\n");
	}
}

# # Manual fix:
# Not anymore, see download.pl for comments
# # Rename cebus_capucinus to cebus_imitator
# if (($rel == 108))
# {
# 	$query = Query("UPDATE ensembl_species SET fullspecies='cebus_imitator' WHERE fullspecies='cebus_capucinus'");
# 	if (Numrows($query) > 0)
# 	{
# 		warn("Warning: DEBUG: Manually renamed 'cebus_capucinus' to 'cebus_imitator' in table 'ensembl_species' because: that's what it's called on the species list at https://useast.ensembl.org/info/about/species.html, and on https://www.uniprot.org/taxonomy/2715852, and it would clash with Cebus_imitator's file (same ENSP identifiers etc.)\n");
# 	}
# }



# Parse Ensembl
run("Read ENSEMBL IDs and sequences", "main.pl $rel");

# Add UniProt species mnemonics (5 letters)
run("Add short unispec species names to 'ensembl' and 'ensembl_species' tables", "species.pl");

# Note: Gene symbols are parsed from the FASTA file now (gene_symbol:... goes into 'symbol' in the table)

# These manual downloads here are not necessary anymore:
# run("Add HGNC symbols and IDs for human", "symbols.pl human");
# run("Add MGI symbols and IDs for mouse", "symbols.pl mouse");
# run("Add ZFIN symbols and IDs for zebrafish", "symbols.pl zebrafish");

# Ensembl doesn't have mappings to UniGene anymore
# run("Add Unigene and NCBI/Entrez Gene IDs for human", "unigene.pl human");



# uniens: Moved from ~/update/uniprot
# # # Fill uniens table (sequence-based) from scratch, so it only contains cases where the UniProt and Ensembl sequences really match.
# # # Another reason not to use UniProt's mapping: it contains NM_... and NP_... IDs for some species (e.g. BOMMO in 2022_04). Ensembl doesn't use these.
# # # UniProt's own (outdated) mapping is retained in uniens_uniprot.
# # Random isoform examples to look into:
# # uniens tests:
# SELECT * FROM uniens WHERE ensp='ENSP00000003583';
# SELECT * FROM uniens WHERE acc='A2AIG8';
# SELECT * FROM uniens WHERE acc='A2AIG8-2';
# SELECT * FROM uniens WHERE name IN ('ECE2_MOUSE', 'EFCE2_MOUSE');
# SELECT * FROM uniens WHERE acc='B2RQR8-2';

# run("Fill uniens table", "uniens.pl -debug");
run("Fill uniens table", "uniens.pl");


starttime();
Optimize('ensembl');
Optimize('ensembl_species');
Optimize('uniens');
stoptime();

run("Test whether cds translates to seq in table 'ensembl' (human only)", "test_cds_seq.pl human");
run("Test whether cds translates to seq in table 'ensembl' (all species)", "test_cds_seq.pl all");

done();
