#!/usr/bin/env perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');


# $rel = 105;				# December 2021 (used in UniProt 2022_01) [probably according to ~/update/uniprot/uniens_seqmatch.pl -debug]
# $rel = 106;				# April 2022 (used in UniProt 2022_02)
# $rel = 107;				# July 2022 (presumably will be used in the upcoming UniProt 2022_03)
$rel = 108;					# October 2022
# $plants_rel = 42;



# Download species annotation (taxon IDs etc.)
# run("Download species annotation (taxon IDs etc.)", "wget -nv -O input-uniprot_report_EnsemblVertebrates.txt 'http://ftp.ensembl.org/pub/release-$rel/uniprot_report_EnsemblVertebrates.txt'");
run("Download species annotation (taxon IDs etc.)", "wget -nv -O 'input-uniprot_report_EnsemblVertebrates-$rel.txt' 'http://ftp.ensembl.org/pub/release-$rel/uniprot_report_EnsemblVertebrates.txt'");


# All species
run("Make 'input_$rel' directory", "mkdir -p 'input_$rel'");
# run("Download proteomes", "wget -nv -P 'input_$rel' -r -N -nd -A '*.pep.all.fa.gz' 'ftp://ftp.ensembl.org/pub/release-$rel/fasta'");
# run("Download CDS", "wget -nv -P 'input_$rel' -r -N -nd -A '*.cds.all.fa.gz' 'ftp://ftp.ensembl.org/pub/release-$rel/fasta'");
run("Download proteomes & CDS", "wget -nv -P 'input_$rel' -r -N -nd -A '*.pep.all.fa.gz' -A '*.cds.all.fa.gz' 'ftp://ftp.ensembl.org/pub/release-$rel/fasta'");
# run("Unpack proteomes", "gunzip -f input/*.pep.all.fa.gz");

if (($rel == 106) or ($rel == 107) or ($rel == 108))
{
	run("Delete 'cyprinus_carpio' if release is 106 or 107 or 108 and possibly higher (superseded by 'cyprinus_carpio_carpio', which it clashes with in terms of ENSPs etc.)", "rm -fv input_$rel/Cyprinus_carpio.common_carp_genome.pep.all.fa.gz");
	run("Delete 'cyprinus_carpio' if release is 106 or 107 or 108 and possibly higher (superseded by 'cyprinus_carpio_carpio', which it clashes with in terms of ENSPs etc.)", "rm -fv input_$rel/Cyprinus_carpio.common_carp_genome.cds.all.fa.gz");
}

# Problem in 107 and 108: Cebus_imitator and Cebus_capucinus seem to be the same species, but are present as two files. Cebus_imitator isn't listed in uniprot_report_EnsemblVertebrates.txt, although it seems to be the more recent one:
# 106 has Cebus_capuchinus only.

# https://en.wikipedia.org/wiki/Panamanian_white-faced_capuchin
# "Until the 21st century the Panamanian white-faced capuchin was considered conspecific with Cebus capucinus, the Colombian white-faced capuchin, but as a separate subspecies C. capucinus imitator"
# "C. imitator and C. capucinus split up to 2 million years ago"
# So: "Cebus capucinus imitator" (= Cebus capucinus here) is now incorrect, and it should be "Cebus imitator".

# Cebus_capucinus.Cebus_imitator-1.0.pep.all.fa.gz and Cebus_imitator.Cebus_imitator-1.0.pep.all.fa.gz contain exactly the same set of ENSPVs, and in the same order.
# Only the sequence headers differ (the sequences are exactly the same, and in the same order).

# The headers in Cebus_imitator have more information than those in Cebus_capucinus.
# The Cebus_capucinus file is one day newer.

# 17585  zcat Cebus_capucinus.Cebus_imitator-1.0.pep.all.fa.gz | g -o "^>ENS\w+\.\d+" > ../tmp_capucinus.txt
# 17586  zcat Cebus_imitator.Cebus_imitator-1.0.pep.all.fa.gz | g -o "^>ENS\w+\.\d+" > ../tmp_imitator.txt
# 17587  diff ../tmp_capucinus.txt ../tmp_imitator.txt 
# 
# 17610  zcat Cebus_capucinus.Cebus_imitator-1.0.pep.all.fa.gz | g -c gene_symbol:
# 17611  zcat Cebus_imitator.Cebus_imitator-1.0.pep.all.fa.gz | g -c gene_symbol:
# 
# 17592  zcat Cebus_capucinus.Cebus_imitator-1.0.pep.all.fa.gz | g -v "^>" > ../tmp_seq_capucinus.txt
# 17593  zcat Cebus_imitator.Cebus_imitator-1.0.pep.all.fa.gz | g -v "^>" > ../tmp_seq_imitator.txt
# 17594  diff ../tmp_seq_capucinus.txt ../tmp_seq_imitator.txt 

# There are no differences except in the titles. In the titles, only Cebus_imitator has "description:" entries, but I don't use these. I only use "gene_symbol:".
# Cebus_imitator only has 2 proteins with "gene_symbol:" annotation. (APOA1, twice)
# Cebus_capucinus has 55 proteins with "gene_symbol:" annotation. (CYP4F12 and other CYPs (all CYPs))

# 17612  zcat Cebus_capucinus.Cebus_imitator-1.0.pep.all.fa.gz | g "^>" | perl -npe 's/ description:.+//; s/ gene_symbol:\S+//g;' > ../tmp_titles_capucinus.txt
# 17613  zcat Cebus_imitator.Cebus_imitator-1.0.pep.all.fa.gz | g "^>" | perl -npe 's/ description:.+//; s/ gene_symbol:\S+//g;' > ../tmp_titles_imitator.txt
# 17614  diff ../tmp_titles_capucinus.txt ../tmp_titles_imitator.txt 


# Finally: Cebus_capucinus is the name that's used in the Compara species tree in Compara 108.

# >> Go with Cebus_capucinus, delete Cebus_imitator.

# Since I do actually use the gene_symbols, but not the descriptions: I'll keep everything untouched and parse Cebus_capucinus. It's also the newer file.
if (($rel == 107) or ($rel == 108))
{
	# Unlike cyprinus_carpio, there is no need to delete cebus_imitator from the ensembl_species table since only cebus_capucinus is present there
	run("Delete 'cebus_imitator' if release is 108 and possibly higher (clashes with 'cebus_capucinus' in terms of ENSPs etc. with all-identical sequences, and it's the same species, cebus_imitator actually being the newer name)", "rm -fv input_$rel/Cebus_imitator.Cebus_imitator-1.0.pep.all.fa.gz");
	run("Delete 'cebus_imitator' if release is 108 and possibly higher (clashes with 'cebus_capucinus' in terms of ENSPs etc. with all-identical sequences, and it's the same species, cebus_imitator actually being the newer name)", "rm -fv input_$rel/Cebus_imitator.Cebus_imitator-1.0.cds.all.fa.gz");
}

# show directory
run("ls", "ls input_$rel");
run("Number of files", "ls -1U input_$rel | wc -l");
state("Downloaded ENSEMBL release '$rel'");

done();
