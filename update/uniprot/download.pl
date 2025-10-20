#!/usr/bin/env perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');




# Swiss-Prot XML (all species)
run("Download", "wget -O input/uniprot_sprot.xml.gz 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz'");
# run("Unpack", "gunzip -f input/uniprot_sprot.xml.gz");

# Swiss-Prot XML (human)
run("Download", "wget -O input/uniprot_sprot_human.xml.gz 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.xml.gz'");
# run("Unpack", "gunzip -f input/uniprot_sprot_human.xml.gz");

# TrEMBL XML (human)
run("Download", "wget -O input/uniprot_trembl_human.xml.gz 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_human.xml.gz'");
# run("Unpack", "gunzip -f input/uniprot_trembl_human.xml.gz");

# ID mapping (human)
# Note: Not entirely sure yet whether to use _selected or the normal one. Previously, I always used the normal idmapping.dat.gz. So, going with normal so it's the correct format (_selected.tab has a different format).
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README has some explanations but they're a bit unclear.
# Size is comparable (_selected is 10 GB, full is 20 GB)
run("Download", "wget -O input/HUMAN_9606_idmapping.dat.gz 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'");
# run("Download", "wget -O input/HUMAN_9606_idmapping_selected.tab.gz 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz'");



# Swiss-Prot splice variants (additional sequences, FASTA only)
# Still need these sequences, even for human-only. The human XML files do not contain the isoforms' sequences.
# ALso note: There is no "varsplic" file for TrEMBL, only for Swiss-Prot.
run("Download", "wget -O input/uniprot_sprot_varsplic.fasta.gz 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz'");
# run("Download", "wget -O input/uniprot_sprot_varsplic.fasta.gz 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz'");
# Unpack it (it's tiny)
run("Unpack", "gunzip -f input/uniprot_sprot_varsplic.fasta.gz");

# ID mapping (includes Swiss-Prot & TrEMBL)
run("Download", "wget -O input/idmapping.dat.gz 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz'");




# show directory
run("ls", "ls input");

done();
