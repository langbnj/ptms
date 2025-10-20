#!/usr/bin/env perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');

# $rel = 80;
# $rel = 106;	# April 2022
# $rel = 107;	# July 2022
$rel = 108;		# October 2022


# PAL2NAL:
# http://www.bork.embl.de/pal2nal/pal2nal.pdf
# http://www.bork.embl.de/pal2nal/
# http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz


# Download
cd("input");

# Get human orthologs (in 199 other species in release 107, i.e. there is a total of 200 Compara species including human)
run("Download", "wget -nv 'ftp://ftp.ensembl.org/pub/release-$rel/mysql/ensembl_mart_$rel/hsapiens_gene_ensembl__homolog_*__dm.txt.gz'");
# run("Download", "wget -nv 'http://ftp.ensembl.org/pub/release-$rel/mysql/ensembl_mart_$rel/hsapiens_gene_ensembl__homolog_*__dm.txt.gz'");
run("Unpack", "gunzip -f hsapiens_gene_ensembl__homolog_*__dm.txt.gz");

# Alternative file:
# http://ftp.ensembl.org/pub/release-108/tsv/ensembl-compara/homologies/Compara.108.protein_default.homologies.tsv.gz
# Example:
# gene_stable_id	protein_stable_id	species	identity	homology_type	homology_gene_stable_id	homology_protein_stable_id	homology_species	homology_identity	dn	ds	goc_score	wga_coverage	is_high_confidence	homology_id
# ENSG00000184895	ENSP00000372547	homo_sapiens	34.3137	ortholog_one2one	ENSMUSG00000069036	ENSMUSP00000088717	mus_musculus	17.7215	NULL	NULL	0	46.62	0	90358196
# ENSG00000114374	ENSP00000342812	homo_sapiens	82.818	ortholog_one2one	ENSMUSG00000069044	ENSMUSP00000088727	mus_musculus	82.7856	NULL	NULL	25	98.85	1	90352720
# ENSG00000242875	ENSP00000372484	homo_sapiens	31.6532	ortholog_many2many	ENSMUSG00000094511	ENSMUSP00000139784	mus_musculus	41.3158	NULL	NULL	0	0.00	0	90343821

# >> This is much cleaner and easier to parse into comparahomology table since it's a single file with a coherent format (whereas the individual files have 4 different formats with slightly different numbers of columns, column order etc.)
# >> However, it's lacking clade information, which I'm interested in.

cd("..");



# show directory
run("ls", "ls input");

done();
