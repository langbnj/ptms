#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# our $superloudmysql = 1;

$table = "comparahomology";

# Keeps one paralog (the best-matching one by query_percent_id) for one2many and many2many cases
our $usage = "$0 [-hc]\n\n -hc: Only keep high-confidence orthologs\n\nExample: $0";
args(0);

$label = "1para";
if (switch('hc'))
{
	$label = "1parahc";
}


# start

run("Create directory", "mkdir -p input/$label", 1);

starttime();
if (switch('hc'))
{
	# hc=1 only
	$mainquery = Query("SELECT DISTINCT ensp1 FROM $table WHERE hc=1 AND ensp1 IS NOT NULL ORDER BY ensp1");
}
else
{
	# all
	$mainquery = Query("SELECT DISTINCT ensp1 FROM $table WHERE ensp1 IS NOT NULL ORDER BY ensp1");
}
startme("Getting orthologs as well as best-matching outparalogs from table '$table' for human proteins in other species, then getting sequences and writing protein and CDS FASTA files", 0, Numrows($mainquery));
while (($ensp) = Fetch($mainquery))
{
	# Get human (reference) protein sequence from comparafasta
	# $query = Query("SELECT REPLACE(seq, '-', '') FROM comparafasta WHERE ensp='$ensp' AND species='homo_sapiens'");
	# Get it from ensembl instead
	$query = Query("SELECT seq, cds FROM ensembl WHERE ensp='$ensp' AND fullspecies='homo_sapiens'");
	if (Numrows($query) == 0)
	{
		# addme("no sequence in comparafasta for ensp", $ensp);
		addme("no sequence in ensembl for ensp", $ensp);
		next;
	}
	# $seq = FetchOne($query);
	($seq, $cds) = FetchOne($query);

	# Print human (reference) protein sequence
	$outfile = "input/$label/$ensp.$label.txt";
	open(OUT, ">$outfile") or die("Error: Couldn't open '$outfile'");
	fastabreak();
	print OUT ">$ensp\n".split60($seq)."\n";

	# Print human (reference) CDS sequence
	$cdsfile = "input/$label/$ensp.$label.cds.txt";
	open(CDS, ">$cdsfile") or die("Error: Couldn't open '$cdsfile'");
	fastabreak();
	print CDS ">$ensp\n".split60($cds)."\n";
	
	# Write homolog sequences (rest of the homology cluster)
	# Homology types documentation: https://useast.ensembl.org/info/genome/compara/homology_types.html
	# Using sorting elaborated in ~/Documents/MySQL/comparahomology.sql
	if (switch('hc'))
	{
		# hc=1 only
		# Order by ENSP to get the same order as before
		$speciesquery = Query("SELECT DISTINCT fullspecies2 FROM comparahomology WHERE ensp1='$ensp' AND hc=1 ORDER BY ensp2");
	}
	else
	{
		# all
		# Order by ENSP to get the same order as before
		$speciesquery = Query("SELECT DISTINCT fullspecies2 FROM comparahomology WHERE ensp1='$ensp' ORDER BY ensp2");
	}
	while (($fullspecies2) = Fetch($speciesquery))
	{
		if (switch('hc'))
		{
			# hc=1 only
			$homoquery = Query("SELECT ensp2 FROM comparahomology WHERE ensp1='$ensp' AND ensp2 IS NOT NULL AND fullspecies2='$fullspecies2' AND hc=1 ORDER BY (gene_order_score + whole_genome_align_score) / 2 DESC, gene_order_score DESC, whole_genome_align_score DESC, query_percent_id DESC, target_percent_id DESC LIMIT 1");
		}
		else
		{
			# all (this is what was used for einsi_tree_1para)
			$homoquery = Query("SELECT ensp2 FROM comparahomology WHERE ensp1='$ensp' AND ensp2 IS NOT NULL AND fullspecies2='$fullspecies2' ORDER BY hc DESC, (gene_order_score + whole_genome_align_score) / 2 DESC, gene_order_score DESC, whole_genome_align_score DESC, query_percent_id DESC, target_percent_id DESC LIMIT 1");
		}
		($homoensp) = FetchOne($homoquery);

		# Protein

		# From comparafasta
		# $query = Query("SELECT REPLACE(seq, '-', '') FROM comparafasta WHERE ensp='$homoensp'");
		# From ensembl
		$query = Query("SELECT seq, cds FROM ensembl WHERE ensp='$homoensp'");
		if (Numrows($query) == 0)
		{
			# addme("no sequence in comparafasta for homoensp", $homoensp);
			addme("no sequence in ensembl for homoensp", $homoensp);
			next;
		}
		# $seq = FetchOne($query);
		($seq, $cds) = FetchOne($query);

		print OUT ">$homoensp\n".split60($seq)."\n";

		# # CDS
		# 
		# $query = Query("SELECT REPLACE(seq, '-', '') FROM comparacds WHERE ensp='$homoensp'");
		# if (Numrows($query) == 0)
		# {
		# 	addme("no sequence in comparacds (but in comparafasta!) for homoensp", $homoensp);
		# 	next;
		# }
		# $cds = FetchOne($query);

		print CDS ">$homoensp\n".split60($cds)."\n";
	}

	close(OUT);
	close(CDS);

	# stepme(100);
	stepme(1);
}
stopme();
stoptime();

# showmesome();
showmeall();

done();
