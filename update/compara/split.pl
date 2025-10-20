#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [Compara release number]\n\nExample: $0 108";
($rel) = args(1);
# if (switch('tree')) { $mode .= "_tree"; }

$fastafile = "input/Compara.$rel.protein_default.aa.fasta.gz";
$cdsfile = "input/Compara.$rel.protein_default.cds.fasta.gz";
$treefile = "input/Compara.$rel.protein_default.nh.emf.gz";
$treexfile = "input/Compara.$rel.protein_default.nhx.emf.gz";
$outclustfile = "output-clusters.txt";
# $outclusthumanfile = "output-clusters-human-ensg.txt";

open(FASTA, "zcat $fastafile|") or die("\nError: Couldn't zcat '$fastafile'\n\n");
open(CDS, "zcat $cdsfile|") or die("\nError: Couldn't zcat '$cdsfile'\n\n");
open(TREES, "zcat $treefile|") or die("\nError: Couldn't zcat '$treefile'\n\n");
open(TREESX, "zcat $treexfile|") or die("\nError: Couldn't zcat '$treexfile'\n\n");
open(OUTCLUST, ">$outclustfile") or die("\nError: Couldn't open '$outclustfile'\n\n");
# open(OUTCLUSTHUMAN, ">$outclusthumanfile") or die("\nError: Couldn't open '$outclusthumanfile'\n\n");


# start

# Get ensembl sequences
# No need to get species since ENSPs are unique. Table 'comparaenspspecies' confirms this:
# SELECT COUNT(DISTINCT ensp), COUNT(DISTINCT species) FROM comparaenspspecies;
# SELECT ensp, COUNT(DISTINCT species) AS c FROM comparaenspspecies GROUP BY ensp ORDER BY c DESC, ensp;
$query = Query("SELECT ensp, seq FROM ensembl");
startme("Getting ENSP-Sequence mapping from table 'ensembl'", 0, Numrows($query));
%seq = ();
while (($ensp, $seq) = Fetch($query))
{
	$seq{$ensp} = $seq;
	stepme(100000);
}
stopme();


# # Hash of versioned ENSPVs (to be trimmed to ENSPs)
# # Currently only caenorhabditis_elegans has versioned ENSPVs (all other species have unversioned ENSPs)
# %ensp = ();

# # List of ENSPs to skip (because of sequence mismatches with table 'ensembl')
# %skip = ();

startme("Splitting '$fastafile' into one-alignment chunks");
starttime();
$/ = "\n//\n\n";
%f = ();
while ($chunk = <FASTA>)
{
	chomp($chunk);
	
	stepme(1);

	# replace spaces at the end of FASTA titles
	$chunk =~ s/ +$//mg;

	# print "\n\n\nCHUNK\n\n\n['$chunk']\n\n\n";


	$outfile = "tmp/".sprintf("%06d", getme()).".fasta.txt";
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

	# list FASTA titles
	$f{getme()} = '';
	# while ($chunk =~ /^>(.+)/mg)
	foreach $subchunk (split(/\n>/, $chunk))
	{
		# print "\n\n\nSUBCHUNK\n\n\n['$subchunk']\n\n\n";

		# $ensp = $1;
		($ensp, $seq) = getfasta($subchunk);

		# Remove gaps from sequence
		$gapseq = $seq;
		$seq =~ s/-//g;

		# print "\n\n\nENSP\n\n\n['$ensp']\n\n\n";
		# print "\n\n\nSEQ\n\n\n['$seq']\n\n\n";

		# Verify sequences with table 'ensembl'
		if (exists($seq{$ensp}))
		{
			if ($seq ne $seq{$ensp})
			{
				# warn("\n\nWARNING: Sequence mismatch between '$fastafile' and table 'ensembl' for ensp (skipped): $ensp\n\nENSEMBL $seq{$ensp}\n\nFOUND   $seq\n\n");
				# addme("sequence mismatch between '$fastafile' and table 'ensembl' for ensp (skipped)", $ensp);
				addme("sequence mismatch between '$fastafile' and table 'ensembl' for ensp (kept)", $ensp);

				# $skip{$ensp} = 1;
				# next;
			}
		}
		else
		{
			addme("no sequence in table 'ensembl' for verification for ensp (kept)", $ensp);
		}
		# Save in array (list of ENSPs in this chunk) for diagnostics (to verify that the ENSPs contained in a given chunk are identical between the FASTA file, the NH tree file, and the NHX tree file)
		$f{getme()} .= "$ensp|";

		# Print to single-alignment output file		
		# print "\n\n\nOUT\n\n\n>$ensp\n".split60($gapseq)."\n\n\n";
		print OUT ">$ensp\n".split60($gapseq)."\n";
	}
	# @enspvs = unique(@enspvs);

	# # Replace ENSPVs with ENSPs in this chunk's sequence titles
	# if (scalar(@enspvs) > 0)
	# {
	# 	# print "\n\n\nENSPVs\n".join("\n", @enspvs)."\n\n\n";
	# 	# print "\n\n\nPRECHUNK  $chunk\n\n\n";
	# 	foreach $enspv (@enspvs)
	# 	{
	# 		# Convert ENSPV to ENSP
	# 		$ensp = $enspv;
	# 		$ensp =~ s/\.\d+$//;
	# 		$chunk =~ s/^>$enspv$/>$ensp/mg;
	# 	}
	# 	# print "\n\n\nPOSTCHUNK $chunk\n\n\n";
	# 	# exit;
	# }

	# $outfile = "tmp/".sprintf("%06d", getme()).".fasta.txt";
	# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	# print OUT $chunk;
	close(OUT);
}
close(FASTA);
stopme();
stoptime();

startme("Splitting '$cdsfile' into one-alignment chunks");
starttime();
$/ = "\n//\n\n";
%c = ();
while (<CDS>)
{
	chomp;
	
	stepme(100);

	# remove spaces at the end of FASTA titles
	s/ +$//mg;

	# list FASTA titles
	$c{getme()} = '';
	while (/^>(.+)/mg)
	{
		$ensp = $1;

		# # Skip this ENSP if it is in the "skip list" (due to sequence mismatch with table 'ensembl')
		# next if (exists($skip{$ensp}));

		$c{getme()} .= "$ensp|";
	}

	$f{getme()} = join("|", unique(split(/\|/, $f{getme()})));
	$c{getme()} = join("|", unique(split(/\|/, $c{getme()})));
	
	if ($f{getme()} ne $c{getme()})
	{
		# die("Error: Proteins not the same in chunk ".getme().":\n\nFASTA\n".$f{getme()}."\nvs.\nTree\n".$c{getme()}."\n\n");
		warn("Warning: Proteins not the same in chunk ".getme().":\n\nFASTA\n".$f{getme()}."\nvs.\nTree\n".$c{getme()}."\n\n");
		addme("ERROR: proteins not the same for fasta and tree in chunk", getme());
	}
	
	# Print to clusters output file (all species, ENSP)
	print OUTCLUST $c{getme()}."\n";
	# # Print to clusters output file (human only, ENSG)
	# @ensgs = ();
	# foreach $ensp (split(/\|/, $c{getme()}))
	# {
	# 	if (exists($ensg{$ensp}))
	# 	{
	# 		push(@ensgs, $ensg{$ensp});
	# 	}
	# }
	# print OUTCLUSTHUMAN join("|", unique(@ensgs))."\n";
	

	$outfile = "tmp/".sprintf("%06d", getme()).".cds.txt";
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	print OUT $_;
	close(OUT);
}
close(CDS);
stopme();
stoptime();

startme("Splitting '$treefile' into one-alignment chunks");
starttime();
$/ = "//\n\n";
%t = ();
while (<TREES>)
{
	chomp;
	
	stepme(100);

	# list FASTA titles
	$t{getme()} = '';
	while (/^SEQ \S+ (\S+) /mg)
	{
		$ensp = $1;

		# # Skip this ENSP if it is in the "skip list" (due to sequence mismatch with table 'ensembl')
		# next if (exists($skip{$ensp}));

		$t{getme()} .= "$ensp|";
	}
	
	# Unique
	$f{getme()} = join("|", unique(split(/\|/, $f{getme()})));
	$t{getme()} = join("|", unique(split(/\|/, $t{getme()})));
	
	if ($f{getme()} ne $t{getme()})
	{
		# die("Error: Proteins not the same in chunk ".getme().":\n\nFASTA\n".$f{getme()}."\nvs.\nTree\n".$t{getme()}."\n\n");
		warn("Warning: Proteins not the same in chunk ".getme().":\n\nFASTA\n".$f{getme()}."\nvs.\nTree\n".$t{getme()}."\n\n");
		addme("ERROR: proteins not the same for fasta and tree in chunk", getme());
	}

	s/^SEQ.+?\n//msg;
	s/^DATA\n//msg;
	
	$outfile = "tmp/".sprintf("%06d", getme()).".tree.txt";
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	print OUT $_;
	close(OUT);
}
close(TREES);
stopme();
stoptime();

startme("Splitting '$treexfile' into one-alignment chunks");
starttime();
$/ = "//\n\n";
%t = ();
while (<TREESX>)
{
	chomp;
	
	stepme(100);

	# list FASTA titles
	$t{getme()} = '';
	while (/^SEQ \S+ (\S+) /mg)
	{
		$ensp = $1;

		# # Skip this ENSP if it is in the "skip list" (due to sequence mismatch with table 'ensembl')
		# next if (exists($skip{$ensp}));

		$t{getme()} .= "$ensp|";
	}
	
	# Unique
	$f{getme()} = join("|", unique(split(/\|/, $f{getme()})));
	$t{getme()} = join("|", unique(split(/\|/, $t{getme()})));
	
	if ($f{getme()} ne $t{getme()})
	{
		# die("Error: Proteins not the same in chunk ".getme().":\n\nFASTA\n".$f{getme()}."\nvs.\nTree\n".$t{getme()}."\n\n");
		warn("Warning: Proteins not the same in chunk ".getme().":\n\nFASTA\n".$f{getme()}."\nvs.\nTree\n".$t{getme()}."\n\n");
		addme("ERROR: proteins not the same for fasta and tree in chunk", getme());
	}

	s/^SEQ.+?\n//msg;
	s/^DATA\n//msg;
	
	$outfile = "tmp/".sprintf("%06d", getme()).".treex.txt";
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	print OUT $_;
	close(OUT);
}
close(TREES);
stopme();
stoptime();

state("Wrote homology clusters (from CDS file) to '$outclustfile'");

showmeall(1);

done();
