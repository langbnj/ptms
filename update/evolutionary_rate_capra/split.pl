#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$rel = 108;		# Compara release (only used for mode 'para')


our $usage = "$0 [mode: e.g. linsi_tree, ginsi]";
($mode) = args(1);

$mode = lc($mode);
$mode =~ s/[^\w]/_/g;

if ($mode eq 'para')
{
	# para
	$fastafile = "../compara/input/Compara.$rel.protein_default.aa.fasta.gz";
}
else
{
	# linsi etc.
	$fastafile = "../comparanopara/output_$mode/Compara.protein_trees.aa.fasta.gz";
}

open(FASTA, "zcat $fastafile|") or die("\nError: Couldn't zcat '$fastafile'\n\n");


# start

# Not anymore! Removed all non-standard-20 AAs from comparafasta_linsi etc. sequences.
# state("Note: Also replaces selenocysteines (U) with cysteines (C) (otherwise Capra can't run, so this is still a better conservation estimate than nothing, and it's entirely at the protein sequence level)");

startme("Splitting '$fastafile' into one-alignment chunks");
starttime();
$/ = "\n//\n\n";
while (<FASTA>)
{
	chomp;

	# next if !/ENSP00000333956/;	# Selenocysteine example (U)
	
	stepme(100);

	# replace spaces at the end of FASTA titles
	s/ +$//mg;
	
	# print "['$_']\n";
	@a = split(/\n/);
	foreach (@a)
	{
			next if (/^>/);
	
			# Do selenocysteines (U) frequently have C in their place in other species?
			# Using linsi_tree as an example:
			# SELECT *, POSITION("U" IN seq) FROM comparafasta_linsi_tree WHERE seq REGEXP 'U';
			# SELECT DISTINCT aln, POSITION("U" IN seq) FROM comparafasta_linsi_tree WHERE seq REGEXP 'U';

			# SELECT c.aln, SUBSTRING(c.seq, 327, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=7520 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 201, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=8491 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 102, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=9301 ORDER BY c.aln, s.id; # >> Mostly C, but this is the only clear case like this
			# SELECT c.aln, SUBSTRING(c.seq, 461, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=9404 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 291, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=9416 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 643, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=9416 ORDER BY c.aln, s.id; # >> Mostly H
			# SELECT c.aln, SUBSTRING(c.seq, 1745, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=12032 ORDER BY c.aln, s.id; # >>  Mostly S
			# SELECT c.aln, SUBSTRING(c.seq, 292, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=13324 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 124, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=14381 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 278, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=14786 ORDER BY c.aln, s.id; # >> Mostly G
			# SELECT c.aln, SUBSTRING(c.seq, 279, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=14786 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 143, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=14792 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 372, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=15280 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 45, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=15307 ORDER BY c.aln, s.id;
			# SELECT c.aln, SUBSTRING(c.seq, 120, 1) FROM comparafasta_linsi_tree c, comparaspecies s WHERE c.species=s.species AND aln=16195 ORDER BY c.aln, s.id;
			# >> There is no clear pattern. However, using X would be a pity since the U is clearly extremely significant. I'll use C to maintain some signal.

			# Not anymore! Removed all non-standard-20 AAs from comparafasta_linsi etc. sequences.
			# # Replace in FASTA because Capra has trouble with these:
			# # # U (selenocysteine, encoded by a UGA stop codon) >> X (unknown)
			# # Positions with an X do not get scored (if the X occurs in the human reference).
			# # U (selenocysteine, encoded by a UGA stop codon) >> C (cysteine, UGU/UGC)
			# s/U/C/g;
			# # Stop codons (lead to KeyError) with X (unknown) (just to allow Capra to run - these sequences are obviously unreliable. They're also quite rare.)
			# # X and - get scored identically by Capra, but X seems more appropriate here.
			# s/\*/X/g;
	}
	$_ = join("\n", @a)."\n";
	# print "['$_']\n";	
	# exit;
	
	$outfile = "tmp/$mode.".sprintf("%06d", getme()).".fasta.txt";
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	print OUT $_;
	close(OUT);
}
close(FASTA);
stopme();
stoptime();

done();
