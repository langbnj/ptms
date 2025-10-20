#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# Minimum alignment size (skip smaller alns):
# $nseqmin = 10;
$nseqmin = 1;

our $usage = "$0 [alignment number] [mode: e.g. linsi_tree, ginsi] [window size]";
($aln, $mode, $window) = args(3);

$type = 'capra'.$window.'_'.$mode;
$table = 'evorate_'.$type;

# start

# # Get ENSPs that are already in $table and skip them below
# Too slow
# %skip_ensps = ();
# $query = Query("SELECT DISTINCT ensp FROM $table");
# startme("Getting ENSPs that are already in $table (to skip them below)", 0, Numrows($query));
# while (($ensp) = Fetch($query))
# {
# 	$skip_ensps{$ensp} = 1;
# 	stepme(1000);
# }
# stopme();




startme("Running Capra (type '$type') on alignment no. $aln from Compara");
starttime();
@ensps = ();

# Get human ENSP IDs from alignment
$infile = $mode.".".sprintf("%06d", $aln).".fasta.txt";
open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
%compara = ();
fastabreak();
while (<IN>)
{
	($ensp, $seq) = getfasta();
	
	# Remove dashes from sequence (it's an alignment right now)
	$seq =~ s/-//g;
	
	# if human sequence: store in hash
	if ($ensp =~ /^ENSP\d+$/)
	{
		die("Error: '$ensp' sequence already exists") if (exists($compara{$ensp}));
		# print " >> $ensp\n";
		addme("total human ensps", $ensp);
		
		# Skip if ENSP is already in $table
		# if (exists($skip_ensps{$ensp}))
		$query = Query("SELECT COUNT(*), COUNT(DISTINCT site), MAX(site) FROM $table WHERE ensp='$ensp'");
		# print "   >> SELECT COUNT(*), COUNT(DISTINCT site), MAX(site) FROM $table WHERE ensp='$ensp'\n";
		# Verify that the counts are consistent (and crash otherwise)
		($tmpcount, $tmpsites, $tmpmaxsite) = FetchOne($query);
		$tmpmaxsite = 0 if (!defined($tmpmaxsite));
		if (($tmpcount != $tmpsites) or ($tmpcount != $tmpmaxsite))
		{
			die("Error: Site counts in table '$table' differ for '$ensp' (should be ".length($seq).", but got $tmpcount $tmpsites $tmpmaxsite)")
		}
		# If there's already something in the table:
		# print "     >> tmpcount $tmpcount\n";
		if ($tmpcount > 0)
		{
			# Verify that the entire protein is there
			if ($tmpcount != length($seq))
			{
				die("Error: Site counts incomplete in table '$table' differ for '$ensp' (should be ".length($seq).", but got $tmpcount $tmpsites $tmpmaxsite)")
			}
			else
			{
				# print "       >> complete results for ensp already existed in $table for ensp (skipped)\n";
				addme("complete results for ensp already existed in $table for ensp (skipped)", $ensp);
				next;
			}
		}

		# Do not run Capra using human sequences that contain X as the reference (it will treat them as gaps, leading to a mismatch in length with the original sequence).
		# SELECT * FROM comparafasta_einsi WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# SELECT * FROM comparafasta_ginsi WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# SELECT * FROM comparafasta_linsi WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# SELECT * FROM comparafasta_einsi_tree WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# SELECT * FROM comparafasta_ginsi_tree WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# SELECT * FROM comparafasta_linsi_tree WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# SELECT * FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# >> Only affects comparafasta (para)

		# SELECT COUNT(DISTINCT aln), COUNT(DISTINCT ensp) FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# >> Only affects 54 alns and 62 ENSPs
		# >> alns affected: 173, 656, 745, 788, 961, 1583, 1808, 4215, 4429, 5184, 5227, 6301, 9879, 10255, 12074, 12122, 12303, 12345, 12634, 12816, 13366, 13501, 14199, 14635, 14980, 15613, 15892, 15917, 16705, 16763, 17213, 17393, 17400, 29905, 30579, 40072, 46621, 49347, 51331, 53031, 53415, 53781, 53803, 53859, 53897, 53907, 53945, 53947, 53950, 53969, 53971, 53973, 53974, 53982
		# SELECT * FROM comparafasta f, uniens ue WHERE f.seq REGEXP 'X' AND f.species='homo_sapiens' AND f.ensp=ue.ensp;
		# >> None of these X-containing ENSPs are mapped to any UniProt accessions, so they will drop out of my PTM analyses anyway.

		# Delete affected ENSPs (all human ENSPs from alns where a human sequence has an X):
		# SELECT DISTINCT aln FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens';
		# SELECT DISTINCT ensp FROM comparafasta WHERE species='homo_sapiens' AND aln IN (SELECT DISTINCT aln FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens');
		# SELECT * FROM evorate_capra0_para WHERE ensp IN (SELECT DISTINCT ensp FROM comparafasta WHERE species='homo_sapiens' AND aln IN (SELECT DISTINCT aln FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens'));
		# DELETE FROM evorate_capra0_para WHERE ensp IN (SELECT DISTINCT ensp FROM comparafasta WHERE species='homo_sapiens' AND aln IN (SELECT DISTINCT aln FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens'));
		# DELETE FROM evorate_capra1_para WHERE ensp IN (SELECT DISTINCT ensp FROM comparafasta WHERE species='homo_sapiens' AND aln IN (SELECT DISTINCT aln FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens'));
		# DELETE FROM evorate_capra3_para WHERE ensp IN (SELECT DISTINCT ensp FROM comparafasta WHERE species='homo_sapiens' AND aln IN (SELECT DISTINCT aln FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens'));
		# DELETE FROM evorate_capra5_para WHERE ensp IN (SELECT DISTINCT ensp FROM comparafasta WHERE species='homo_sapiens' AND aln IN (SELECT DISTINCT aln FROM comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens'));
		

		$compara{$ensp} = $seq;

		if ($seq =~ /X/)
		{
			addme("skipped human ENSP as reference sequence because its sequence contained X (would break Capra residue numbers - only occurs in comparafasta (para)) for ensp (skipped)", $ensp);
			next;
		}

		# Add ENSP to be used as Capra reference below
		push(@ensps, $ensp);
	}
}
normalbreak();
close(IN);

# Happens automatically below
# # Skip this alignment if it doesn't contain any human sequences (or if they're all already in the evorate table)
# if (scalar(@ensps) == 0)
# {
# 	state("No human sequences found in this alignment, exiting");
# 	exit;
# }

# Run Capra for each ENSP as the reference sequence + parse
foreach $ensp (@ensps)
{
	$aln = sprintf("%06d", $aln);
	
	@a = `grep $aln.fasta.txt ../output-check-$mode-$window.txt` or die("\n\nError: Command failed (something is wrong with check.pl's output, it's missing $aln.fasta.txt):\n\ngrep $aln.fasta.txt ../output-check-$mode-$window.txt\n\n");
	die if (scalar(@a) != 1);
	chomp($a[0]);
	@a = split(/\t/, $a[0]);
	$nseq = $a[1];
	# $lseq = $a[3];
	# $score = $nseq * $lseq;
	# $expect = $score * $factor;
	# print "     >> $nseq sequences * $lseq alignment length = $score\n";
	
	if ($nseq < $nseqmin)
	{
		print "       >> too few sequences (<$nseqmin), skipping!\n";
		next;
	}

# 	# This would skip any alns that don't contain a modified protein
# 	$mod = 0;
# 	$query = Query("SELECT acc FROM uniens WHERE ensp='$ensp' AND species='human'");
# 	if (Numrows($query) == 0)
# 	{
# 		addme("no accs in uniens for ensp", $ensp);
# 		next;
# 	}
# 	$query = Query("SELECT m.id FROM unimod m, uniens e WHERE e.ensp='$ensp' AND e.acc=m.acc AND e.species='human' AND m.ptm!=''");
# 	if (Numrows($query) > 0)
# 	{
# 		$mod = 1;
# 	}
# 
# 	if ($mod == 0)
# 	{
# 		print "     >> Doesn't have any PTMs according to unimod 'ptm' field\n";
# 		next;
# 	}
# 	else
# 	{
# 		print "     >> Has PTMs\n";
# 	}

	# Get Compara sequence
	if (!exists($compara{$ensp}))
	{
		die("Error: No sequence in Compara for ENSP '$ensp'");
		# addme("no sequence in compara for ensp", $ensp);
		# next;
	}
	
	# Try to get matching UniProt sequence
	# Not really necessary anymore, since the UniProt sequence will always match Ensembl's now. Also, this will fail for selenocysteine-containing sequences (in which U got replaced with C for Capra to run)
	# $query = Query("SELECT s.seq FROM uniseq s, uniens e WHERE e.ensp='$ensp' AND e.acc=s.acc AND e.species='human' AND s.type IN ('UniProt', 'UniIso') AND s.seq='$compara{$ensp}'");
	# if (Numrows($query) == 0)
	# {
	# 	addme("couldn't find a uniprot acc with matching sequence for ensp", $ensp);
	# 	warn("Warning: Couldn't find a uniprot acc with matching sequence for ENSP '$ensp':\n\nSEQ $compara{$ensp}");
	# 	next;
	# }
	$query = Query("SELECT acc FROM uniens WHERE ensp='$ensp' AND species='human'");
	if (Numrows($query) == 0)
	{
		addme("no acc found in uniens for human ensp (not mapped to any Uniprot protein, which is okay) (kept)", $ensp);
		# warn("Warning: Couldn't find acc via uniens for ENSP '$ensp'");
		# next;
	}
	else
	{
		addme("acc successfully found in uniens for human ensp (mapped to a Uniprot protein) (kept)", $ensp);
	}
	
	# print " >> $ensp\n";
	$outfile = "../output/$mode.".sprintf("%06d", $aln).".$ensp.capra$window.txt";

	#     -a	reference sequence. Print scores in reference to a specific sequence (ignoring gaps). Default prints the entire column. [sequence name]
	# 	>> Each human ENSP in the aln (one after another)
	#     -b	lambda for window heuristic linear combination. Default=.5 [real in [0,1]]
	# 	>> default
	#     -d	background distribution file, e.g., swissprot.distribution. Default=BLOSUM62 background [filename]
	# 	>> BLOSUM62
	#     -g	gap cutoff. Do not score columns that contain more than gap cutoff fraction gaps. Default=.3 [real in [0, 1)]
	# 	>> 0.999999 (allow any number of gaps in a column) (up to 999,999 in a million)
	#     -h	help. Print this message.
	# 
	#     -l	use sequence weighting. Default=True [True|False]
	# 	>> default (True). "Sequence weighting: An alignment will often contain sequences at a range of evolutionary distances. If an alignment consists of several very similar sequences, all columns may look conserved, and it will be difficult to discriminate positions under evolutionary pressure from those that are not. We implemented the sequence weighting method proposed in Henikoff and Henikoff (1994) that rewards sequences that are ‘surprising’. Sequence weighting is used with all methods and results given subsequently, except for R4S, which builds an evolutionary tree as the first part of its analysis."
	# 	>> I think leaving this on is a good idea. Different alns have different evolutionary distance coverage, so rewarding "surprisingly conserved" columns should be a good idea.
	#     -m	similarity matrix file, e.g., matrix/blosum62.bla or .qij. Default=identity matrix [filename]
	# 	>> BLOSUM62
	#     -n	normalize scores. Print the z-score (over the alignment) of each column raw score. Default=False
	# 	>> default (False). No normalisation.
	#     -o	name of output file. Default=output to screen [filename]
	# 	
	#     -p	use gap penalty. Lower the score of columns that contain gaps. Default=True [True|False]
	# 	>> default (True).
	#     -s	conservation estimation method. 
	# 		Options: shannon_entropy, property_entropy, property_relative_entropy, vn_entropy, relative_entropy, js_divergence, sum_of_pairs. Default=js_divergence
	# 	>> default (Jensen-Shannon divergence)
	#     -w	window size. Number of residues on either side included in the window. Default=3 [int]
	# 	>> 0, 1, 3 and 5 (to see which works best). 0 seems to be best for PTM sites.


	run("capra", "python2 ../conservation_code/score_conservation.py -g 0.999999 -w $window -m ../conservation_code/matrix/blosum62.bla -a $ensp -o $outfile $infile");
	
	addme("ran human ensps", $ensp);
	
	# Parse output
	open(TMP, $outfile) or die("\nError: Couldn't open '$outfile'\n\n");
	$seq = '';
	while(<TMP>)
	{
		chomp;
		next if /^#/;
		
		@a = split(/\t/, $_);
		
		$seq .= $a[1];
		$rate = $a[2];
		
		if ($rate ne '-1000')
		{
			Query("INSERT INTO $table SET ensp='$ensp', site=".length($seq).", rate=$rate");
		}
		else
		{
			addme("Warning: rate -1000 (should not occur with -g 0.999999) for ensp|site", $ensp.'|'.length($seq));
			addme("Warning: rate -1000 (should not occur with -g 0.999999) for ensp", $ensp);
		}
	}
	
	# Check for mismatch between the sequence in the output and the one from 'uniseq'
	#  (not anymore) ...but correct for the fact that leading "X"-es (any amino acid) are dropped in the output (not anymore)
	# $uniseq =~ s/^(X+)//;
	# if (defined($1) and (length($1) > 0))
	# {
	# 	# adjust positions according to 'ensembl' sequence (which may start with leading "X"-es)
	# 	Query("UPDATE `evorate` SET site=(site+".length($1).") WHERE name='$ensp' AND type='capra$window' ORDER BY site");
	# }
	# if ($seq ne $uniseq)
	# {
	# 	Query("DELETE FROM evorate WHERE name='$ensp' AND type='capra$window'");
	# 	die("Error: Sequence mismatch - deleted '$ensp' 'capra$window' results again from table 'evorate'\n\n'uniseq' for '$ensp': ['$uniseq']\n\nOutput: ['$seq']\n\n");
	# }
	
	stepme(1);
}
stopme();

# showme("ran ensps");
# print "\n";
# showme("no accs in uniens for ensp", 1);
# showme("couldn't find a uniprot acc with matching sequence for ensp", 1);
showmeall(1);

done();
