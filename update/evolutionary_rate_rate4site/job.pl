#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
# our $usage = "$0 [species] [mode: 'fast', 'slow'] [alignment number] [ensp id as reference sequence] [-clear]\n\n -clear: Overwrite values in table 'evorate_rate4site_type' for the given ENSP if they already exist.\n\nExample: $0 human linsi_tree fast 1 ENSP0000001";
our $usage = "$0 [species] [type] [alignment number] [ensp id as reference sequence] [-clear]\n\n -clear: Overwrite values in table 'evorate_rate4site_type' for the given ENSP if they already exist.\n\nExample: $0 human linsi_tree 1 ENSP0000001";
($species, $type, $n, $ensp) = args(4);
# die("Error: mode must be 'fast' or 'slow'") if (($mode ne 'fast') and ($mode ne 'slow'));

# start

# startme("Running rate4site on alignment no. $n from Compara");
state("Running rate4site on alignment no. $n, ENSP ID $ensp from Compara");
starttime();

@ensps = ();

# Get Ensembl species 
$query = Query("SELECT species FROM comparaspecies WHERE unispec='$species'");
($comparaspec) = FetchOne($query);
state("UniProt species '$species' >> Ensembl species '$comparaspec'");

# Get $species ENSP IDs from alignment
$fastafile = sprintf("%06d", $n).".fasta.txt";
open(FASTA, "$fastafile") or die("\nError: Couldn't open '$fastafile'\n\n");
while (<FASTA>)
{
	if (/^>(\S+)$/)
	{
		# push(@ensps, $1);
		
		$title = $1;
		$query = Query("SELECT species FROM comparaenspspecies WHERE ensp='$title'");
		($tmpspec) = FetchOne($query);
		# If this is the desired species..
		if ($tmpspec eq $comparaspec)
		{
			push(@ensps, $title);
		}
	}
}
close(FASTA);

# Skip this alignment if it doesn't contain any $species sequences
if (scalar(@ensps) != scalar(unique(@ensps)))
{
	die("Error: Duplicate sequences in '$fastafile'!");
}

# Skip this alignment if it doesn't contain any $species sequences
if (scalar(@ensps) == 0)
{
	print " >> Error: Exiting: No $species sequences in this alignment!\n";
	exit;
}

# Run rate4site for each ENSP as the reference sequence + parse
# foreach $ensp (@ensps)
# {
	$success = 0;

	print " >> $ensp\n";
	$treefile = sprintf("%06d", $n).".tree.txt";
	$outfile = "../output/$ensp.txt";
	$treeoutfile = "tmp-$ensp-treeout.txt";
	$normalizedoutfile = "tmp-$ensp-normalizedout.txt";
	$logfile = "tmp-$ensp-log.txt";
	
	if (switch('clear'))
	{
		# $query = Query("DELETE FROM `evorate` WHERE ensp='$ensp' AND type='rate4site_$mode'");
		$query = Query("DELETE FROM `evorate_rate4site_$type` WHERE ensp='$ensp'");
		state("Cleared old 'evorate_rate4site_$type' values for ensp '$ensp'");
	}

	# $query = Query("SELECT id FROM `evorate` WHERE ensp='$ensp' AND type='rate4site_$mode'");
	$query = Query("SELECT id FROM `evorate_rate4site_$type` WHERE ensp='$ensp'");
	if (Numrows($query) > 0)
	{
		# print "     >> Error: Already have 'rate4site_$mode' values for '$ensp' in table 'evorate'! Skipping (use -clear switch on job.pl or on submit.pl to clear old values)\n";
		print "     >> Error: Already have values for '$ensp' in table 'evorate_rate4site_$type'! Skipping (use -clear switch on job.pl or on submit.pl to clear old values)\n";
		exit;
	}
	
	print "     >> Running rate4site...\n";
	
	starttime();
	# if ($mode eq 'fast')
	# {
        # system("../rate4site_sources_3.3_fast/programs/rate4site/rate4site -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $outfile >& $logfile");
	# }
	# elsif ($mode eq 'slow')
	# {
        system("../rate4site_sources_3.3_slow/programs/rate4site/rate4site.doubleRep -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $outfile >& $logfile");
	# }
	
	print "       >> ";
	stoptime(1);
	
	if (`grep -c "likelihoodComputation::getLofPos: likelihood of pos was zero!" $logfile` > 0)
	{
		print "         >> Error: rate4site crashed! ('likelihood of pos was zero')\n";
		# next;
		exit;
	}
	
	if (`grep -c "Amino acid was not one of the following" $logfile` > 0)
	{
		print "         >> Error: rate4site crashed! ('Amino acid was not one of the following')\n";
		# next;
		exit;
	}
	
	# Parse output
	open(TMP, $outfile) or die("\nError: Couldn't open '$outfile'\n\n");
	$seq = '';
	while (<TMP>)
	{
		#       1     M  0.1426   [7.164e-322,0.2178]  0.1727    3/48
		# 	    2     S  0.6591   [0.3887,0.8783]  0.2995   45/48
		# 	    3     N   0.966   [0.4888, 1.264]  0.4619   45/48

		chomp;
		next if /^#/;
		next if /^$/;
		# Grab bit before [
		/^(.+)\[/ or die("Error: Couldn't match line '$_'");
		$s = $1;
		# Strip whitespace
		$s =~ s/^\s+//;
		$s =~ s/\s+$//;
		
		# Split
		@a = split(/\s+/, $s);
		# shift(@a);
		
		# print "$_ >> 0='".$a[0]."', 1='".$a[1]."', 1='".$a[2]."'\n";
		
		$pos = $a[0];
		$aa = $a[1];
		$rate = $a[2];

		$seq .= $aa;
		
		die("Error: Couldn't match position '$pos' in line '$_'") if ($pos !~ /^\d+$/);
		die("Error: Couldn't match amino acid '$aa' in line '$_'") if (!aa($aa));
        # die("Error: Couldn't match evorate '$rate' in line '$_'") if ($rate !~ /^[\d\.]+$/);
		die("Error: Couldn't match evorate '$rate' in line '$_'") if ($rate !~ /^[\d\.e\-]+$/);
		
		die("Error: Position mismatch") if ($pos != length($seq));
		
		# Insert into evorate
		# Query("INSERT INTO `evorate` SET ensp='$ensp', type='rate4site_$mode', site=$pos, rate=$rate");
		Query("INSERT INTO `evorate_rate4site_$type` SET ensp='$ensp', site='$pos', rate='$rate'");
		$success = 1;
	}
	
	if ($success == 1)
	{
		print "       >> Success!\n";
	}
	else
	{
		print "       >> Error: Failed! (No rates inserted into evorate.)\n";
	}
	
	# Check for mismatch between the sequence in the output and the one from 'uniseq'
	#  (not anymore) ...but correct for the fact that leading "X"-es (any amino acid) are dropped in the output (not anymore)
	# $uniseq =~ s/^(X+)//;
	# if (defined($1) and (length($1) > 0))
	# {
	# 	# adjust positions according to 'ensembl' sequence (which may start with leading "X"-es)
	# 	Query("UPDATE `evorate` SET site=(site+".length($1).") WHERE ensp='$ensp' AND type='Capra$window' ORDER BY site");
	# }
	# if ($seq ne $uniseq)
	# {
	# 	Query("DELETE FROM evorate WHERE ensp='$ensp' AND type='rate4site_$mode'");
	# 	die("Error: Sequence mismatch - deleted '$ensp' 'rate4site_$mode' results again from table 'evorate'\n\n'uniseq' for '$ensp': ['$uniseq']\n\nOutput: ['$seq']\n\n");
	# }
	
	# stepme(1);
# }
# stopme();

done();
