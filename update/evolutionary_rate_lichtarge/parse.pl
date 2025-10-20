#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [path] [species] [type]\n\nExample: $0 output human linsi_tree";
($path, $species, $type) = args(3);

# start
if (switch('clear'))
{
    Clear("evorate_lichtarge_$type");
}

startme("Parsing '$species' '$type' output in '$path'");
startme2();
starttime();

foreach $outfile (`ls $path/*.ranks`)
{
	chomp($outfile);
	
	# Get ENSP
	$outfile =~ /\/(ENSP\d+)\.ranks$/ or next;
	$ensp = $1;
	
	# Check if ENSP already has evorates in 'evorate' for this type
	$query = Query("SELECT id FROM `evorate_lichtarge_$type` WHERE ensp='$ensp'");
	if (Numrows($query) > 0)
	{
		addme("ensp already in evorate table", $ensp);
		next;
	}
	
	# First: Parse only sequence from output file
	open(TMP, $outfile) or die("\nError: Couldn't open '$outfile'\n\n");
	$seq = '';
	while(<TMP>)
	{
		chomp;
		next if /^%/;
		next if /^$/;
	
		@a = split(/\s+/, $_);
		shift(@a);
		
		$pos = $a[1];
		$seq .= $a[2];
	
		# die("Error: Position mismatch (position $pos != seq length ".length($seq).", sequence is '$seq')") if ($pos != length($seq));
	}
	close(TMP);

	# # Compare output sequence to ENSEMBL sequence
	# $query = Query("SELECT sequence FROM ensembl WHERE ensp='$ensp'");
	# if (Numrows($query) == 0)
	# {
	# 	addme("no sequence in 'ensembl' for ensp", $ensp);
	# 	next;
	# }
	# ($enseq) = FetchOne($query);
	# if ($seq ne $enseq)
	# {
	# 	addme("'ensembl' sequence doesn't match output sequence for ensp", $ensp);
	# 	next;
	# }
	
	# Compare output sequence to UniProt sequence
	# $query = Query("SELECT name FROM uniseq WHERE seq='$seq' AND type='UniProt'");
	# $query = Query("SELECT name FROM uniseq WHERE seq='$seq' AND type IN ('UniProt', 'UniIso')");
	$query = Query("SELECT name FROM uniens WHERE ensp='$ensp'");
	if (Numrows($query) == 0)
	{
		# addme("no sequence match in 'uniseq' for output sequence for ensp", $ensp);
		addme("ensp not found in table 'uniens' for ensp (kept)", $ensp);
		# next;
	}
	
	# Second: Parse rates from output (same file)
	open(TMP, $outfile) or die("\nError: Couldn't open '$outfile'\n\n");
	while(<TMP>)
	{
		chomp;
		next if /^%/;
		next if /^$/;
	
		@a = split(/\s+/, $_);
		shift(@a);
		
		$pos = $a[1];
		$rate = $a[6];
		
		die("Error: Position > sequence length") if ($pos > length($seq));

		Query("INSERT INTO `evorate_lichtarge_$type` SET ensp='$ensp', site=$pos, rate=$rate");
		stepme2();
		addme("proteins inserted", $ensp);
	}
	close(TMP);
	stepme(100);
}
stopme();

print "Inserted ".getme2()." rates for ".getme()." proteins\n";

# showme("ensp already in evorate table", !switch('debug'));
# showme("proteins inserted", !switch('debug'));
# # showme("no sequence in 'ensembl' for ensp", 1);
# # showme("'ensembl' sequence doesn't match output sequence for ensp", 1);
# showme("no name in 'uniens' for ensp", !switch('debug'));
# showme("no sequence in 'uniseq' for name", !switch('debug'));
# showme("no sequence match in 'uniseq' for output sequence for ensp", !switch('debug'));
showmesome(20);

done();
