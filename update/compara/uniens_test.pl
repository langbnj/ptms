#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [Compara release number]\n\nExample: $0 108";
($rel) = args(1);

$infile = "input/Compara.$rel.protein_default.aa.fasta.gz";
open(IN, "zcat $infile|") or die("\nError: Couldn't zcat '$infile'\n\n");


# start

startme("Reading human ENSP IDs from '$infile'");
starttime();

@ensp = ();
while (<IN>)
{
	chomp;
	
	next if !/^>/;
	
	if (/^>(ENSP\d+)/)
	{
		# human protein
		push(@ensp, $1);

		stepme(1000);
	}
}
stopme();
stoptime();

@ensp = unique(@ensp);

state("Read ".getme()." human ENSP IDs (".scalar(@ensp)." unique)");

startme("Verifying that these human ENSPs exist in table 'uniens'");
@found = ();
foreach $ensp (@ensp)
{
	$query = Query("SELECT acc FROM uniens WHERE ensp='$ensp' AND species='human'");
	if (Numrows($query) > 0)
	{
		push(@found, $ensp);
	}
	stepme(1000);
}
stopme();

state("Found ".scalar(@found)." / ".scalar(@ensp)." human ENSP IDs in table 'uniens' (that's ".sprintf("%.1f", (scalar(@found) / scalar(@ensp)) * 100)."%)");

done();
