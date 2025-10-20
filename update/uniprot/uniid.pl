#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [-humanonly]\n\n -humanonly: Read 'human-only' file, instead of the mapping for all species.";
args(0);

if (!switch('humanonly'))
{
	$infile = "input/idmapping.dat.gz";
}
else
{
	$infile = "input/HUMAN_9606_idmapping.dat.gz";
}
open(IN, "zcat $infile|") or die("\nError: Couldn't zcat '$infile'\n\n");


# start
startme("Getting Swiss-Prot accessions from 'uniacc'");
starttime();
$query = Query("SELECT DISTINCT acc FROM uniacc");
%accs = ();
while (($acc) = Fetch($query))
{
	$accs{$acc} = 1;
	stepme(10000);
}
stopme();


startme("Reading ID mappings from '$infile' and inserting into table 'uniid'");
$skipped = 0;
while (<IN>)
{
	chomp;
	
	@a = split(/\t/);
	die if (scalar(@a) != 3);
	($acc, $type, $value) = @a;
	
	if (!exists($accs{$acc}))
	{
		# Accession isn't in uniacc (i.e. it's not in any of the XMLs I parse)
		$skipped++;
		next;
	}
	
	$query = Query("SELECT name, species FROM uniacc WHERE acc='$acc'");
	if (Numrows($query) == 0)
	{
		addme("uniacc: acc not found in table 'uniacc'", $acc);
		next;
	}
	while (($name, $species) = FetchOne($query))
	{
		if (switch('debug'))
		{
			print "INSERT INTO uniid SET acc='$acc', name='$name', species='$species', type='$type', value='$value'\n" unless (switch('quiet'));
		}
		else
		{
			Query("INSERT INTO uniid SET acc='$acc', name='".esc($name)."', species='".esc($species)."', type='".esc($type)."', value='".esc($value)."'");
		}

		addme("uniid: accs in uniid", $acc);
		addme("uniid: names in uniid", $name);
		addme("uniid: species in uniid", $species);
		addme("uniid: types in uniid", $type);
		
		stepme(10000);
	}
}
stopme();
stoptime();

state("Inserted ".commify(getme())." ID mappings ($skipped entries skipped because they weren't in uniacc)");
showmeall(1);
# state("For release 2012_01 (for example), these numbers should look like: \"Show only reviewed (534,242)  (UniProtKB/Swiss-Prot) or unreviewed (19,434,245)  (UniProtKB/TrEMBL) entries\"");
# state("For release 2013_10 (for example), these numbers should look like: \"Show only reviewed (541,561)  (UniProtKB/Swiss-Prot) or unreviewed (44,746,523)  (UniProtKB/TrEMBL) entries\"");
# state("For release 2013_11 (for example), these numbers should look like: \"Show only reviewed (541,762)  (UniProtKB/Swiss-Prot) or unreviewed (48,180,424)  (UniProtKB/TrEMBL) entries\"");
state("For release 2021_01 (for example), these numbers should look like (human): \"Show only reviewed (20,396)  (UniProtKB/Swiss-Prot) or unreviewed (174,126)  (UniProtKB/TrEMBL) entries\"");

done();
