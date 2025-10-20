#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [table] [-debug]\n\n -debug: Don't actually make any changes.\n\n$0 -debug\n$0 ensembl -debug";
if (scalar(@ARGV) == 0)
{
	$ensembl = 'ensembl';
}
elsif (scalar(@ARGV) == 1)
{
	if ($ARGV[1] =~ /^-/)
	{
		# Only a switch (must be -debug)
		args(0);
	}
}
else
{
	($ensembl) = args(1);
	if (($ensembl ne 'ensembl') and !switch('debug'))
	{
		die("Error: Switch -debug must be active when specifying a nonstandard table (only 'ensembl' will actually be filtered)");
	}
}
print "Using source table: $ensembl\n";
# exit;

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

# # $tmphumanonlyand = '';
# $tmphumanonlywhere = '';
# if (switch('humanonly'))
# {
# 	# $tmphumanonlyand = " AND species='human'";
# 	$tmphumanonlywhere = " WHERE species='human'";
# }


# start

# # Show counts (uniens)
# showcounts();




# Simply clear the table now that I'm filling it from scratch (and using uniens_uniprot for getting species mnemonics in species.pl)
Clear('uniens');



# # Get unmapped accs and sequences from uniseq and add ENSP mappings to uniens based on sequence identity
# $mainquery = Query("SELECT DISTINCT s.name, s.acc, s.canon, s.species, s.seq FROM uniseq s LEFT OUTER JOIN uniens ue ON s.acc=ue.acc WHERE s.type IN ('UniProt', 'UniIso') AND ue.id IS NULL");
# Get all accs and sequences from uniseq and add ENSP mappings to uniens based on sequence identity (unless they already existed)
$mainquery = Query("SELECT DISTINCT name, acc, canon, species, seq FROM uniseq WHERE type IN ('UniProt', 'UniIso')");
startme("Sequence matching: Getting all acc/species/sequences from table 'uniseq', and adding ENSP mappings to table 'uniens' based on species/sequence identity in table '$ensembl'", 0, Numrows($mainquery));
starttime();
$inserted = 0;
while (($name, $acc, $canon, $species, $uniseq) = Fetch($mainquery))
{
	stepme(1000);

	$canonical = 0;
	$canonical = 1 if ($acc eq $canon);

	# addme("seqmatch: total names", $name);
	addme("seqmatch: total accs", $acc);
	# addme("seqmatch: total canons", $canon);
	# addme("seqmatch: total species", $species);
	# addme("seqmatch: total uniprot sequences", $uniseq);

	# Check if this UniProt sequence is in table 'ensembl' for this species
	$ensemblquery = Query("SELECT DISTINCT ensgv, ensg, enstv, enst, enspv, ensp FROM $ensembl WHERE species='$species' AND seq='$uniseq'");
	if (Numrows($ensemblquery) == 0)
	{
		# addme("seqmatch: no sequence match found in table $ensembl for name (skipped)", $name);
		addme("seqmatch: no sequence match found in table $ensembl for acc (skipped)", $acc);
		# addme("seqmatch: no sequence match found in table $ensembl for canon (skipped)", $canon);
		# addme("seqmatch: no sequence match found in table $ensembl for species (skipped)", $species);
		# addme("seqmatch: no sequence match found in table $ensembl for uniprot sequence (skipped)", $uniseq);
		next;
	}
	while (($ensgv, $ensg, $enstv, $enst, $enspv, $ensp) = Fetch($ensemblquery))
	{
		# # This was to avoid duplicate acc-ensp mappings, but skipping these means that for e.g. ECE2_MOUSE and EFCE2_MOUSE (which both have the acc B2RQR8-2) only one of them gets inserted, which is not useful.
		# # Since this table in any case has lots of duplicates (accs map to many ENSPs), I'll leave these "duplicates" in. They're still unique name-ENSP mappings.
		# # Queries:
		# # SELECT * FROM uniseq WHERE acc='B2RQR8-2';
		# # SELECT * FROM uniseq WHERE acc='B2RQR8-2' GROUP BY acc, seq;
		# # SELECT * FROM uniens WHERE acc='B2RQR8-2' AND ensp='ENSMUSP00000119693';
		# $query = Query("SELECT id FROM uniens WHERE acc='$acc' AND ensp='$ensp'");
		# if (Numrows($query) == 0)
		# {
			if (!switch('debug'))
			{
				Query("INSERT INTO uniens SET name='$name', acc='$acc', canon='$canon', canonical='$canonical', species='$species', ensgv='$ensgv', ensg='$ensg', enstv='$enstv', enst='$enst', enspv='$enspv', ensp='$ensp'");
			}
			addme("seqmatch: sequence match found for acc|ensp", "$acc|$ensp");
			addme("seqmatch: sequence match found for acc", $acc);

			# addme("seqmatch: sequence match successfully inserted for ensgv", $ensgv);
			# addme("seqmatch: sequence match successfully inserted for enstv", $enstv);
			# addme("seqmatch: sequence match successfully inserted for enspv", $enspv);
			# addme("seqmatch: sequence match successfully inserted for name", $name);
			addme("seqmatch: sequence match successfully inserted for acc|ensp", "$acc|$ensp");
			addme("seqmatch: sequence match successfully inserted for acc", $acc);
			# addme("seqmatch: sequence match successfully inserted for canon", $canon);
			# addme("seqmatch: sequence match successfully inserted for species", $species);
			# addme("seqmatch: sequence match successfully inserted for uniprot sequence", $uniseq);
			$inserted++;
	}
}
stopme();
stoptime();

state("Rows inserted: ".commify($inserted));

showmeall(1);
# showmeallsorted(1);

# Show counts (uniens)
showcounts();

Optimize('uniens');

done();





sub showcounts
{
	# Show counts (uniens)
	nl();
	state("Distinct counts in table 'uniens':", 1);
	$query = Query("SELECT COUNT(DISTINCT acc, ensp), COUNT(DISTINCT name), COUNT(DISTINCT acc), COUNT(DISTINCT canon), COUNT(DISTINCT ensp), COUNT(DISTINCT species) FROM uniens");
	($acc_enspcount, $namecount, $acccount, $canoncount, $enspcount, $speciescount) = Fetch($query);
	nl();
	state("All:   ".commify($acc_enspcount)." acc|ensps, ".commify($namecount)." names, ".commify($acccount)." accs, ".commify($canoncount)." canons, ".commify($enspcount)." ENSPs, ".commify($speciescount)." species", 1);
	$query = Query("SELECT COUNT(DISTINCT acc, ensp), COUNT(DISTINCT name), COUNT(DISTINCT acc), COUNT(DISTINCT canon), COUNT(DISTINCT ensp) FROM uniens WHERE species='HUMAN'");
	($acc_enspcount, $namecount, $acccount, $canoncount, $enspcount) = Fetch($query);
	state("Human: ".commify($acc_enspcount)." acc|ensps, ".commify($namecount)." names, ".commify($acccount)." accs, ".commify($canoncount)." canons, ".commify($enspcount)." ENSPs", 1);
	nl();
}
