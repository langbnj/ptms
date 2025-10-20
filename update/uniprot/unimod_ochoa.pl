#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$ochoa_evid = 'ECO:0000312';		# "Another database"
$ochoa_pmid = '31819260';			# Ochoa et al. phosphorylation paper
$intable = 'unimod_ochoa_full';		# Input table (imported from Ochoa et al. paper by converting tables to CSV and then running my import.pl)

# our $superloudmysql = 1;

our $usage = "$0";
# ($var) = args(1);
args(0);



# start

nl();
state("Deleting existing source='Ochoa' sites from table 'unimod'...", 1);
if (!switch('debug'))
{
	$query = Query("DELETE FROM unimod WHERE source='Ochoa'");
}
else
{
	$query = Query("SELECT id FROM unimod WHERE source='Ochoa'");
}
state("Rows affected: ".commify(Numrows($query)), 1);


# Verified that uniprot|position is unique:
# SELECT uniprot, position, COUNT(DISTINCT residue) AS c FROM unimod_ochoa_full GROUP BY uniprot, position ORDER BY c DESC;
# SELECT uniprot, position, COUNT(*) AS c FROM unimod_ochoa_full GROUP BY uniprot, position ORDER BY c DESC;
# >> c is max 1
# >> No need to GROUP or do DISTINCT

$mainquery = Query("SELECT uniprot AS acc, position AS site, residue AS aa FROM $intable ORDER BY uniprot, position");

startme("Getting PTM sites from '$intable', verifying accessions and amino acids, and inserting into table 'unimod'", 0, Numrows($mainquery));
starttime();

while (($thisacc, $site, $aa) = Fetch($mainquery))
{
	# Manual workaround for an outdated acc where the current acc is also present in the input, resulting in 3 duplicate sites. This is the only case where this occurs. Since this is a static dataset, I'll simply skip the outdated acc.
	# Diagnostic queries:
	# SELECT COUNT(*), COUNT(DISTINCT acc, site), COUNT(DISTINCT name, site) FROM unimod WHERE source='Ochoa';
	# SELECT acc, site, COUNT(id) AS c FROM unimod WHERE source='Ochoa' GROUP BY acc, site ORDER BY c DESC;
	# SELECT * FROM unimod WHERE source='Ochoa' AND acc='Q2VIR3';
	# SELECT acc FROM uniacc WHERE canon='Q2VIR3';
	if ($thisacc eq 'F8W810')
	{
		addme("DEBUG: manually skipped outdated acc for Q2VIR3 (which is also present) for acc (skipped)", $thisacc);
		next;
	}

	# Set type and description
	if ($aa eq 'S')
	{
		# S-p
		$type = 'modified residue';
		$desc = 'Phosphoserine';
	}
	elsif ($aa eq 'T')
	{
		# T-p
		$type = 'modified residue';
		$desc = 'Phosphothreonine';
	}
	elsif ($aa eq 'Y')
	{
		# Y-p
		$type = 'modified residue';
		$desc = 'Phosphotyrosine';
	}

	# No isoforms:
	# SELECT * FROM unimod_ochoa_full WHERE uniprot LIKE '%-%';
	# >> 0 rows

	# Validate listed accession and update it to primary acc if it is outdated
	# Also get protein name and species
	# Non-isoform
	$namequery = Query("SELECT DISTINCT name, canon, canon AS acc, species, trembl FROM uniacc WHERE canon='$thisacc'");
	if (Numrows($namequery) == 0)
	{
		# If this acc isn't a canonical, primary acc:
		# Get updated acc
		# This is okay to do since uniseq only ever contains "canon" (latest, primary) accs anyway:
		# SELECT * FROM uniacc a LEFT OUTER JOIN uniseq s ON a.acc=s.acc AND s.type='UniProt' WHERE s.acc IS NOT NULL AND a.acc!=a.canon;
		# >> 0 rows
		$namequery = Query("SELECT DISTINCT name, canon, canon AS acc, species, trembl FROM uniacc WHERE acc='$thisacc'");
		if (Numrows($namequery) == 0)
		{
			# If this still didn't work:
			addme("listed non-isoform acc skipped (no match in uniacc) for listed acc (skipped)", $thisacc);
			next;
		}
	}

	# Get sequence (via canonical acc, since some are isoforms)
	while (($name, $canon, $acc, $species, $trembl) = Fetch($namequery))
	{
		# Verify site in uniseq (is it within the protein?)
		die("Error: Couldn't match site '$site'") if ($site !~ /^\d+$/);
		$query = Query("SELECT DISTINCT seq FROM uniseq WHERE acc='$acc' AND type IN ('UniProt', 'UniIso')");
		($uniseq) = FetchOne($query);
		if ($site >= length($uniseq))
		{
			addme("site outside of uniseq for name (skipped)", $name);
			addme("site outside of uniseq for description|name|site (skipped)", "$desc|$name|$site");
			next;
		}

		if (substr($uniseq, $site-1, 1) ne $aa)
		{
			addme("aa mismatch in uniseq for name (skipped)", $name);
			addme("aa mismatch in uniseq for description|name|site (skipped)", "$desc|$name|$site");
			next;
		}

		# Check UniProt name
		if ($name =~ /^(\w{1,10})_(\w{3,5})$/)
		{
			$species = $2;
			addme("total species", $species);
		}
		else
		{
			die("Error: Couldn't match name '$name'") ;
		}

		die("Error: Couldn't match site '$site'") if ($site !~ /^\d+$/);

		addme("total names inserted", $name);
		addme("total name|sites inserted", "$name|$site");
		addme("total accs inserted", $acc);
		addme("total acc|sites inserted", "$acc|$site");
		addme("total descriptions inserted", $desc);

		# TrEMBL / Swiss-Prot
		if ($trembl == 0)
		{
			addme("total swiss-prot accs inserted", $acc);
			addme("total swiss-prot acc|sites inserted", "$acc|$site");
		}
		else
		{
			addme("total trembl accs inserted", $acc);
			addme("total trembl acc|sites inserted", "$acc|$site");
		}

		# Insert into unimod
		$q = "INSERT INTO unimod SET name='$name', acc='$acc', canon='$canon', species='$species', site='$site', ptm='', type='$type', description='$desc', source='Ochoa', subset='', scale='', evid='$ochoa_evid', pmids='$ochoa_pmid'";
		$q =~ s/=''/=NULL/g;
		if (switch('debug'))
		{
			print "$q\n" unless switch('quiet');
		}
		else
		{
			Query($q);
		}

		addme("added PTM site for acc|site", "$acc|$site");

	}

	stepme(1000);
}
stopme();
stoptime();

# showmeall(1);
showmesome(10);

Optimize('unimod');

done();
