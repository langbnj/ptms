#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# Manual mapping (not used):
# I could try to fix these species here manually, but it's not worth it since they have so few PTM sites (48):
# SELECT * FROM unitax WHERE species IN ('NOTSC', 'NOTSN', 'OVIMU', 'SHEEP');
# 3485	70142	NOTSC
# 5384	1027870	NOTSN
# 4818	9938	OVIMU
# 261	9940	SHEEP
# SELECT * FROM unimod WHERE species IN ('NOTSC', 'NOTSN', 'OVIMU', 'SHEEP') AND ptm IS NOT NULL;
# >> 48 sites

# $manual_mapping{""} = "";
# $manual_mapping{""} = "";
# $manual_mapping{""} = "";
# $manual_mapping{""} = "";
# # $manual_mapping{""} = "";
# 
# $manual_tax_fix{} = 70142;
# $manual_tax_fix{} = ;
# $manual_tax_fix{} = ;
# $manual_tax_fix{} = ;
# # $manual_tax_fix{} = ;

our $usage = "$0 [-debug] [-noensembl]\n -debug: Diagnostic only, don't touch any of the tables.\n -noensembl: Don't update the ensembl table (only ensembl_species) (fast, for debugging).";
args(0);



# start

if (!switch('debug'))
{
	Query("UPDATE ensembl SET species=NULL") if (!switch('noensembl'));
	Query("UPDATE ensembl_species SET species=NULL");
	state("Cleared 'species' field in tables 'ensembl' and 'ensembl_species'");
}

state("Filling 'species' field in tables 'ensembl' and 'ensembl_species' via tables 'unitax' and 'uniens_uniprot'...");
starttime();

$mainquery = Query("SELECT DISTINCT fullspecies, tax FROM ensembl_species ORDER BY fullspecies, tax");
%tax_to_unispec = ();
%unispec_to_fullspecies = ();
$affected_ensembl = 0;
$affected_ensembl_species = 0;
while (($fullspecies, $tax) = Fetch($mainquery))
{
    addme("total ensembl fullspecies", $fullspecies);

	# e.g. Acanthochromis_polyacanthus
    
   	# Get UniProt species, if available
   	$species = '';

    # via NCBI taxon ID
    print " >> $fullspecies";

	# # Workarounds here:
	# # Ensembl's yeast is S288c, so YEAST, not YEASX - UniProt uses the taxon ID 559292 though. Ensembl uses the non-specific 4932, so YEASX.
	# In Ensembl 108 (and probably earlier), this has been fixed. Ensembl now also uses the S288c-specific taxon ID 559292, so Ensembl and UniProt match.
	# if ($tax eq '4932') { $tax = '559292'; }
    print " (taxon $tax)\n";


	# Skip this fullspecies/tax if one of the fullspecies with this taxon ID already has a unispec attached.
	# This is to avoid e.g. assigning the same unispec to multiple mouse strains. Only the first strain will receive a unispec mapping.
	# The canonical strain/breed is always the first listed, and it's the one that is in Compara. (This might not be the case in future releases, but it is in Ensembl 108.)
	if (exists($tax_to_unispec{$tax}))
	{
		# addme("already assigned unispec to a different tax for unispec (skipped)", $species) if ($species ne '');
		addme("already assigned unispec to a different tax for fullspecies (skipped)", $fullspecies);
		addme("already assigned unispec to a different tax for tax (skipped)", $tax);
		next;
	}
	
	

	# Get UniProt species mnemonic (unispec) via unitax
    $query = Query("SELECT DISTINCT species FROM unitax WHERE tax='$tax'");
    # $query = Query("SELECT DISTINCT species FROM unitax_from_mapping_file WHERE tax='$tax'");
    if (Numrows($query) == 0)
    {
		# die("\n\nError: No unispec found in unitax for fullspecies '$fullspecies' (tax '$tax')\n\n");
        addme("no unispec found via unitax for fullspecies (kept)", $fullspecies);
        addme("no unispec found via unitax for tax (kept)", $tax);
	    print "   >> unitax >> -\n";
	    # print "\n";
        # next;
		$unitax_species = '';
    }
	else
	{
	    ($unitax_species) = FetchOne($query);
	    print "   >> unitax >> $unitax_species\n";
	    addme("uniprot species found via unitax for fullspecies", $fullspecies);
	    addme("uniprot species found via unitax for uniprot species", $unitax_species);
	}


	# Also get UniProt species mnemonic (unispec) via:
	# ensembl's fullspecies -> ensembl's ensg -> uniens_uniprot's ensg -> uniens_uniprot's unispec
	# Note: This needs to be ensg, not ensp or even enspv, since uniens_uniprot is for an outdated Ensembl release (perhaps 105). I have Ensembl 108 here.
	$query = Query("SELECT DISTINCT ue.species FROM ensembl e, uniens_uniprot ue WHERE e.fullspecies='$fullspecies' AND e.ensg=ue.ensg");
    if (Numrows($query) == 0)
    {
		# die("\n\nError: No unispec found in unitax for fullspecies '$fullspecies' (tax '$tax')\n\n");
        addme("no unispec found via uniens_uniprot for fullspecies (kept)", $fullspecies);
        addme("no unispec found via uniens_uniprot for tax (kept)", $tax);
	    print "   >> uniens_uniprot >> -\n";
	    # print "\n";
        # next;
		$uniens_uniprot_species = '';
    }
	else
	{
	    ($uniens_uniprot_species) = FetchOne($query);
	    print "   >> uniens_uniprot >> $uniens_uniprot_species\n";
	    addme("uniprot species found via uniens_uniprot for fullspecies", $fullspecies);
	    addme("uniprot species found via uniens_uniprot for uniprot species", $uniens_uniprot_species);
	}
	

	# Check if uniprot species obtained via unitax and uniens_uniprot match (if both are available)
	# If both defined:
	if (($unitax_species ne '') and ($uniens_uniprot_species ne ''))
	{
		if ($unitax_species eq $uniens_uniprot_species)
		{
			addme("uniens_uniprot and unitax species both found and they match for fullspecies", $fullspecies);
			addme("uniens_uniprot and unitax species both found and they match for tax", $tax);
			addme("uniens_uniprot and unitax species both found and they match for uniprot species", $unitax_species);
			$species = $unitax_species;
		}
		else
		{
			warn("Error/WARNING: uniens_uniprot and unitax species MISMATCH for fullspecies '$fullspecies', unitax '$unitax_species', and uniens_uniprot '$uniens_uniprot_species' (skipped)");
			addme("WARNING: uniens_uniprot and unitax species both found and they MISMATCH for fullspecies (skipped)", $fullspecies);
			addme("WARNING: uniens_uniprot and unitax species both found and they MISMATCH for tax (skipped)", $tax);
			addme("WARNING: uniens_uniprot and unitax species both found and they MISMATCH for unitax|uniens_uniprot species (skipped)", "$unitax_species|$uniens_uniprot_species");
			next;
		}
	}
	elsif ($unitax_species ne '')
	{
		# Only unitax is defined
		addme("only unitax species found (not in uniens_uniprot) for fullspecies", $fullspecies);
		addme("only unitax species found (not in uniens_uniprot) for tax", $tax);
		$species = $unitax_species;
	}
	elsif ($uniens_uniprot_species ne '')
	{
		# Only uniens_uniprot is defined
		addme("only uniens_uniprot species found (not in unitax) for fullspecies", $fullspecies);
		addme("only uniens_uniprot species found (not in unitax) for tax", $tax);
		$species = $uniens_uniprot_species;
	}
	else
	{
		# # Try manual mapping
		# if exists($manual_mapping{$fullspecies})
		# {
		# 	$species = $manual_mapping{$fullspecies};
		# 	addme("manually assigned species for fullspecies (kept)", $fullspecies);
		# 	addme("manually assigned species for tax (kept)", $tax);
		# }
		# else
		# {
			# Try to get UniProt species mnemonic (unispec) via:
			# ensembl's fullspecies -> tax's taxon name -> tax's taxon ID -> unitax's unispec
			# Note: This needs to be ensg, not ensp or even enspv, since uniens_uniprot is for an outdated Ensembl release (perhaps 105). I have Ensembl 108 here.
			$tmp_fullspecies = $fullspecies;
			$tmp_fullspecies =~ s/_/ /g;
			$query = Query("SELECT DISTINCT ut.species FROM unitax ut, tax t WHERE t.name='$tmp_fullspecies' AND t.tax=ut.tax");
			if (Numrows($query) == 0)
			{
				# die("\n\nError: No unispec found in unitax for fullspecies '$fullspecies' (tax '$tax')\n\n");
				addme("no unispec found via unitax nor uniens_uniprot nor tax for fullspecies (skipped)", $fullspecies);
				addme("no unispec found via unitax nor uniens_uniprot nor tax for tax (skipped)", $tax);
				addme("no unispec found via tax for fullspecies (skipped)", $fullspecies);
				addme("no unispec found via tax for tax (skipped)", $tax);
				print "   >> tax >> -\n";
				# print "\n";
				next;
				$tax_species = '';
			}
			else
			{
				($tax_species) = FetchOne($query);
				print "   >> tax >> $uniens_uniprot_species\n";
				addme("no unispec found via unitax nor uniens_uniprot, but found one via tax for fullspecies (kept)", $fullspecies);
				addme("no unispec found via unitax nor uniens_uniprot, but found one via tax for tax (kept)", $tax);
				addme("uniprot species found via tax for fullspecies (kept)", $fullspecies);
				addme("uniprot species found via tax for uniprot species (kept)", $uniens_uniprot_species);
			}
		# }
	}
	
    addme("uniprot species successfully found for fullspecies", $fullspecies) if ($species ne '');
    addme("uniprot species successfully found for uniprot species", $species) if ($species ne '');

	$tax_to_unispec{$tax} = $species;

	# Skip this mapping if this unispec has already been attached to a fullspecies.
	# This is to avoid e.g. assigning the same unispec to multiple mouse strains. Only the first strain will receive a unispec mapping.
	# The canonical strain/breed is always the first listed, and it's the one that is in Compara. (This might not be the case in future releases, but it is in Ensembl 108.)
	if (exists($unispec_to_fullspecies{$species}))
	{
		addme("already assigned unispec to a different fullspecies for unispec (skipped)", $species) if ($species ne '');
		addme("already assigned unispec to a different fullspecies for fullspecies (skipped)", $fullspecies);
		addme("already assigned unispec to a different fullspecies for tax (skipped)", $tax);
		next;
	}
	$unispec_to_fullspecies{$species} = $fullspecies;
	


	if (!switch('noensembl'))
	{
		if (!switch('debug'))
		{
		    $query = Query("UPDATE ensembl SET species='$species' WHERE fullspecies='$fullspecies'");
		}
		else
		{
		    $query = Query("SELECT id FROM ensembl WHERE fullspecies='$fullspecies'");
		}
		die("Error: No matches for update query:\n\nUPDATE ensembl SET species='$species' WHERE fullspecies='$fullspecies'\n\n") if (Numrows($query) == 0);
	    $affected_ensembl += Numrows($query);
	}

	if (!switch('debug'))
	{
		$query = Query("UPDATE ensembl_species SET species='$species' WHERE fullspecies='$fullspecies' AND tax='$tax'");
	}
	else
	{
		$query = Query("SELECT id FROM ensembl_species WHERE fullspecies='$fullspecies' AND tax='$tax'");
	}
	die if (Numrows($query) != 1);
    $affected_ensembl_species += Numrows($query);
}
nl();
nl();
stoptime();

starttime();
Optimize('ensembl_species');
Optimize('ensembl');
stoptime();

# showme("no unispec found via uniens_uniprot for fullspecies (kept)");
# showme("no unispec found via unitax or uniens_uniprot for fullspecies (skipped)");
# showme("only unitax species found (not in uniens_uniprot) for fullspecies");

showmeall(1);
# showmesome(5);
# showmesomesorted(5);

# Show rows affected
nl();
$query = Query("SELECT id FROM ensembl");
$all = Numrows($query);
state("Table 'ensembl': ".commify($affected_ensembl)." of ".commify($all)." rows affected", 1);
$query = Query("SELECT id FROM ensembl_species");
$all = Numrows($query);
state("Table 'ensembl_species': ".commify($affected_ensembl_species)." of ".commify($all)." rows affected", 1);
nl();


# Get fullspecies that are mapped to multiple unispecs and show them (shouldn't happen)
$query = Query("SELECT fullspecies, COUNT(DISTINCT species) AS c, GROUP_CONCAT(DISTINCT species ORDER BY species) FROM ensembl_species GROUP BY fullspecies HAVING c > 1 ORDER BY fullspecies");
if (Numrows($query) > 0)
{
	warn("WARNING: There are fullspecies that are mapped to multiple unispecs in table 'ensembl_species' (shouldn't happen)");
	state("Getting fullspecies that are mapped to multiple unispecs:");
	while (($fullspecies, $unispeccount, $unispecs) = Fetch($query))
	{
		print " >> $fullspecies >> $unispeccount unispecs >> $unispecs\n";
	}
}


# It's normal that not all Ensembl species have a UniProt mnemonic yet (there are dozens of cases)
# $query = Query("SELECT DISTINCT fullspecies FROM ensembl_species WHERE species IS NULL");
# # while (($fullspecies) = Fetch($query))
# # {
# # 	warn("Error: '$fullspecies' still doesn't have a UniProt species assigned to it");
# # }
# $notyet = Numrows($query);
# if ($notyet < $all)
# {
# 	warn("Warning: $notyet out of $all species still don't have a UniProt mnemonic assigned to them. This is probably because they're not in Swiss-Prot yet.");
# }

done();
