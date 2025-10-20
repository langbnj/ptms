#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$speciestable = 'comparaspecies';
$enspspeciestable = 'comparaenspspecies';
$ensemblspeciestable = 'ensembl_species';

# Previously: just topology (no branch lengths)
# $speciestreefile = "input/species_tree.ensembl.topology.nw";
# $speciestreefile = "input/species_tree.branch_len.nw"; # Has molecular branch lengths (see download.pl for more details on these trees)
# $speciestreefile = "input/default_protein-trees_cafe.nh";	# This tree has exactly 200 species, and fairly correct species names, but TimeTree branch lengths
# $speciestreefile = "input/vertebrates_species-tree_Ensembl.nh";
# $speciestreefile = "input/species_tree.branch_len.reordered_in_dendroscope.nw";
# Ensembl Compara 108:
# Reordered so human is at the top:
# species_tree.branch_len.reordered_in_dendroscope.nw
# Reordering steps:
# - Click "Homo sapiens"
# - Click the "reroot" button in the top toolbar
# - Click the "ladderize right" button in the top toolbar
# - Flip additional branches to put Ciona, Drosophila and yeast (all the "other" species in table "comparaspecies") at the bottom
# - Click File > Export to Newick format (saving as .nw rather than .tree)
# $speciestreefile = "input/species_tree.branch_len.reordered_in_dendroscope.nw";
$speciestreefile = "input/species_tree_final.nw";

open(SPECIESTREE, $speciestreefile) or die("\nError: Couldn't open '$speciestreefile'\n\n");

our $usage = "";
args(0);

# start

Clear($speciestable);

# Get Compara species list:
# zcat input/Compara.protein_trees.nh.emf.gz | perl -ne '/^SEQ (\S+) /; print "$1\n";' | sort | uniq
# Longest species name:
# zcat input/Compara.protein_trees.nh.emf.gz | perl -ne '/^SEQ (\S+) /; print "$1\n";' | sort | uniq | wc -L


# The tree above 
# Read 


startme("Getting Compara Species Tree species order from '$speciestreefile', retrieving UniProt species mnemonic and taxon ID from table '$ensemblspeciestable', and inserting into table '$speciestable'");
$_ = <SPECIESTREE>;
chomp();
# Strip everything but species order & insert into table
# while (/([A-Za-z_]+)/g)
foreach $s (split(/,/, $_))
{
	# Strip tree formatting
	$s =~ s/:\d\.\d+(e-\d+)?//g;
	$s =~ s/[();]//g;
	$species = $s;
	
	addme("total ensembl fullspecies", $species);
	
	# Now removed in the tree itself (v3)
	# # Fix gallus_gallus and sus_scrofa, which have GCA assembly IDs attached to their "reference breeds":
	# # gallus_gallus_gca016699485v1
	# # sus_scrofa_gca000003025v6
	# $species =~ s/_gca\d+v\d+$//;
	
	# process
	# $species = lc($species);
	
	# Display name (e.g. Homo sapiens)
	$display = $species;
	# Make first letter uppercase
	substr($display, 0, 1) = uc(substr($display, 0, 1));
	# Convert underscores to spaces
	$display =~ s/_/ /g;
	
	# Shortened display name (e.g. H. sapiens)
	$shortdisplay = $species;
	# Make first letter uppercase
	substr($shortdisplay, 0, 1) = uc(substr($shortdisplay, 0, 1));
	# Abbreviate genus name
	$shortdisplay =~ s/^(\w)(\w+)_/$1. /;
	# Convert underscores to spaces
	$shortdisplay =~ s/_/ /g;
	
	# Check if this is really a Compara species, or a taxon like 'Fungi_Metazoa_group'
	$query = Query("SELECT ensp FROM `$enspspeciestable` WHERE species='$species' LIMIT 1");
	# $query = Query("SELECT id FROM `$enspspeciestable` WHERE species='$species' LIMIT 1");
	if (Numrows($query) > 0)
	{
		# # Get UniProt species, if available
		# $unispec = '';
		# print "SELECT e.species FROM `$enspspeciestable` s, uniens e WHERE s.species='$species' AND e.ensp=s.ensp GROUP BY e.species\n" if (switch('debug'));
		# $query = Query("SELECT e.species FROM `$enspspeciestable` s, uniens e WHERE s.species='$species' AND e.ensp=s.ensp GROUP BY e.species");
		# if (Numrows($query) > 0)
		# {
		# 	# Workaround for yeast strain problem (YEAST / YEASX are both S. cerevisiae)
		# 	if (Numrows($query) > 1)
		# 	{
		# 		if ($species eq 'saccharomyces_cerevisiae') { $unispec = 'YEAST'; }
		# 	}
		# 	else
		# 	{
		# 		($unispec) = FetchOne($query);
		# 		$unispec = uc($unispec);
		# 	}
		# }
		# # else
		# # {
		# # 	addme("species that aren't cross-mapped from UniProt yet (in table 'uniens')", $species);
		# # 	next;
		# # }
		#
		#         # Get NCBI taxon IDs
		# $tmpspec = $species;
		# substr($tmpspec, 0, 1) = uc(substr($tmpspec, 0, 1));
		# $tmpspec =~ s/_/ /g;
		# $query = Query("SELECT tax, speciestype, primate, mammal, vertebrate, metazoan, eukaryote FROM tax WHERE name='".esc($tmpspec)."'");
		# if (Numrows($query) == 0)
		# {
		# 	die("Error: Couldn't find NCBI Taxon ID in table 'tax' for species '$species' >> '$tmpspec'");
		# }
		# ($tax, $type, $primate, $mammal, $vertebrate, $metazoan, $eukaryote) = FetchOne($query);
		#         addme("total species types", $type);
		#         addme("species type '$type'", $species);
		#
		#
		# addme("unique uniprot species ids before non-uniens additions", $unispec);
		# # ...or set by hand if not
		# # SELECT * FROM comparaspecies ORDER BY unispec, display;a
		# # SELECT * FROM comparaspecies s, comparafasta f WHERE f.species=s.species AND s.unispec='' GROUP BY s.species;
		#
		# if ($unispec eq '')
		# {
		# 	# $query = Query("SELECT DISTINCT species FROM uniid WHERE type='NCBI_TaxID' AND value='$tax'");
		# 	$query = Query("SELECT DISTINCT species FROM unitax WHERE tax='$tax'");
		# 	if (Numrows($query) == 0)
		# 	{
		# 	    # Some manual overrides
		#                 if ($tax eq '59894')
		#                 {
		#                     $unispec = 'FICAL';
		#                     warn("Warning: Manually assigned TrEMBL-only UniProt species ID '$unispec' for NCBI_TaxID '$tax' for species '$species' >> '$tmpspec' (kept)");
		#                 }
		#                 elsif ($tax eq '48698')
		#                 {
		#                     $unispec = 'POEFO';
		#                     warn("Warning: Manually assigned TrEMBL-only UniProt species ID '$unispec' for NCBI_TaxID '$tax' for species '$species' >> '$tmpspec' (kept)");
		#                 }
		#                 elsif ($tax eq '7994')
		#                 {
		#                     $unispec = 'ASTMX';
		#                     warn("Warning: Manually assigned TrEMBL-only UniProt species ID '$unispec' for NCBI_TaxID '$tax' for species '$species' >> '$tmpspec' (kept)");
		#                 }
		#                 else
		# 	    {
		#                     die("Error: Couldn't find UniProt species ID in table 'uniid' for NCBI_TaxID '$tax' for species '$species' >> '$tmpspec'");
		#                     # warn("Warning: Couldn't find UniProt species ID in table 'uniid' for NCBI_TaxID '$tax' for species '$species' >> '$tmpspec' (this is normal to occur a few times - some aren't in Swiss-Prot, they're TrEMBL only) (kept)");
		#                     # $unispec = '';
		# 	    }
		# 	}
		# 	else
		# 	{
		#     			($unispec) = FetchOne($query);
		# 	}
		#
		# 	addme("species that aren't cross-mapped from UniProt yet (in table 'uniens')", "$species >> $tmpspec >> $tax >> $unispec");
		#             # die("Error: Can't get unispec for '$species'") if (!defined($unispec) or ($unispec eq ''));
		#             warn("Warning: Can't get unispec for '$species'") if (!defined($unispec) or ($unispec eq ''));
		#
		# 	# Old manual approach:
		# 	# if ($species eq 'dipodomys_ordii') { $unispec = 'DIPOR'; }
		# 	# elsif ($species eq 'gadus_morhua') { $unispec = 'GADMO'; }
		# 	# elsif ($species eq 'gasterosteus_aculeatus') { $unispec = 'GASAC'; }
		# 	# elsif ($species eq 'nomascus_leucogenys') { $unispec = 'NOMLE'; }
		# 	# elsif ($species eq 'petromyzon_marinus') { $unispec = 'PETMA'; }
		# 	# elsif ($species eq 'pteropus_vampyrus') { $unispec = 'PTEVA'; }
		# 	# else
		# 	# {
		# 	# 	warn("WARNING: There were some species without a UniProt species ID. Check list above and fix them via uniprot.org/taxonomy!");
		# 	# 	addme("species that I can't find a UniProt species ID for (fix me manually, via uniprot.org/taxonomy!)", $species);
		# 	# }
		# }
		
		# Get unispec and tax from ensembl_species
		$query = Query("SELECT species, tax FROM $ensemblspeciestable WHERE fullspecies='$species'");
		($unispec, $tax) = FetchOne($query);

		# Get species 'type' field (primate, mammal, vertebrate, other) and species categories ('primate', 'mammal', 'vertebrate', 'metazoan', 'other' fields)
		$query = Query("SELECT DISTINCT speciestype, homininae, greatape, primate, mammal, vertebrate, metazoan, eukaryote FROM tax WHERE tax='$tax'");
		if (Numrows($query) == 0)
		{
			die("Error: Couldn't find NCBI Taxon ID in table 'tax' for species '$species' (unispec '$unispec')'");
		}
		($type, $homininae, $greatape, $primate, $mammal, $vertebrate, $metazoan, $eukaryote) = FetchOne($query);

        addme("total types of species", $type);
		if (!defined($unispec))
		{
			$unispec = '';
	        addme("ensembl species with unispec=NULL for fullspecies", $species);
	        addme("ensembl species with unispec=NULL for '$type' for fullspecies", $species);
		}
		else
		{
			addme("total uniprot species", $species);
		
			addme("uniprot species successfully assigned for fullspecies", $species);
			addme("uniprot species successfully assigned for unispec", $unispec);
	        addme("uniprot species successfully assigned for type '$type' for unispec", $unispec);
		}
		
		
		
		# Insert into table 'comparaspecies'
		$q = "INSERT INTO `$speciestable` SET species='$species', display='$display', shortdisplay='$shortdisplay', unispec='$unispec', tax='$tax', type='$type', homininae='$homininae', greatape='$greatape', primate='$primate', mammal='$mammal', vertebrate='$vertebrate', metazoan='$metazoan', eukaryote='$eukaryote'";
		$q =~ s/=''/=NULL/g;
		Query($q);
		
		addme("unique species successfully inserted", $species);

		stepme(10);
	}
	else
	{
		addme("things that don't look like a species or aren't in any compara alignments (checked comparaenspspecies) (skipped)", $species);
		next;
	}
}
close(SPECIESTREE);
stopme();

# showmesomesorted(10);
showmeallsorted(1);
# showmesome(10);
# showmeall(1);

showme("things that don't look like a species or aren't in any compara alignments (checked comparaenspspecies) (skipped)");
showme("unique species successfully inserted");

# # Set types
# Now obtained from table 'tax'
# # Primates
# # Query("UPDATE `$speciestable` SET type='primate' WHERE species IN ('homo_sapiens', 'pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii', 'nomascus_leucogenys', 'macaca_mulatta', 'callithrix_jacchus', 'tarsius_syrichta', 'microcebus_murinus', 'otolemur_garnettii')");
# 
# Query("UPDATE `$speciestable` SET type='other'");
# Query("UPDATE `$speciestable` SET type='primate' WHERE id>=1 AND id<=10");
# Query("UPDATE `$speciestable` SET type='mammal' WHERE id>=11 AND id<=41");
# Query("UPDATE `$speciestable` SET type='vertebrate' WHERE id>=42 AND id<=61");
# 
# Query("UPDATE `$speciestable` SET type='' WHERE type IS NULL");
# 
# # Carry over more precise classification from older table (in this case blang_2012_01). I've done this mostly manually via NCBI Taxonomy.
# Query("UPDATE blang.comparaspecies s, blang_2012_01.comparaspecies p SET s.primate=p.primate, s.mammal=p.mammal, s.vertebrate=p.vertebrate, s.metazoan=p.metazoan, s.eukaryote=p.eukaryote WHERE s.species=p.species");

state("Inserted Compara Species Tree order for ".getme()." species into '$speciestable'");

done();
