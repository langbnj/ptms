#!/usr/bin/env perl -w

# initialize

# our $superloudmysql = 1;

require('functions.inc.pl');
require('mysql.inc.pl');

starttime();

our $usage = "$0 [release number]\n\nExample: $0 107";
($rel) = args(1);

Clear('ensembl');
print "\n>> Cleared table 'ensembl'\n";


# run

# Count total
$totalspecies = chompme(`ls input_$rel/*.pep.all.fa.gz -1 | wc -l`);

open(DIR, "ls input_$rel -1|");
$myspecies = 0;
while (<DIR>)
{
	chomp;
	
	$pepfile = $_;

	# Skip directories
	next if (-d "input_$rel/$_");
	# Skip anything except .pep.all.fa.gz files
	next if (!/\.pep\.all\.fa\.gz$/);
	
    /^(\w+)\.[\w\.\-]+\.pep\.all\.fa\.gz$/ or die("Error: Couldn't match file name '$_'");
    # $species = lc($1);
	$species = $1;

	$myspecies++;

	# Now fixing this in download.pl
	# # Manual fix for an Ensembl 106 and 107 (but not 105) problem: Cyprinus_carpio clashes with the newer Cyprinus_carpio_carpio (and only the latter is included on the website, for a total of 314 species, one less than is on the FTP server and in ensembl_species)
	# if ($species eq 'Cyprinus_carpio')
	# {
	# 	addme("DEBUG: manually skipped species 'Cyprinus_carpio' since it wasn't on the species list at https://useast.ensembl.org/info/about/species.html, and it clashes with 'Cyprinus_carpio_carpio' (same ENSP identifiers etc.)", $species);
	# 	warn("Warning: DEBUG >> $myspecies / $totalspecies >> Skipping '$species' proteome from '$_' because: manually skipped species 'Cyprinus_carpio' since it wasn't on the species list at https://useast.ensembl.org/info/about/species.html, and it clashes with 'Cyprinus_carpio_carpio' (same ENSP identifiers etc.)\n");
	# 	next;
	# }

	# Check if this species is "permitted" (if it is in table 'ensembl_species')
	$query = Query("SELECT species FROM ensembl_species WHERE fullspecies='$species'");
	if (Numrows($query) == 0)
	{
		addme("skipped species because it wasn't in ensembl_species for fullspecies", $species);
		print " >> $myspecies / $totalspecies >> Skipping '$species' proteome from '$_' because this species isn't included in table 'ensembl_species'\n";
		next;
	}
	($unispec) = FetchOne($query);
	
	# No UniProt species mnemonic for this Ensembl species (it's not in Swiss-Prot yet, only in TrEMBL)
	if (!defined($unispec))
	{
		$unispec = '';
	}


	# Read CDS file for this species
	$cdsfile = $_;
	$cdsfile =~ s/\.pep\.all\.fa\.gz$/.cds.all.fa.gz/;
	print " >> $myspecies / $totalspecies >> Reading '$species' CDSs from '$cdsfile'\n";
	open(CDS, "zcat input_$rel/$cdsfile|") or die("Error: Couldn't zcat CDS file 'input_$rel/$cdsfile'");
	%cds = ();
	fastabreak();

	print "0..";
	$i = 0;
	$| = 1;

	while (<CDS>)
	{
	    ($cdstitle, $cds) = getfasta();
		
        # $title =~ /^(\S+) .*gene:(\S+)/ or die("Error: Couldn't match title '$title'");
        $cdstitle =~ /^(\S+) .*/ or die("Error: Couldn't match title '$cdstitle'");
		$enstv = $1; 
		# $ensgv = $2;
		
		# # Get gene symbol
		# $symbol = '';
        # # $title =~ / gene_symbol:(\S+) / or die("Error: Couldn't match gene_symbol in title '$title'");
		# if ($title =~ / gene_symbol:(\S+) /)
		# {
		# 	$symbol = $1;
		# }
		
        die("Error: $_") if (/Sequence unavailable/);
		
		# # Old-school alternative to getfasta():
		# # This also removes internal stop codons (*) (anything that's not \w).
		# >> Not necessary, no stop codons here
		# # $seq = '';
		# # @a = split(/\n/); shift(@a);
		# # foreach $s (@a) { $seq .= $s; }
		# # $seq =~ s/\W//g;
		# These CDSs contain just about every possible weird DNA FASTA character, except U, at least:
		# ~/update/ensembl/input_108 >> zcat *.cds.all.fa.gz | g -v "^>" | g "[^ACGTNRYMKSWBHVD]"
		# >> Returns nothing
		# ...but so did the Compara protein_default CDS file, actually:
		# ~/update/compara/input >> zcat Compara.108.protein_default.cds.fasta.gz | g -v "^>" | g "[^ACGT\-NRYWSMKHVDB/]"
		# >> Exactly the same non-standard DNA FASTA characters.
		# >> Should not cause any issues.
		
        # # $ensg = $ensgv;
        # $enst = $enstv;
		
		# # Strip version number
		# No need, not using ENST here (only ENSTV)
        # if (($species ne 'Drosophila_melanogaster') and ($species ne 'Saccharomyces_cerevisiae') and ($ensgv !~ /^MGP_/))
        # {
        #     # $ensg =~ s/\.\d+$//;
        #     $enst =~ s/\.\d+$//;
        # }

		$cds{$enstv} = $cds;
		
		$i++;
		if ($i%1000 == 0)
		{
			print "$i..";
		}
	}
	normalbreak();
	close(CDS);
	print "$i\n";


	

	print " >> $myspecies / $totalspecies >> Inserting '$species' proteome from '$pepfile'\n";
    # Query("INSERT INTO `species` SET species='$species'");
	# open(IN, "input/$pepfile") or die("Error: Couldn't open '$pepfile'\n");
	open(IN, "zcat input_$rel/$pepfile|") or die("Error: Couldn't zcat '$pepfile'\n");

	print "0..";
	$i = 0;
	$| = 1;

	fastabreak();
	while (<IN>)
	{
	    ($title, $seq) = getfasta();
		
        $title =~ /^(\S+) .*gene:(\S+) .*transcript:(\S+)/ or die("Error: Couldn't match title '$title'");
		$ensgv = $2; $enstv = $3; $enspv = $1;
		
		# Get gene symbol
		$symbol = '';
        # $title =~ / gene_symbol:(\S+) / or die("Error: Couldn't match gene_symbol in title '$title'");
		if ($title =~ / gene_symbol:(\S+) /)
		{
			$symbol = $1;
		}
		
        die("Error: '$_' in '$pepfile'") if (/Sequence unavailable/);
		
		# Old-school alternative to getfasta():
		# This also removed internal stop codons (*) (anything that's not \w).
		# $seq = '';
		# @a = split(/\n/); shift(@a);
		# foreach $s (@a) { $seq .= $s; }

		# This would remove internal stop codons (no longer used):
		# $seq =~ s/\W//g;

		# The sequences contain all of these: XUZB*
		# That's it though:
		# ~/update/ensembl/input_108 >> zcat *.pep.all.fa.gz | g -v "^>" | g "[^ACDEFGHIKLMNPQRSTVWYXUZB\*]"
		# >> Returns nothing.
		# UniProt contains XUZB as well, plus O (pyrrolysine), but no stop codons:
		# SELECT * FROM uniseq WHERE type IN ('UniProt', 'UniIso') AND seq REGEXP '[^ACDEFGHIKLMNPQRSTVWYXUBZO]';
		# >> Returns nothing.
		# >> Can't really filter out any characters, no point.
		# die("Error: Sequence contains unexpected non-AA characters") if ($seq =~ /[^]/)
		# Also, leaving stop codons in would not improve agreement with UniProt since UniProt never has these.
		
        $ensg = $ensgv;
        $enst = $enstv;
        $ensp = $enspv;
		
        # # Note - Release 84 - CAEEL, DROME and YEASX don't have versioned ENSGs
        # if (($species ne 'Caenorhabditis_elegans') and ($species ne 'Drosophila_melanogaster') and ($species ne 'Saccharomyces_cerevisiae'))
        # Note - Release 107 - Drosophila, Yeast, and various mouse strains (but not Mus_musculus itself, nor Mus_spicilegus) and Mus caroli/pahari/spretus don't have versioned ENSG/T/Ps

		# Release 108:

		# # ENS%.1
		# SELECT * FROM ensembl WHERE enspv NOT REGEXP '^ENS[A-Z]+[0-9]+\.[0-9]+$';
		# SELECT DISTINCT species FROM ensembl WHERE enspv NOT REGEXP '^ENS[A-Z]+[0-9]+\.[0-9]+$';
		# SELECT DISTINCT fullspecies FROM ensembl WHERE enspv NOT REGEXP '^ENS[A-Z]+[0-9]+\.[0-9]+$';
		# 
		# # %.1
		# SELECT * FROM ensembl WHERE enspv NOT REGEXP '\.[0-9]+$';
		# SELECT DISTINCT species FROM ensembl WHERE enspv NOT REGEXP '^\.[0-9]+$';
		# SELECT DISTINCT fullspecies FROM ensembl WHERE enspv NOT REGEXP '^\.[0-9]+$';
		
		# CAEEL
		# DROME
		# MUSCR
		# MUSMC
		# MUSPA
		# MUSSP
		# YEAST

		# Caenorhabditis_elegans
		# Drosophila_melanogaster
		# Mus_caroli
		# Mus_musculus_129s1svimj
		# Mus_musculus_aj
		# Mus_musculus_akrj
		# Mus_musculus_balbcj
		# Mus_musculus_c3hhej
		# Mus_musculus_c57bl6nj
		# Mus_musculus_casteij
		# Mus_musculus_cbaj
		# Mus_musculus_dba2j
		# Mus_musculus_fvbnj
		# Mus_musculus_lpj
		# Mus_musculus_nodshiltj
		# Mus_musculus_nzohlltj
		# Mus_musculus_pwkphj
		# Mus_musculus_wsbeij
		# Mus_pahari
		# Mus_spretus
		# Saccharomyces_cerevisiae

		# >> Still the same as in release 107

		# # Run this to verify the ensembl table:
		# # SELECT fullspecies, COUNT(DISTINCT ensgv) AS ensgvs, COUNT(DISTINCT REGEXP_REPLACE(ensgv, '\\.\\d$', '')) AS ensgs, COUNT(DISTINCT enstv) AS enstvs, COUNT(DISTINCT REGEXP_REPLACE(enstv, '\\.\\d$', '')) AS ensts, COUNT(DISTINCT enspv) AS enspvs, COUNT(DISTINCT REGEXP_REPLACE(enspv, '\\.\\d$', '')) AS ensps FROM blang.ensembl GROUP BY fullspecies HAVING ensgvs!=ensgs OR enstvs!=ensts OR enspvs!=ensps;
		# SELECT fullspecies, COUNT(DISTINCT ensgv) AS ensgvs, COUNT(DISTINCT ensg) AS ensgs, COUNT(DISTINCT enstv) AS enstvs, COUNT(DISTINCT enst) AS ensts, COUNT(DISTINCT enspv) AS enspvs, COUNT(DISTINCT ensp) AS ensps FROM blang.ensembl GROUP BY fullspecies HAVING ensgvs!=ensgs OR enstvs!=ensts OR enspvs!=ensps;
		# >> returns nothing
		# This query should return nothing, otherwise unversioned ensgs, ensts or ensps have lost uniqueness due to \.\d+$ stripping.
		# 
		# To see unversioned species (actually also includes ENSPV for the MGP_... Mus species, but those are not in Compara, so I don't mind):
		# SELECT * FROM ensembl GROUP BY fullspecies HAVING ensg NOT LIKE 'ENS%';
        # 
		# if (($species ne 'drosophila_melanogaster') and ($species ne 'saccharomyces_cerevisiae') and ($ensgv !~ /^MGP_/))

		# CAEEL, DROME, and YEAST are the unversioned species.
        if (($species ne 'Caenorhabditis_elegans') and ($species ne 'Drosophila_melanogaster') and ($species ne 'Saccharomyces_cerevisiae'))
        {
            $ensg =~ s/\.\d+$//;
            $enst =~ s/\.\d+$//;
            $ensp =~ s/\.\d+$//;
        }
		
		die("Error: No CDS in '$cdsfile' for ENSTV '$enstv'") if (!exists($cds{$enstv}));
		$cds = $cds{$enstv};

		# Skip any rows where the AA sequence contains non-standard-20 characters (X/B/Z/U/*)
		if (!aa($seq))
		{
			addme("AA sequence contains non-standard-20 characters for species (skipped)", $species);
			addme("AA sequence contains non-standard-20 characters for ensp (skipped)", $ensp);
			next;
		}
		# Skip any rows where the CDS contains non-ACGT characters
		if (!dna($cds))
		{
			addme("AA sequence is standard-20, but CDS contains non-ACGT characters for species (skipped)", $species);
			addme("AA sequence is standard-20, CDS contains non-ACGT characters for ensp (skipped)", $ensp);
			next;
		}


		# Check if CDS translates to AA sequence, and if not, try to trim it to repair it, otherwise leave it as NULL
		# if (aa($seq) and dna($cds))
		# {
			# if (length($seq) != length($cds) / 3)
			# if (98 != 271 / 3)
			if (length($cds) != length($seq) * 3)
			{
				# addme("cds length / 3 doesn't match expected seq length for expected seq (skipped)", $seq);
				# $mismatch++;
				# next;

				$overhang = length($cds) - (length($seq) * 3);
				# print "\nOVERHANG >> $overhang\n\n";
				if (($overhang >= 1) and ($overhang <= 3))
				{
					# If the overhang is between +1 and +3: simply remove it
					addme("cds length / 3 doesn't match expected seq length for species (kept, trimmed)", $species);
					addme("cds length / 3 doesn't match expected seq length for ensp (kept, trimmed)", $ensp);
					addme("cds length / 3 doesn't match expected seq length for ensp (kept, trimmed by $overhang)", $ensp);
					
					$cds = substr($cds, 0, length($seq) * 3);
				}
				else
				{
					# Overhang not fixable
					addme("cds length / 3 doesn't match expected seq length for species (skipped, couldn't be trimmed)", $species);
					addme("cds length / 3 doesn't match expected seq length for ensp (skipped, couldn't be trimmed)", $ensp);
					next;
				}
			}

			# Translate CDS
			# $transeq = translate_loosely($cds);
			$transeq = translate($cds);

			# # Remove terminal stop codon from translation (the CDS doesn't always include this, and the expected sequence never does)
			# $transeq =~ s/\*$//;

			# Check if CDS translation matches expected seq
			if ($transeq ne $seq)
			{
				# addme("translation mismatch for species (set CDS to NULL)", $species);
				# addme("translation mismatch for ensp (set CDS to NULL)", $ensp);
				addme("translation mismatch for species (skipped)", $species);
				addme("translation mismatch for ensp (skipped)", $ensp);
				# $mismatch++;
				# $cds = '';
				next;
			}
			# else
			# {
			# 	addme("translation match for expected seq", $seq);
			# 	$match++;
			# }
		# }
		# else
		# {
		# 	# Either the AA sequence or the CDS contain weird characters (non-standard-20 aa, or non-ACGT)
		# 	# NULL
		# 	$cds = '';
		# }

		$q = "INSERT INTO `ensembl` SET fullspecies='$species', ensgv='$ensgv', ensg='$ensg', enstv='$enstv', enst='$enst', enspv='$enspv', ensp='$ensp', seq='$seq', cds='$cds', symbol='".esc($symbol)."'";
		$q =~ s/=''/=NULL/g;
		Query($q);
		$i++;
		if ($i%1000 == 0)
		{
			print "$i..";
		}
	}
	close(IN);
	normalbreak();
	
	addme("inserted proteins for species", $species);

	print "$i\n";
	print "  >> ".commify($i)." '$species' proteins inserted into 'ensembl'\n";
}

stoptime();

showmeall(1);

# Doesn't seem to make any speed difference
# state("Enabling keys...");
# starttime();
# Query("ALTER TABLE ensembl ENABLE KEYS");
# stoptime();

done();
