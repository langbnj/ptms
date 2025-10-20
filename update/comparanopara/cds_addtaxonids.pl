#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species] [MAFFT mode] [-tree] [-1para]\n\n Species: Reference species.\n MAFFT mode: linsi, einsi or ginsi (for L-INS-i, E-INS-i or G-INS-i)\n -tree: Use MAFFT results with species tree.\n -1para: Keeps the best-matching outparalog (by various metrics) for one2many and many2many cases.\n\nExample: $0 human linsi -tree";
($unispec, $mafftmode) = args(2);

if (switch('tree'))
{
	$mafftmode .= ".tree";
}

# Handling for 1para
if (switch('1parahc'))
{
	$label = '1parahc';
	$mafftmode .= '.'.$label;
}
elsif (switch('1para'))
{
	$label = '1para';
	$mafftmode .= '.'.$label;
}



# start

$query = Query("SELECT species FROM comparaspecies WHERE unispec='$unispec'");
($species) = FetchOne($query);

if (!switch('1para') and !switch('1parahc'))
{
	# normal
	# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e, unimod m WHERE c.species='$species' AND e.ensp=c.ensp AND m.name=e.name AND m.ptm!='' GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
	# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e WHERE c.species='$species' AND e.ensp=c.ensp GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
	$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND homology='ortholog_one2one' AND hc=1 AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1");
}
elsif (switch('1parahc'))
{
	# 1parahc
	# 1para
	if (!switch('debug'))
	{
		$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND hc=1 AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1");
	}
	else
	{
		$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND hc=1 AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1 LIMIT 1");
	}
}
else
{
	# 1para
	if (!switch('debug'))
	{
		$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1");
	}
	else
	{
		$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1 LIMIT 1");
	}
}

startme("Adding taxon IDs from table 'comparaspecies' to CDS for '$unispec' ($species) mode '$mafftmode' to cDNA alignments", 0, Numrows($mainquery));
starttime();

while (($ensp) = Fetch($mainquery))
{
	stepme(10);
	
	$cdsinfile = "output/cds.$ensp.$mafftmode.txt";
	$cdsoutfile = "output/cdstax.$ensp.$mafftmode.txt";

	print(" >> $cdsinfile\n") if (switch('debug'));
	print("   >> $cdsoutfile\n") if (switch('debug'));

	if (!-s $cdsinfile)
	{
		addme("cds file not found for ensp (skipped)", $ensp);
		next;
	}

    if (-s $cdsoutfile)
    {
	    # Skip this ensp if output file already exists and is non-zero
        addme("output file already exists for ensp (skipped)", $ensp);
        next;
		# # Crash instead
		# die("Error: Output file '$cdsoutfile' already exists for ensp '$ensp' (should run clean.pl first)");
    }

	open(CDS, $cdsinfile) or die("Error: Couldn't open '$cdsinfile' for reading");
	fastabreak();
	$cds = '';
	while (<CDS>)
	{
		($title, $seq) = getfasta();
		
		# remove any existing ncbi ids from title (e.g. if partially run already)
		$title =~ s/_\d+$//;
		
		# $query = Query("SELECT t.tax FROM comparaenspspecies e, comparaspecies s, tax t WHERE e.species=s.species AND t.name=s.display AND e.ensp='$title'");
		$query = Query("SELECT s.tax FROM comparaenspspecies e, comparaspecies s WHERE e.species=s.species AND e.ensp='$title'");
		if (Numrows($query) == 0)
		{
		    die("Error: No taxon id in 'tax' for ensp '$title'");
            # addme("no taxon id in tax for title (kept)", $title);
            # addme("no taxon id in tax for ensp (kept)", $ensp);
		}
		else
		{
			($tax) = FetchOne($query);
			$title .= "_$tax";
		}
	
		$cds .= ">".$title."\n".split60($seq)."\n";

		addme("wrote output for ensp (ok)", $ensp);
	}
	close(CDS);
	open(CDS, ">$cdsoutfile") or die("Error: Couldn't open '$cdsoutfile' for writing");
	normalbreak();
	print CDS $cds;
	close(CDS);
}
stopme();
stoptime();

showmeall();

print("\nNote: ENSP00000467141 is titin and always crashes MAFFT.\n\n");

done();
