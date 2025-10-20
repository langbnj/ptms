#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$maxbranch = '100';

our $usage = "$0 [species] [MAFFT mode] [-tree]\n\n Species: Reference species.\n MAFFT mode: linsi, einsi or ginsi (for L-INS-i, E-INS-i or G-INS-i)\n".
" -tree: Use MAFFT results with species tree.\n -1para: Keeps the best-matching outparalog (by various metrics) for one2many and many2many cases.\n\nExample: $0 human linsi -tree";
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



# Build output directory name
$tmp = $mafftmode;
$tmp =~ s/\./_/g;
$outdir = "./output_".$tmp."/";

if (!-d $outdir)
{
	mkdir($outdir) or die("Error: Couldn't create output directory '$outdir'");
}

# Open output files
$fastaout = $outdir."Compara.protein_trees.aa.fasta";
$treeout  = $outdir."Compara.protein_trees.nh.emf";

# Open sequence output file
open(FASTA, ">$fastaout") or die("Error: Couldn't open FASTA output file '$fastaout'");

# Open tree output file
open(TREE, ">$treeout") or die("Error: Couldn't open tree output file '$treeout'");



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

startme("Parsing results for '$unispec' ($species) using MAFFT alignments of type '$mafftmode'", 0, Numrows($mainquery));
starttime();

while (($ensp) = Fetch($mainquery))
{
	# Open input file for this ENSP
	if (switch('1para') or switch('1parahc'))
	{
		$infile = "input/$label/$ensp.$label.txt";
	}
	else
	{
		$infile = "input/$ensp.txt";
	}
	print("   >> $infile\n") if (switch('debug'));
	if (!-s $infile)
	{
		addme("FASTA file not found for ensp (skipped)", $ensp);
		next;
	}
	
	# Check MAFFT alignment file
	$mafftfile = "output/mafft.$ensp.$mafftmode.txt";
	if (!-s $mafftfile)
	{
		addme("MAFFT file not found for ensp (skipped)", $ensp);
		next;
	}

	# Check output CDS alignment file
	$cdsfile = "output/cds.$ensp.$mafftmode.txt";
	if (!-s $cdsfile)
	{
		addme("back-translated CDS file not found for ensp (skipped)", $ensp);
		next;
	}
	
	# Check output CDS alignment file
	$cdstaxfile = "output/cdstax.$ensp.$mafftmode.txt";
	if (!-s $cdstaxfile)
	{
		addme("back-translated CDS file with taxon IDs not found for ensp (skipped)", $ensp);
		next;
	}
	
	# Check TreeBeST tree file
	$treefile = "output/tree.$ensp.$mafftmode.txt";
	if (!-s $treefile)
	{
		addme("TreeBeST output tree not found for ensp (skipped)", $ensp);
		next;
	}
    # Have to do fixtree()
    #
    # Dendroscope doesn't like it otherwise, though Bio::Phylo seems quite okay with it.
    # Branch lengths and constituents remain intact, I checked.
    # NHX bootstrap info gets dropped this way (not needed for evolutionary rates)
    #
	# state("TREEFILE $treefile");
	open(TREEIN, $treefile) or die("Error: Couldn't open '$treefile'");
	$tree = '';
	while (<TREEIN>)
	{
		$tree .= $_;
	}
	close(TREEIN);
	# state("TREE $tree");
	$tree = fixtree_compara($tree);
	while ($tree =~ /:(\d+\.\d+)/g)
	{
		if ($1 >= $maxbranch)
		{
			# This would replace branch lengths greater than $maxbranch with $maxbranch
			# s/:$1/:$maxbranch/;
			
			# die("Error: Branch length $1 is greater than $maxbranch in tree '$tree' for ENSP '$ensp'");
            warn("Warning: Branch length $1 is greater than $maxbranch in tree '$tree' for ensp '$ensp' (kept)");

            addme("branch length greater than $maxbranch for ensp (kept)", $ensp);
            # next;
		}
	}
	
	
	# Write to sequence output file
	open(MAFFT, $mafftfile) or die("Error: Couldn't open '$mafftfile'");
	fastabreak();
	@titles = ();
	while (<MAFFT>)
	{
		($title, $seq) = getfasta();
		
		# Remove NCBI taxon ID at end of fasta title
		$title =~ s/_\d+$//;
		push(@titles, $title);
		
		print FASTA ">$title\n".split60($seq)."\n";
	}
	normalbreak();
	print FASTA "\n//\n\n";
	
	# Write sequence titles to tree output file (only to make /evorate/split.pl happy, all it checks for is: "while (/^SEQ \S+ (\S+) /mg)")
	foreach $title (@titles)
	{
		print TREE "SEQ tmp $title tmp\n";
	}
	# Write tree to tree output file
	# 
	# Remove quotes around protein names (introduced by Bio::Phylo)
	$tree =~ s/'//g;
	# Write tree to file
	print TREE "DATA\n";
	print TREE "$tree\n";
	print TREE "//\n\n";
	
	stepme(100);
}
stopme();
stoptime();

print "Wrote to '$fastaout'\n";
print "Wrote to '$treeout'\n\n";

run("Compress output", "gzip $fastaout");
run("Compress output", "gzip $treeout");

showmeall();

done();
