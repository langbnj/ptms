#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'tax';

# our $usage = "";
# ($var) = args(1);
args(0);

$infile = "input/nodes.dmp";
#$outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
#open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


Query("UPDATE tax SET `rank`=NULL");
state("Cleared 'rank' field in table 'tax'");



# start
$lines = `cat '$infile' | wc -l`;
chomp($lines);

startme("Reading taxon ID rank (family, genus etc.) and its immediate parent taxon ID from '$infile'", 0, $lines);
starttime();

%rank = ();
%parent = ();
while (<IN>)
{
	chomp;
	
	@a = split(/\t/);
	
	$tax = $a[0];
	$parent = $a[2];
	$rank = $a[4];
	
	addme("total taxon IDs", $tax);
	addme("total parent taxon IDs", $parent);
	addme("total ranks", $rank);

	# next if ($type ne 'scientific name');
	
    if (exists($parent{$tax}))
    {
        die("Error: Multiple parents for taxon ID '$tax' (this means the taxon ID was observed multiple times in nodes.dmp)");
    }
    else
    {
        $parent{$tax} = $parent;
    }

    if (exists($rank{$tax}))
    {
        die("Error: Multiple ranks for taxon ID '$tax' (this means the taxon ID was observed multiple times in nodes.dmp)");
    }
    else
    {
        $rank{$tax} = $rank;
    }
    
    # last if ((getme() % 10000 == 0) and (switch('debug')));

	stepme(10000);
}
stopme();
stoptime();

# Follow up & collapse parent lineages
state("Following up & collapsing immediate parents into lineages until convergence on taxon ID '1' is reached:");
%parents = ();
$convergence = 0;
$round = 1;
$longest = 0;
while ($convergence == 0)
{
    state(" >> Round $round:", 1);
    startme("   >> Adding parents' parents", 1, scalar(keys(%parent)));
    foreach $tax (keys(%parent))
    {
        state("\n     >> Looking at '$tax'", 1) if (switch('debug2'));
        if (exists($parents{$tax}))
        {
            $parents{$tax} .= "|".$parent{$tax};
            state("       >> Immediate parent '".$parent{$tax}."'") if (switch('debug2'));
            state("         >> Child '$tax''s collected parents: '".$parents{$tax}."'", 1) if (switch('debug2'));
            
            $parent = $parent{$tax};
            next if ($parent eq '1');
            while (exists($parent{$parent}))
            {
                $parents{$tax} .= "|".$parent{$parent};
                state("       >> Parent '$parent''s parent: '".$parent{$parent}."'", 1) if (switch('debug2'));
                state("         >> Child '$tax''s collected parents: '".$parents{$tax}."'", 1) if (switch('debug2'));
                $parent = $parent{$parent};
                last if ($parent eq '1');
            }
        }
        else
        {
            $parents{$tax} = $parent{$tax};
            state("       >> Immediate parent '".$parent{$tax}."'") if (switch('debug2'));
            state("         >> Child '$tax''s collected parents: '".$parents{$tax}."'", 1) if (switch('debug2'));

            $parent = $parent{$tax};
            next if ($parent eq '1');
            while (exists($parent{$parent}))
            {
                $parents{$tax} .= "|".$parent{$parent};
                state("       >> Parent '$parent''s parent: '".$parent{$parent}."'", 1) if (switch('debug2'));
                state("         >> Child '$tax''s collected parents: '".$parents{$tax}."'", 1) if (switch('debug2'));
                $parent = $parent{$parent};
                last if ($parent eq '1');
            }
        }

        if (length($parents{$tax}) > $longest)
        {
            $longest = length($parents{$tax});
        }

        stepme(10000, 1);
    }
    stopme(1);
    
    # test convergence
    startme("    >> Testing convergence", 1, scalar(keys(%parent)));
    $convergence = 1;
    foreach $tax (keys(%parent))
    {
        stepme(10000, 1);
        if (!contains('1', split(/\|/, $parents{$tax})))
        {
            stopme(1);
            state("     >> Convergence broken by taxon id '$tax' (parents '".$parents{$tax}."') (lineage not ending in 1), starting new round", 1);
            $convergence = 0;
            $round++;
            last;
        }
    }
}
stopme(1);
state("Reached convergence (parents always end at taxon id 1) after $round rounds!");
state("Longest parents string: $longest characters");



startme("Filling 'rank' and 'parents' fields in table '$table'", 0, scalar(keys(%parents)));
starttime();
$affected = 0;
foreach $tax (keys(%parents))
{
    # $tax = 9606;
    
    # Assign species type
    # 
    # SELECT * FROM blang.tax WHERE name='Eukaryota';
    # SELECT * FROM blang.tax WHERE name='Metazoa';
    # SELECT * FROM blang.tax WHERE name='Vertebrata';
    # SELECT * FROM blang.tax WHERE name='Mammalia';
    # SELECT * FROM blang.tax WHERE name='Primate';
    # 
    # SELECT * FROM tax WHERE tax='2759';       # eukaryote
    # SELECT * FROM tax WHERE tax='33208';      # metazoan
    # SELECT * FROM tax WHERE tax='7742'; -- SELECT * FROM tax WHERE tax='1261581';     # vertebrate (1261581 is a genus of sponges that are unfortunately called "Vertebrata", ignore it)
    # SELECT * FROM tax WHERE tax='40674';      # mammal
    # SELECT * FROM tax WHERE tax='9443';       # primate

    $speciestype = 'other';
    @parents = split(/\|/, $parents{$tax});
    # show(@parents);
    
    # Manually assign speciestype:

    # SELECT * FROM tax WHERE tax=2759;		# eukaryote		# Superkingdom (according to https://www.uniprot.org/taxonomy/2759)
    # SELECT * FROM tax WHERE tax=33208;	# metazoan		# Kingdom
    # SELECT * FROM tax WHERE tax=7742;		# vertebrate	# "No rank", but directly below (and the only descendant of) Craniata, which is a Subphylum
    # SELECT * FROM tax WHERE tax=40674;	# mammal		# Class
    # SELECT * FROM tax WHERE tax=9443;		# primate		# Order
    # SELECT * FROM tax WHERE tax=9604;		# greatape		# Family
    # SELECT * FROM tax WHERE tax=207598;	# homininae		# Subfamily (also called the "Homo/Pan/Gorilla group", excludes Pongo (Orang-utan))

    if (contains('2759', @parents)) { $speciestype = 'eukaryote'; };
    if (contains('33208', @parents)) { $speciestype = 'metazoan'; };
    if (contains('7742', @parents)) { $speciestype = 'vertebrate'; }
    if (contains('40674', @parents)) { $speciestype = 'mammal'; }
    if (contains('9443', @parents)) { $speciestype = 'primate'; }
    if (contains('9604', @parents)) { $speciestype = 'greatape'; }
    if (contains('207598', @parents)) { $speciestype = 'homininae'; }
    
    addme("before simplifying: taxon ids for speciestype '$speciestype'", $tax);
    
    $eukaryote = 0;
    $metazoan = 0;
    $vertebrate = 0;
    $mammal = 0;
    $primate = 0;
    $greatape = 0;
    $homininae = 0;
    
    if ($speciestype eq 'eukaryote')    { $eukaryote = 1; $speciestype = 'other'; }
    if ($speciestype eq 'metazoan')     { $eukaryote = 1; $metazoan = 1; $speciestype = 'other'; }
    if ($speciestype eq 'vertebrate')   { $eukaryote = 1; $metazoan = 1; $vertebrate = 1; }
    if ($speciestype eq 'mammal')       { $eukaryote = 1; $metazoan = 1; $vertebrate = 1; $mammal = 1; }
    if ($speciestype eq 'primate')      { $eukaryote = 1; $metazoan = 1; $vertebrate = 1; $mammal = 1; $primate = 1; }
    if ($speciestype eq 'greatape')      { $eukaryote = 1; $metazoan = 1; $vertebrate = 1; $mammal = 1; $primate = 1; $greatape = 1; }
    if ($speciestype eq 'homininae')      { $eukaryote = 1; $metazoan = 1; $vertebrate = 1; $mammal = 1; $primate = 1; $greatape = 1; $homininae = 1; }

    addme("after simplifying: taxon ids for speciestype '$speciestype'", $tax);
    addme("total species types (after simplifying)", $speciestype);
    

    # Update tax table
    if (switch('debug'))
    {
        state("UPDATE `$table` SET `rank`='".$rank{$tax}."', parents='".$parents{$tax}."', speciestype='$speciestype', homininae=$homininae, greatape=$greatape, primate=$primate, mammal=$mammal, vertebrate=$vertebrate, metazoan=$metazoan, eukaryote=$eukaryote WHERE tax='$tax'", 1);
        $query = Query("SELECT id FROM `$table` WHERE tax='$tax'");
        $affected += Numrows($query);
    }
    else
    {
        $query = Query("UPDATE `$table` SET `rank`='".$rank{$tax}."', parents='".$parents{$tax}."', speciestype='$speciestype', homininae=$homininae, greatape=$greatape, primate=$primate, mammal=$mammal, vertebrate=$vertebrate, metazoan=$metazoan, eukaryote=$eukaryote WHERE tax='$tax'");
        $affected += Numrows($query);
    }

    # last;

    stepme(10000);
}
stopme();
stoptime();

# Get total number of rows in tax table for display
$query = Query("SELECT COUNT(id) FROM `$table`");
($count) = FetchOne($query);

nl();
state("$count total rows in table '$table'", 1);
state("$affected rows affected", 1);
nl();

Optimize($table);

done();
