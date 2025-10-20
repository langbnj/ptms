#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [input tree]\n\nExample: $0 input/species_tree.branch_len.reordered_in_dendroscope.yeastrooted.taxids.nw";
($infile) = args(1);

$outfile = $infile;
$outfile =~ s/\.nw$/.internal_nodes_labelled.nw/;

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

state("Reading '$infile', labelling internal nodes (at closing parentheses, ')') and writing to '$outfile':");
$nl = $/;
$/ = '';
$tree = <IN>;
$/ = $nl;
close(IN);

print "Unlabelled tree:\n\n$tree\n\n";

# Labels will be aa, ab, ac, etc.
# ASCII value of a: 97
# ASCII value of z: 122
$c1 = 'a';
$c2 = 'a';
$newtree = '';
$added = 0;
while ($tree =~ /(.+?\))/g)
{
	# print " >> ".pos($tree)."\n";
	# print "   >> ".substr($tree, pos($tree)-1, 1)."\n";
	$newtree .= substr($1, $added);
	$added = length($1);
	# print "   >> ".$1."\n";
	substr($tree, pos($tree)-1, 1) = "}";
	# print "  >> $c1$c2\n";
	$newtree .= "$c1$c2";
	# print "    >> $newtree\n";

	# Step second character
	$c2 = chr(ord($c2) + 1);

	# Step first character and reset second to 'a'
	# If z has been passed:
	if (ord($c2) > 122)
	{
		$c1 = chr(ord($c1) + 1);
		$c2 = 'a';
	}
}
$newtree .= substr($tree, $added);

print "Labelled tree:\n\n$newtree\n\n";

print OUT $newtree;
close(OUT);

state("Wrote to '$outfile'");

done();
