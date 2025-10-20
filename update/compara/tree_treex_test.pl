#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
our $usage = "$0";
args(0);

# start

startme("Comparing Compara NHX trees vs. NH trees (Are they identical except for additional NHX information?)");
starttime();
chdir("tmp");
while (-e sprintf("%06d", getme() + 1).".tree.txt")
{
	$treefile = sprintf("%06d", getme() + 1).".tree.txt";
	$treexfile = sprintf("%06d", getme() + 1).".treex.txt";
	
	stepme(100);	
	
	open(TREE, $treefile) or die("Error: Couldn't open '$treefile'");
	$tree = '';
	while (<TREE>)
	{
		$tree .= $_;
	}
	close(TREE);

	open(TREEX, $treexfile) or die("Error: Couldn't open '$treexfile'");
	$treex = '';
	while (<TREEX>)
	{
		$treex .= $_;
	}
	close(TREEX);
	
	# Strip extra information from treex
	# [&&NHX:D=N:G=ENSMGAG00000008201:T=9103]
	$treex =~ s/\[&&NHX.+?\]//g;
	
	# Tiny 0 at the very end:
	# tree  ENSCVAP00000023609:0.696084):0.095707):0.032159);
	# treex ENSCVAP00000023609:0.696084):0.095707):0.032159):0;
	$treex =~ s/\):0;$/);/;
	
	# YBR289W:1.6117)Fungi/Metazoa_group:0.0000;
	# $treex =~ s/\)\S+?:/):/g;

	# # Compara 65: Change 0.0000 (NHX) to just 0 (NH)
	# $treex =~ s/:0.0000/:0/g;
	# # Compara 65: Change 100000.0000 (NHX) to just 100000 (NH)
	# $treex =~ s/:100000.0000/:100000/g;
	# # Compara 65: Change 100000.0000 (NHX) to just 100000 (NH)
	# $treex =~ s/:(\d{5,}).0000/:$1/g;
	
	# Compare tree vs treex
	if ($tree eq $treex)
	{
		addme("trees match for alignment", "$treefile");
	}
	else
	{
		# Compara 65: Change 0.0001 (NH) to 0 (NHX) (rounding error)
		# Compara 65: Change 0.0001 (NHX) to 0 (NH) (rounding error)
		# Try this:
		$tree =~ s/:0.0001/:0/g;
		$treex =~ s/:0.0001/:0/g;

		if ($tree eq $treex)
		{
			addme("trees match for alignment", "$treefile");
		}
		else
		{
			$treex =~ s/:1.0000/:1/g;

			if ($tree eq $treex)
			{
				addme("trees match for alignment", "$treefile");
			}
			else
			{
				addme("trees don't match for alignment", "$treefile");
		
				if (switch('debug'))
				{
					print "Mismatch:\n";
					print "\n";
					print "$treefile:  ['$tree']\n\n";
					print "$treexfile: ['$treex']\n\n";
			
					exit;
				}
			}
		}
	}
}
chdir("..");
stopme();
stoptime();

showme("trees match for alignment", 1);
showme("trees don't match for alignment", 1);

done();
