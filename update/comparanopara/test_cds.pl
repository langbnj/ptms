#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species]\n\n Species: Reference species.\n\nExample: $0 human";
($unispec) = args(1);



# start

$query = Query("SELECT species FROM comparaspecies WHERE unispec='$unispec'");
($species) = FetchOne($query);

# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e, unimod m WHERE c.species='$species' AND e.ensp=c.ensp AND m.name=e.name AND m.ptm!='' GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
$mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e WHERE c.species='$species' AND e.ensp=c.ensp GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
startme("Testing CDS agreement of the API CDSes with the /update/compara/input CDSes for '$unispec' ($species)");
starttime();

while (($ensp, $aln) = Fetch($mainquery))
{
	stepme(100);
	
	addme("total ensps", $ensp);
	addme("total alns", $aln);

	$apifile = "input/$ensp.cds.txt";
	$comparafile = "../compara/tmp/".sprintf("%06d", $aln).".cds.txt";
	
	# system("cat $apifile");
	# system("cat $comparafile");
	# exit;
	
	if (!-e $apifile)
	{
		addme("apifile not found for ensp", $ensp);
		addme("apifile not found for aln", $aln);
		next;
	}
	
	if (!-e $comparafile)
	{
		addme("comparafile not found for ensp", $ensp);
		addme("comparafile not found for aln", $aln);
		next;
	}
	
	open(API, $apifile) or die("Error: Couldn't open '$apifile'");
	fastabreak();
	%api = ();
	while (<API>)
	{
		($title, $seq) = getfasta();
		
		# Remove stop codons
		# $seq =~ s/(TAA|TAG|TGA)$//;
				
		$api{$title} = $seq;
	}
	close(API);
	
	open(COMP, $comparafile) or die("Error: Couldn't open '$comparafile'");
	fastabreak();
	%comp = ();
	while (<COMP>)
	{
		($title, $seq) = getfasta();
		
		$comp{$title} = $seq;
	}
	close(COMP);
	
	@apikeys = sort(keys(%api));
	@compkeys = sort(keys(%comp));
	
	# Obviously there'd be a count mismatch (since I've pruned out the paralogs)
	# if (scalar(@apikeys) != scalar(@compkeys))
	# {
	# 	addme("count mismatch for ensp", $ensp);
	# 	addme("count mismatch for aln", $aln);
	# }
	
	foreach $key (@apikeys)
	{
		if (!contains($key, @compkeys))
		{
			addme("api title not found in compara file for ensp", $ensp);
			addme("api title not found in compara file for aln", $aln);
		}
		else
		{
			addme("api title successfully found in compara file for ensp", $ensp);
			addme("api title successfully found in compara file for aln", $aln);
		}
	}
	
	foreach $key (@apikeys)
	{
		if (!exists($comp{$key}))
		{
			addme("api ensp not found in comparafasta for ensp (skipped)", $ensp);
			addme("api ensp not found in comparafasta for aln (skipped)", $aln);
			next;
		}

		$comp{$key} =~ s/-//g;
		
		if ($api{$key} ne $comp{$key})
		{
			# See if the difference is just a stop codon in the API CDS:
			$tmp = $api{$key};
			$tmp =~ s/(TAA|TAG|TGA)$//;
			if ($tmp eq $comp{$key})
			{
				addme("sequence match for ensp", $ensp);
				addme("sequence match for aln", $aln);
			}
			else
			{
				addme("sequence mismatch for ensp", $ensp);
				addme("sequence mismatch for aln", $aln);

				state("API sequence (".length($api{$key})." bp):\n".$api{$key});
				state("Compara sequence (".length($comp{$key})." bp):\n".$comp{$key});
			}
		}
		else
		{
			addme("sequence match for ensp", $ensp);
			addme("sequence match for aln", $aln);
		}
	}
	
	# last if (getme() == 100);
}
stopme();
stoptime();

# showmesome(0);
showmeall(1);

done();
