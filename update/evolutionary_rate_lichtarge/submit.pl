#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');
# use List::Util 'shuffle';

# initialize
# $nseqmin = 10;
$nseqmin = 1;
# $type = "lichtarge";

our $usage = "$0 [species] [type] [-all]\n\n -all: run on ALL ENSPs, rather than just the ones that have a non-NULL ptm field in unimod\n\nExample: $0 human linsi_tree";
($species, $type) = args(2);

startme("Running lichtarge on '$species' '$type' alignments from Compara");
starttime();
chdir("tmp");
@aln = `ls -1 *.msf.txt`;
# @aln = shuffle(@aln);
$free = 0;
while (-e sprintf("%06d", getme() + 1).".msf.txt")
{
	$msffile = sprintf("%06d", getme() + 1).".msf.txt";
	$aln = sprintf("%06d", getme() + 1);

	stepme(1);
	
	print " >> Alignment #$aln\n";
	
	# look for human ENSP IDs in msffile
	@ensps = `grep -E "^ENSP[0-9]" $msffile | perl -ane 'print "\$F[0]\n";'` or next;
	# shuffle(@a);
	
	foreach $ensp (unique(@ensps))
	{
		chomp($ensp);
		$ensp =~ s/^>//;
		
		print "   >> Reference sequence $ensp\n";

		# $query = Query("SELECT id FROM `evorate` WHERE name='$ensp' AND type='$type'");
		# if (Numrows($query) > 0)
		# {
		# 	print "     >> Already have '$type' values for '$ensp' in table 'evorate'! Skipping (use -clear switch on main.pl or on submit.pl to clear old values)\n";
		# 	next;
		# }

		@a = `grep $aln.fasta.txt ../output-check.txt` or die;
		die if (scalar(@a) != 1);
		chomp($a[0]);
		@a = split(/\t/, $a[0]);
		$nseq = $a[1];
		$lseq = $a[3];
		$score = $nseq * $lseq;
		# $expect = $score * $factor;
		print "     >> $nseq sequences * $lseq alignment length = $score\n";
		
		if ($nseq < $nseqmin)
		{
			print "       >> too few sequences (<$nseqmin), skipping!\n";
			next;
		}
		
		if (!switch('all'))
		{
    		$mod = 0;
    		$query = Query("SELECT name FROM uniens WHERE ensp='$ensp' AND species='$species'");
    		if (Numrows($query) == 0)
    		{
    			addme("no names in uniens for ensp", $ensp);
    			next;
    		}
    		$query = Query("SELECT m.id FROM unimod m, uniens e WHERE e.ensp='$ensp' AND e.name=m.name AND e.species='$species' AND m.ptm!=''");
    		if (Numrows($query) > 0)
    		{
    			$mod = 1;
    		}

    		if ($mod == 0)
    		{
    			if (switch('all'))
    			{
    				print "     >> Doesn't have any PTMs according to unimod 'ptm' field, but switch -all is enabled\n";
    			}
    			else
    			{
    				print "     >> Doesn't have any PTMs according to unimod 'ptm' field (use -all to submit these anyway)\n";
    				next;
    			}
    		}
    		else
    		{
    			print "     >> Has PTMs\n";
    		}
		}

        # # Get Compara sequence
        # $query = Query("SELECT REPLACE(seq, '-', '') FROM comparafasta WHERE ensp='$ensp'");
        # if (Numrows($query) > 1)
        # {
        #   die("Error: Multiple sequences in Compara for ENSP '$ensp'");
        #   # addme("multiple sequences in compara for ensp", $ensp);
        #   # next;
        # }
        # if (Numrows($query) == 0)
        # {
        #   die("Error: No sequence in Compara for ENSP '$ensp'");
        #   # addme("no sequence in compara for ensp", $ensp);
        #   next;
        # }
        # ($comparaseq) = FetchOne($query);
        # 
        # # Try to get matching UniProt sequence
        # $query = Query("SELECT s.seq FROM uniseq s, uniens e WHERE e.ensp='$ensp' AND e.name=s.name AND e.species='$species' AND s.type='UniProt' AND s.seq='$comparaseq'");
        # if (Numrows($query) == 0)
        # {
        #   addme("couldn't find a uniprot name with matching sequence for ensp", $ensp);
        #   next;
        # }

		$treefile = $msffile;
		$treefile =~ s/\.msf\./.tree./;

		# $logfile = "../output/log-$aln-$ensp.txt";

		print "       >> Running lichtarge...\n";

		# starttime();
		while ($free <= 0)
		{
			sleep(5);
			$free = freenodes(1);
		}

		run("qsub", "~/scripts/qsub.sh ../ETcode/wetc -p $msffile -x $ensp -readtree $treefile -o ../output/$ensp", 1);
		
		addme("ran ensps", $ensp);
		
		# print "       >> ";
		# stoptime(1);
		
		# $success = 0;
		# @tmp = `ls -1 output/$ensp*`;
		# if (scalar(@tmp) > 0)
		# {
		# 	$success = 1;
		# 	print "         >> Success! (output exists)\n";
		# }		
		# else
		# {
		# 	print "         >> Lichtarge crashed! (no output)\n";
		# }
		
		# if (`grep -c "likelihoodComputation::getLofPos: likelihood of pos was zero!" $logfile` > 0)
		# {
		# 	print "         >> rate4site crashed! ('likelihood of pos was zero')\n";
		# 	next;
		# }
		# 
		# if (`grep -c "Amino acid was not one of the following" $logfile` > 0)
		# {
		# 	print "         >> rate4site crashed! ('Amino acid was not one of the following')\n";
		# 	next;
		# }



		# print "         >> Logging to '$outfile'\n";
		# 
		# print OUT "$aln\t$nseq\t$lseq\t$score\t$ensp\t$success\n";
	}
}
chdir("..");
stopme();
stoptime();

showme("ran ensps", 1);
print "\n";
# showme("no names in uniens for ensp", 1);
# showme("couldn't find a uniprot name with matching sequence for ensp", 1);

done();