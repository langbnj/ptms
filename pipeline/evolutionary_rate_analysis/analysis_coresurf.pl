#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species] [source: UniProt/Ochoa/PhosphoSitePlus/dbPTM/all] [evorate type] [disorder predictor] [minimum number of sites for a PTM type] [-nolog]\n\n".
"\n -nolog: Passes on -nolog argument to draw.pl (draw figures with linear scale, values are unaffected)".
"\nExample: $0 human UniProt lichtarge_linsi_tree MULTICOM 1000";
($species, $source, $evorate, $predictor, $minsites) = args(5);

# $infile = "input.txt";
$outfile = "tmp/tmp-dataframe-all-$species-$source-$evorate-$predictor-coresurf.txt";

open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
print OUT "ptm\tacc\tsite\ttype\tdis\trate\n";

# PTM source
$tmpsource = '';
$tmpsourcem = '';
if ($source ne 'all')
{
	$tmpsource = " AND source='$source'";
	$tmpsourcem = " AND m.source='$source'";
}

$tmp_plddt = '';
if ($predictor =~ /_pLDDT70$/)
{
	$tmp_plddt = '_pLDDT70';
}
elsif ($predictor =~ /_pLDDT90$/)
{
	$tmp_plddt = '_pLDDT90';
}
# Note - even though this doesn't explicitly handle _PAE1 and _PAE2, $predictor will lead to correct filtering).




# start
$ptmquery = Query("SELECT ptm FROM (SELECT ptm, COUNT(DISTINCT acc, site) AS c FROM unimod WHERE species='$species'$tmpsource AND ptm IS NOT NULL GROUP BY ptm ORDER BY ptm) AS t WHERE c >= $minsites");
startme2();
startme3();
while (($ptm) = Fetch($ptmquery))
{
	# clearme("no control sites for acc");

	# Control
	startme("Looking at '$ptm' control sites");
	starttime();
	$proteinquery = Query("SELECT 
	    m.acc
	FROM
	    unimod m
	WHERE
	    m.species = '$species'$tmpsourcem AND m.ptm = '$ptm'
	GROUP BY m.acc");
	# while (($acc, $disseq, $seq, $count) = Fetch($proteinquery))
	while (($acc) = Fetch($proteinquery))
	{
		# evorate must exist for at least one PTM site
		$query = Query("SELECT 
		    m.acc
		FROM
		    unimod m,
		    uniens e,
		    `evorate_$evorate` r
		WHERE
		    m.acc = '$acc' AND m.species = '$species'$tmpsourcem AND m.ptm = '$ptm' AND e.acc = m.acc AND e.ensp = r.ensp AND m.site = r.site
		GROUP BY m.acc");
		if (Numrows($query) == 0)
		{
			addme("PTM: no evorate information for acc", $acc) if switch('debug');
			next;
		}
		
		# disorder must exist for at least one PTM site
		$query = Query("SELECT 
		    m.acc
		FROM
		    unimod m,
		    uniseq s
		WHERE
		    m.acc = '$acc' AND m.acc = s.acc AND m.species = '$species'$tmpsourcem AND m.ptm = '$ptm' AND s.type = '$predictor' AND SUBSTRING(s.seq, m.site, 1) IN ('*','.')
		GROUP BY m.acc");
		if (Numrows($query) == 0)
		{
			addme("PTM: no disorder information for acc", $acc) if switch('debug');
			next;
		}
		
		# coresurf must exist for at least one PTM site
		$query = Query("SELECT 
		    m.acc
		FROM
		    unimod m,
		    uniseq s
		WHERE
		    m.acc = '$acc' AND m.acc = s.acc AND m.species = '$species'$tmpsourcem AND m.ptm = '$ptm' AND s.type = 'CoreSurf$tmp_plddt' AND SUBSTRING(s.seq, m.site, 1) IN ('C','S')
		GROUP BY m.acc");
		if (Numrows($query) == 0)
		{
			addme("PTM: no disorder information for acc", $acc) if switch('debug');
			next;
		}
		
		# get disseq and surfseq and seq for protein
		$query = Query("SELECT 
		    s.seq, t.seq, u.seq
		FROM
		    unimod m,
		    uniseq s,
		    uniseq t,
		    uniseq u
		WHERE
		    m.acc = '$acc' AND m.acc = s.acc AND m.acc = t.acc AND m.acc = u.acc AND m.species = '$species'$tmpsourcem AND m.ptm = '$ptm' AND s.type = '$predictor' AND t.type = 'CoreSurf$tmp_plddt' AND u.type = 'UniProt' AND SUBSTRING(s.seq, m.site, 1) IN ('*','.') AND SUBSTRING(t.seq, m.site, 1) IN ('C','S')  AND SUBSTRING(u.seq, m.site, 1) = '".substr($ptm, 0, 1)."'
		GROUP BY m.acc, s.seq, t.seq, u.seq");
		if (Numrows($query) == 0)
		{
			addme("PTM: no disorder information for acc", $acc) if switch('debug');
			next;
		}
		($disseq, $surfseq, $seq) = FetchOne($query);
		
		
		
		# get $count (number of PTM sites with evorate and disorder info)
		$count = 0;
		$countstr = 0;
		$countstrcore = 0;
		$countstrsurf = 0;
		$countdis = 0;
		$countdiscore = 0;
		$countdissurf = 0;
		$sitequery = Query("SELECT 
		    m.site
		FROM
		    unimod m
		WHERE
		    m.acc = '$acc' AND m.species = '$species'$tmpsourcem AND m.ptm = '$ptm'
		GROUP BY m.site");
		while (($site) = Fetch($sitequery))
		{
			# Get disorder information
			$query = Query("SELECT SUBSTRING(seq, $site, 1) FROM uniseq WHERE acc='$acc' AND type='$predictor' AND SUBSTRING(seq, $site, 1) IN ('*','.')");
			if (Numrows($query) == 0)
			{
				addme("Control/PTM: no disorder information for acc", $acc) if switch('debug');
				addme("Control/PTM: no disorder information for acc|site", "$acc|$site") if switch('debug');
				next;
			}
			($dis) = FetchOne($query);

			# Get core/surface information
			$query = Query("SELECT SUBSTRING(seq, $site, 1) FROM uniseq WHERE acc='$acc' AND type='CoreSurf$tmp_plddt' AND SUBSTRING(seq, $site, 1) IN ('C','S')");
			if (Numrows($query) == 0)
			{
				addme("Control/PTM: no surface/core information for acc", $acc) if switch('debug');
				addme("Control/PTM: no surface/core information for acc|site", "$acc|$site") if switch('debug');
				next;
			}
			($surf) = FetchOne($query);

			# Get evorate information
			$query = Query("SELECT DISTINCT r.rate FROM uniens e, `evorate_$evorate` r WHERE e.acc='$acc' AND e.ensp=r.ensp AND r.site=$site GROUP BY r.rate");
			if (Numrows($query) == 0)
			{
				addme("Control/PTM: no evorate information for acc", $acc) if switch('debug');
				addme("Control/PTM: no evorate information for acc|site", "$acc|$site") if switch('debug');
				next;
			}
			if (Numrows($query) > 1)
			{
				addme("Control/PTM: multiple evorates for acc (kept, using minimum as one gene copy can compensate for the other(s))", $acc) if switch('debug');
				addme("Control/PTM: multiple evorates for acc|site (kept, using minimum as one gene copy can compensate for the other(s))", "$acc|$site") if switch('debug');
				# next;
			}

			# Count PTM sites
			if ($dis eq '*')
			{
				$count++;
				$countdis++;
				$countdissurf++ if ($surf eq 'S');
				$countdiscore++ if ($surf eq 'C');
			}
			elsif ($dis eq '.')
			{
				$count++;
				$countstr++;
				$countstrsurf++ if ($surf eq 'S');
				$countstrcore++ if ($surf eq 'C');
			}
			else
			{
				# blank (e.g. in consensus)
				die;
			}
		}
		

		# always '-zeroptms' now: Control sites can't have any PTM sites annotated, not even 'by similarity' etc.
		# $evoquery = Query("SELECT 
		#     r.site, r.rate, COUNT(DISTINCT r.rate), SUBSTRING(s.seq, r.site, 1), SUBSTRING(t.seq, r.site, 1)
		$evoquery = Query("SELECT 
		    r.site, MIN(r.rate), COUNT(DISTINCT r.rate), SUBSTRING(s.seq, r.site, 1), SUBSTRING(t.seq, r.site, 1)
		FROM
		    `evorate_$evorate` r,
		    uniens u,
		    uniseq s,
		    uniseq t
		WHERE
		    u.acc = '$acc' AND u.ensp = r.ensp AND r.site IN ('".join("','", positions(substr($ptm, 0, 1), $seq))."') AND r.site NOT IN (SELECT 
		        m.site
		    FROM
		        unimod m
		    WHERE
		        m.acc = '$acc' AND m.species = '$species'$tmpsourcem) AND s.acc = u.acc AND t.acc = u.acc AND s.type = '$predictor' AND t.type = 'CoreSurf$tmp_plddt' AND SUBSTRING(s.seq, r.site, 1) IN ('*','.') AND SUBSTRING(t.seq, r.site, 1) IN ('C','S')
		GROUP BY r.site");
		%rates = ();
		%ratesdis = ();
		%ratesdissurf = ();
		%ratesdiscore = ();
		%ratesstr = ();
		%ratesstrsurf = ();
		%ratesstrcore = ();
		while (($site, $rate, $ratecount, $dis, $surf) = Fetch($evoquery))
		{
			die if ($ratecount == 0);	# Seems pointless since this is by definition > 0 if something is being fetched
			if ($ratecount > 1)
			{
				# addme("PTM: multiple evorates for acc", $acc) if switch('debug');
				# addme("PTM: multiple evorates for acc|site", "$acc|$site") if switch('debug');
				addme("PTM: multiple evorates for acc (kept, using minimum as one gene copy can compensate for the other(s))", $acc) if switch('debug');
				addme("PTM: multiple evorates for acc|site (kept, using minimum as one gene copy can compensate for the other(s))", "$acc|$site") if switch('debug');
				# next;
			}
			addme("Control: Total available acc|sites", "$acc|$site") if switch('debug');
			$rates{$site} = $rate;
			if ($dis eq '*') { $ratesdis{$site} = $rate; if ($surf eq 'S') { $ratesdissurf{$site} = $rate; } if ($surf eq 'C') { $ratesdiscore{$site} = $rate; } }
			if ($dis eq '.') { $ratesstr{$site} = $rate; if ($surf eq 'S') { $ratesstrsurf{$site} = $rate; } if ($surf eq 'C') { $ratesstrcore{$site} = $rate; } }
		}

		# all
		foreach $site (keys(%rates))
		{
			$rate = $rates{$site};
		
			if (substr($seq, $site-1, 1) ne substr($ptm, 0, 1))
			{
				die("Error: '$acc' site '$site' is '".substr($seq, $site-1, 1)."' instead of ".substr($ptm, 0, 1));		# Throw error if this AA is not what I expected for this PTM
			}
			
			# Write Control evorates
			print OUT "$ptm\t$acc\t$site\tC\tA\t$rate\n";
			addme("Control: Control accs", $acc) if switch('debug');
			addme("Control: Control acc|sites", "$acc|$site") if switch('debug');
			stepme(100);
		}
		#discore
		if ($countdiscore > 0)
		{
			foreach $site (keys(%ratesdiscore))
			{
				$rate = $ratesdiscore{$site};

				if (substr($seq, $site-1, 1) ne substr($ptm, 0, 1))
				{
					die("Error: '$acc' site '$site' is '".substr($seq, $site-1, 1)."' instead of ".substr($ptm, 0, 1));		# Throw error if this AA is not what I expected for this PTM
				}

				# Write Control evorates
				if ((substr($disseq, $site-1, 1) eq '*') and (substr($surfseq, $site-1, 1) eq 'C'))
				{
					print OUT "$ptm\t$acc\t$site\tC\tdiscore\t$rate\n";

					addme("Control: Disordered Core Control accs", $acc) if switch('debug');
					addme("Control: Disordered Core Control acc|sites", "$acc|$site") if switch('debug');
				}
			}
		}
		#dissurf
		if ($countdissurf > 0)
		{
			foreach $site (keys(%ratesdissurf))
			{
				$rate = $ratesdissurf{$site};
		
				if (substr($seq, $site-1, 1) ne substr($ptm, 0, 1))
				{
					die("Error: '$acc' site '$site' is '".substr($seq, $site-1, 1)."' instead of ".substr($ptm, 0, 1));		# Throw error if this AA is not what I expected for this PTM
				}
			
				# Write Control evorates
				if ((substr($disseq, $site-1, 1) eq '*') and (substr($surfseq, $site-1, 1) eq 'S'))
				{
					print OUT "$ptm\t$acc\t$site\tC\tdissurf\t$rate\n";

					addme("Control: Disordered Surface Control accs", $acc) if switch('debug');
					addme("Control: Disordered Surface Control acc|sites", "$acc|$site") if switch('debug');
				}
			}
		}
		# strcore
		if ($countstrcore > 0)
		{
			foreach $site (keys(%ratesstrcore))
			{
				$rate = $ratesstrcore{$site};

				if (substr($seq, $site-1, 1) ne substr($ptm, 0, 1))
				{
					die("Error: '$acc' site '$site' is '".substr($seq, $site-1, 1)."' instead of ".substr($ptm, 0, 1));		# Throw error if this AA is not what I expected for this PTM
				}

				# Write Control evorates
				if ((substr($disseq, $site-1, 1) eq '.') and (substr($surfseq, $site-1, 1) eq 'C'))
				{
					print OUT "$ptm\t$acc\t$site\tC\tstrcore\t$rate\n";

					addme("Control: Structured Core Control accs", $acc) if switch('debug');
					addme("Control: Structured Core Control acc|sites", "$acc|$site") if switch('debug');
				}
			}
		}
		# strsurf
		if ($countstrsurf > 0)
		{
			foreach $site (keys(%ratesstrsurf))
			{
				$rate = $ratesstrsurf{$site};

				if (substr($seq, $site-1, 1) ne substr($ptm, 0, 1))
				{
					die("Error: '$acc' site '$site' is '".substr($seq, $site-1, 1)."' instead of ".substr($ptm, 0, 1));		# Throw error if this AA is not what I expected for this PTM
				}

				# Write Control evorates
				if ((substr($disseq, $site-1, 1) eq '.') and (substr($surfseq, $site-1, 1) eq 'S'))
				{
					print OUT "$ptm\t$acc\t$site\tC\tstrsurf\t$rate\n";

					addme("Control: Structured Surface Control accs", $acc) if switch('debug');
					addme("Control: Structured Surface Control acc|sites", "$acc|$site") if switch('debug');
				}
			}
		}
		
		# stepme(100);
	}
	stopme();
	# # uniqueme("Control: no control sites for ptm|acc");
	# uniqueme("no control sites for acc");
	stoptime();

	# PTM sites
	startme("Looking at '$ptm' PTM sites");
	starttime();
	$sitequery = Query("SELECT 
	    m.acc, m.site
	FROM
	    unimod m
	WHERE
	    m.species = '$species'$tmpsourcem AND m.ptm = '$ptm'
	GROUP BY m.acc,m.site;");
	while (($acc, $site) = Fetch($sitequery))
	{
		# if (contains($acc, returnme("no control sites for acc", 1)))
		# {
		# 	addme("PTM: ptm site skipped because Control: no control sites for ptm|acc for acc", $acc) if switch('debug');
		# 	addme("PTM: ptm site skipped because Control: no control sites for ptm|acc for acc|site", "$acc|$site") if switch('debug');
		# 	next;
		# }
		
		# Note: the "next"s below cause control sites to be in the output without corresponding PTM sites (since the control sites get written first, above). draw.pl takes care of this by filtering out these cases, but it's a bit confusing.

		# Get disorder information
		$query = Query("SELECT SUBSTRING(seq, $site, 1) FROM uniseq WHERE acc='$acc' AND type='$predictor' AND SUBSTRING(seq, $site, 1) IN ('*','.')");
		if (Numrows($query) == 0)
		{
			addme("PTM: no disorder information for acc", $acc) if switch('debug');
			addme("PTM: no disorder information for acc|site", "$acc|$site") if switch('debug');
			next;
		}
		($dis) = FetchOne($query);

		# Get surface/core information
		$query = Query("SELECT SUBSTRING(seq, $site, 1) FROM uniseq WHERE acc='$acc' AND type='CoreSurf$tmp_plddt' AND SUBSTRING(seq, $site, 1) IN ('C','S')");
		if (Numrows($query) == 0)
		{
			addme("PTM: no surface/core information for acc", $acc) if switch('debug');
			addme("PTM: no surface/core information for acc|site", "$acc|$site") if switch('debug');
			next;
		}
		($surf) = FetchOne($query);

		# Get evorate information
		$query = Query("SELECT DISTINCT r.rate FROM uniens e, `evorate_$evorate` r WHERE e.acc='$acc' AND e.ensp=r.ensp AND r.site=$site GROUP BY r.rate");
		if (Numrows($query) == 0)
		{
			addme("PTM: no evorate information for acc", $acc) if switch('debug');
			addme("PTM: no evorate information for acc|site", "$acc|$site") if switch('debug');
			next;
		}
		if (Numrows($query) > 1)
		{
			$query = Query("SELECT MIN(r.rate) FROM uniens e, `evorate_$evorate` r WHERE e.acc='$acc' AND e.ensp=r.ensp AND r.site=$site");
			addme("PTM: multiple evorates for acc (kept, using minimum as one gene copy can compensate for the other(s))", $acc) if switch('debug');
			addme("PTM: multiple evorates for acc|site (kept, using minimum as one gene copy can compensate for the other(s))", "$acc|$site") if switch('debug');
			# next;
		}
		($rate) = FetchOne($query);

		# print "$acc\t$site\t$dis\t$aa\t$ensp\t$rate\n";
	
		# Write PTM site evorates
		if (($dis eq '*') and ($surf eq 'S'))
		{
			print OUT "$ptm\t$acc\t$site\tP\tA\t$rate\n";
			print OUT "$ptm\t$acc\t$site\tP\tdissurf\t$rate\n";

			addme("PTM: accs", $acc) if switch('debug');
			addme("PTM: acc|sites", "$acc|$site") if switch('debug');
			addme("PTM: ptm|acc|sites", "$ptm|$acc|$site") if switch('debug');
		}
		elsif (($dis eq '*') and ($surf eq 'C'))
		{
			print OUT "$ptm\t$acc\t$site\tP\tA\t$rate\n";
			print OUT "$ptm\t$acc\t$site\tP\tdiscore\t$rate\n";
			
			addme("PTM: accs", $acc) if switch('debug');
			addme("PTM: acc|sites", "$acc|$site") if switch('debug');
			addme("PTM: ptm|acc|sites", "$ptm|$acc|$site") if switch('debug');
		}
		elsif (($dis eq '.') and ($surf eq 'S'))
		{
			print OUT "$ptm\t$acc\t$site\tP\tA\t$rate\n";
			print OUT "$ptm\t$acc\t$site\tP\tstrsurf\t$rate\n";
			
			addme("PTM: accs", $acc) if switch('debug');
			addme("PTM: acc|sites", "$acc|$site") if switch('debug');
			addme("PTM: ptm|acc|sites", "$ptm|$acc|$site") if switch('debug');
		}
		elsif (($dis eq '.') and ($surf eq 'C'))
		{
			print OUT "$ptm\t$acc\t$site\tP\tA\t$rate\n";
			print OUT "$ptm\t$acc\t$site\tP\tstrcore\t$rate\n";
			
			addme("PTM: accs", $acc) if switch('debug');
			addme("PTM: acc|sites", "$acc|$site") if switch('debug');
			addme("PTM: ptm|acc|sites", "$ptm|$acc|$site") if switch('debug');
		}
		else
		{
			# blank (e.g. in consensus)
			addme("PTM: no disorder or surface/core information for acc", $acc) if switch('debug');
			addme("PTM: no disorder or surface/core information for acc|site", "$acc|$site") if switch('debug');
			# next;
			die;
		}
		
		stepme(100);
	}
	stopme();
	stoptime();

	if (switch('debug'))
	{
		showmeall(1);
		clearall();
	}

	# last;
}

close(OUT);

# clearme("no control sites for acc");
showmeall(1);

state("Wrote to '$outfile'");

done();

# if (!switch('debug'))
# {
# 	if (switch('nolog'))
# 	{
# 		run("Draw plots", "draw.pl $species $evorate $predictor all -nolog");
# 		run("Draw plots", "draw.pl $species $evorate $predictor all -dismerge -nolog");
# 		run("Draw plots", "draw.pl $species $evorate $predictor all -surfmerge -nolog");
# 		# run("Draw plots", "draw.pl $species $evorate $predictor all -nolog -bonferroni");
# 	}
# 	else
# 	{
# 		run("Draw plots", "draw.pl $species $evorate $predictor all");
# 		run("Draw plots", "draw.pl $species $evorate $predictor all -dismerge");
# 		run("Draw plots", "draw.pl $species $evorate $predictor all -surfmerge");
# 		# run("Draw plots", "draw.pl $species $evorate $predictor all -bonferroni");
# 	}
# }
