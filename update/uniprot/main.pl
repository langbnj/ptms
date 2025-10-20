#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

require XML::Twig;
require IO::Uncompress::Gunzip;

# initialize

# $superloudmysql = 1;
$varsplicfile = "input/uniprot_sprot_varsplic.fasta";


our $usage = "$0 [chunk file name] [-humanonly]\n\n -humanonly: Only retain proteins that either are human, or that are from organisms with human hosts.\n\nExample: $0 input/uniprot_sprot_0001.xml";
($infile) = args(1);

$sizefile = $infile.".size.txt";
if (!-s $sizefile)
{
	if ($infile =~ /\.gz$/)
	{
		# Compressed XML file (.xml.gz)
		run("Get input file size", "zcat '$infile' | grep -c '<entry' > $sizefile", 1);
	}
	else 
	{
		# Uncompressed XML file (.xml)
		run("Get input file size", "cat '$infile' | grep -c '<entry' > $sizefile", 1);
	}
}
$size = chompme(`cat $sizefile`);

# Load varsplic isoform sequences (since these won't be in any table yet)
startme("Reading isoform sequences from '$varsplicfile'");
open(IN, $varsplicfile) or die("Error: Couldn't open '$varsplicfile");
fastabreak();
%varsplic = ();
%varsplic_trembl = ();
while (<IN>)
{
	($title, $seq) = getfasta();

	# >sp|P48347-2|14310_ARATH Isoform 2 of 14-3-3-like protein GF14 epsilon OS=Arabidopsis thaliana OX=3702 GN=GRF10
	# MENEREKQVYLAKLSEQTERYDEMVEAMKKVAQLDVELTVEERNLVSVGYKNVIGARRAS
	# WRILSSIEQKEESKGNDENVKRLKNYRKRVEDELAKVCNDILSVIDKHLIPSSNAVESTV
	# FFYKMKGDYYRYLAEFSSGAERKEAADQSLEAYKAAVAAAENGLAPTHPVRLGLALNFSV
	# FYYEILNSPESACQLAKQAFDDAIAELDSLNEESYKDSTLIMQLLRDNLTLWTSDLNEEG
	# DERTKGADEPQDEV
	
	$title =~ /^(sp|tr)\|([^\|]+)\|/ or die("Error: Couldn't match isoform acc in title '$title' in '$varsplicfile'");
	$trembl = $1;
	$acc = $2;
	
	$varsplic{$acc} = $seq;
	
	if ($trembl eq 'sp')
	{
		$trembl = 0;
	}
	elsif ($trembl eq 'tr')
	{
		$trembl = 1;
	}
	else
	{
		die("Error: Unhandled value '$trembl' for sp|tr in title '$title'");
	}
	$varsplic_trembl{$acc} = $trembl;
	
	
	stepme(10000);
}
stopme();
normalbreak();
close(IN);



$xml = XML::Twig -> new
(
	TwigRoots => {entry => 1},
	twig_handlers => { entry => \&uniprot_entry },
);

# @allowedspecies = ();
# push(@allowedspecies, 'arath');
# push(@allowedspecies, 'ecoli');
# push(@allowedspecies, 'human');
# push(@allowedspecies, 'mouse');
# push(@allowedspecies, 'salty');
# push(@allowedspecies, 'yeast');

# start

if ($infile =~ /\.gz$/)
{
	# Compressed XML file (.xml.gz)
	startme("Parsing compressed XML file '$infile'", 0, $size);
	starttime();
	$gzipfile = IO::Uncompress::Gunzip->new($infile) or die("Error: Couldn't gunzip '$infile'");
	$xml -> parse($gzipfile);
}
else
{
	# Uncompressed XML file (.xml)
	startme("Parsing XML file '$infile'", 0, $size);
	starttime();
	$xml -> parsefile($infile);
}


sub uniprot_entry
{
	my($t, $x) = @_;      # arguments for all twig_handlers
	



	# Parse UniProt entry
	
	# UniProt name
    die("Error: Duplicates for UniProt name") if ($x -> children_count('name') != 1);
    my $name = $x -> children_text('name');
	# print " >> $name\n" if (switch('debug'));

	# Species (UniProt species, taken from name)
	$name =~ /\w+_(\w+)/ or die("Error: Couldn't parse UniProt species from UniProt name '$name'");
	my $species = $1;
	
	# Get primary accession
    my $primary_acc = $x -> first_child('accession') -> text;
	
	# Dataset: Swiss-Prot or TrEMBL
	# <entry dataset="Swiss-Prot" created="1993-07-01" modified="2021-02-10" version="232" xmlns="http://uniprot.org/uniprot">
    my $dataset = $x->{'att'}->{'dataset'};
	my $trembl = 0;
	if ($dataset eq 'Swiss-Prot')
	{
		$trembl = 0;
	}
	elsif ($dataset eq 'TrEMBL')
	{
		$trembl = 1;
	}
	else
	{
		die("Error: Unhandled dataset '$dataset'");
	}
	
	
	
	
	# Check if this is a human protein
	my $human_protein = 0;
	if ($species eq 'HUMAN')
	{
		$human_protein = 1;
	}
	
	# Check if this is a viral protein with a human host (bacteria don't currently have host information)
	# (Abbreviated version of unihost parsing below)
	my $human_host = 0;
    foreach $organismhost ($x -> children('organismHost'))
	{
		$dbreference = $organismhost->first_child('dbReference');
		if ($dbreference->{'att'}->{'type'} eq 'NCBI Taxonomy')
		{
			# my $type = $dbreference->{'att'}->{'type'};
			my $hosttax = $dbreference->{'att'}->{'id'};
			my $host = '';
			
			# Skip non-human-hosted proteins if switch -humanonly is active
			next if (switch('humanonly') and ($hosttax ne '9606'));
			$human_host = 1 if ($hosttax eq '9606');
		}
	}

	# Skip proteins that are not human and that aren't viral proteins with a human host if switch -humanonly is active
	if (switch('humanonly'))
	{
		if (($human_protein == 0) and ($human_host == 0))
		{
			$x -> purge;              # frees memory
			stepme(100);
			return;
		}
	}

	
	
	
	# Full name
	# Example:
	# <protein>
	#   <recommendedName ref="1">  
	#     <fullName>Putative pyrroline-5-carboxylate reductase</fullName>
	#     <shortName>P5C reductase</shortName>
	#     <shortName>P5CR</shortName>
	#   </recommendedName>
	# </protein>
	my $fullname = '';
	if (defined($x->first_child('protein')->first_child('recommendedName')))
	{
		# fullname = recommended name (available for all of Swiss-Prot, but not TrEMBL)
	    $fullname = $x->first_child('protein')->first_child('recommendedName')->first_child('fullName')->text;
	}
	elsif (defined($x->first_child('protein')->first_child('submittedName')))
	{
		# fullname = submitted name
	    $fullname = $x->first_child('protein')->first_child('submittedName')->first_child('fullName')->text;
	}
	else
	{
		# none defined
	    $fullname = '';
	}
	
	# Primary gene / ordered locus id / orf id
	# <gene>
	#   <name type="primary">PAU4</name>
	#   <name type="ordered locus">YLR461W</name>
	#   <name type="ORF">L9122.1</name>
	# </gene>
    my $symbol = '';
    my $orderedlocus = '';
    my $orf = '';
	if ($x->first_child('gene'))
	{
		foreach $genename ($x->first_child('gene')->children('name'))
		{
			if ($genename->{'att'}->{'type'} eq 'primary')
			{
				$symbol = $genename->text;
			}
			if ($genename->{'att'}->{'type'} eq 'ordered locus')
			{
				$orderedlocus = $genename->text;
			}
			if ($genename->{'att'}->{'type'} eq 'ORF')
			{
				$orf = $genename->text;
			}
		}
	}

	# NCBI Taxonomy (NCBI Taxon IDs)
	# Example:
	# 
	# <organism>
	#   <name type="scientific">Pasteurella multocida (strain Pm70)</name>
	#   <dbReference id="272843" type="NCBI Taxonomy"/>
	#   <lineage>
	#     <taxon>Bacteria</taxon>
	#     <taxon>Proteobacteria</taxon>
	#     <taxon>Gammaproteobacteria</taxon>
	#     <taxon>Pasteurellales</taxon>
	#     <taxon>Pasteurellaceae</taxon>
	#     <taxon>Pasteurella</taxon>
	#   </lineage>
	# </organism>
	# if ($x->first_child('organism'))
	foreach $organism ($x->children('organism'))
	{
		# $dbreference = $x->first_child('organism')->first_child('dbReference');
		foreach $dbreference ($organism->children('dbReference'))
		{
			if ($dbreference->{'att'}->{'type'} eq 'NCBI Taxonomy')
			{
				# my $type = $dbreference->{'att'}->{'type'};
				my $tax = $dbreference->{'att'}->{'id'};
				
				# Update unitax, if this taxon ID isn't in there yet
				# Query("LOCK TABLE unitax WRITE") unless switch('debug');		# to prevent duplicate inserts due to concurrent access (faster alternative would be to just say "DISTINCT species" below)		# Not needed anymore since I'm using a single job now
				$query = Query("SELECT species FROM unitax WHERE tax='$tax'");
				if (Numrows($query) == 0)
				{
					Query("INSERT INTO unitax SET tax='$tax', species='$species'") unless switch('debug2');
				}
				else
				{
					($thisspec) = FetchOne($query);
					die("Error: Taxon ID '$tax' already exists in unitax, but it's linked to '$thisspec' instead of '$species'") if ($thisspec ne $species);
				}
				# Query("UNLOCK TABLES") unless switch('debug');
			}
		}
	}
	
	# Host organisms (for viruses)
	# Example:
	# 
    # <organismHost>
    #   <name type="scientific">Homo sapiens</name>
    #   <name type="common">Human</name>
    #   <dbReference type="NCBI Taxonomy" id="9606"/>
    # </organismHost>
    # <organismHost>
    #   <name type="scientific">Macaca</name>
    #   <name type="common">macaques</name>
    #   <dbReference type="NCBI Taxonomy" id="9539"/>
    # </organismHost>
	# my $human_host = 0;
    foreach $organismhost ($x -> children('organismHost'))
	{
		$dbreference = $organismhost->first_child('dbReference');
		if ($dbreference->{'att'}->{'type'} eq 'NCBI Taxonomy')
		{
			# my $type = $dbreference->{'att'}->{'type'};
			my $hosttax = $dbreference->{'att'}->{'id'};
			my $host = '';
			
			# # Skip non-human-hosted proteins if switch -humanonly is active
			# next if (switch('humanonly') and ($hosttax ne '9606'));
			# $human_host = 1 if ($hosttax eq '9606');
			
			# Get host species mnemonic if it exists
			$query = Query("SELECT species FROM unitax WHERE tax='$hosttax'");
			if (Numrows($query) == 0)
			{
				# addme("no species mnemonic in table unitax for hosttax (kept)", $hosttax);
				
				# Filter out viruses that do not have a host species that's in unitax (i.e. human)
				addme("no species mnemonic in table unitax for hosttax (kept)", $hosttax);
				# next;
				$host = '';
			}
			else
			{
				($host) = FetchOne($query);
			}
			
			$q = "INSERT INTO unihost SET name='$name', acc='$primary_acc', species='$species', hosttax='$hosttax', host='$host'";
			$q =~ s/=''/=NULL/g;
			Query($q) unless switch('debug');
		}
	}

	# Accessions
    my @acc = $x -> children_text('accession');
    die("Error: UniProt accessions present multiple times for entry '$name'") if (scalar(@acc) != scalar(unique(@acc)));
	$primary_acc = $acc[0];
	foreach $acc (@acc)
	{
		# print "INSERT INTO uniacc SET name='$name', species='$species', acc='$acc'\n";
		Query("INSERT INTO uniacc SET name='$name', species='$species', acc='$acc', canon='$primary_acc', trembl='$trembl'") unless (switch('debug'));
	}
	
	# Evidence for protein existence
	#   <proteinExistence type="predicted"/>
	
	my $existence = '';
	if (defined($x->first_child('proteinExistence')->{'att'}->{'type'}))
	{
	    $existence = $x->first_child('proteinExistence')->{'att'}->{'type'};
	}
	else
	{
		# none defined
	    $existence = '';
	}
	
	
	# Keywords
	
	# Complete proteome / Reference proteome membership
	# Complete:
	#   <keyword id="KW-0181">Complete proteome</keyword>
	# Reference:
  	#   <keyword id="KW-1185">Reference proteome</keyword>
	
	my $completeproteome = '0';
	my $referenceproteome = '0';
    foreach $keyword ($x -> children('keyword'))
	{
		if ($keyword->{'att'}->{'id'} eq 'KW-0181')
		{
			$completeproteome = '1';
		}
		if ($keyword->{'att'}->{'id'} eq 'KW-1185')
		{
			$referenceproteome = '1';
		}
		
		my $keyid = $keyword->{'att'}->{'id'};
		my $tmpkeyword = $keyword -> text;
		
		# Also insert into table 'unikey'
		Query("INSERT INTO unikey SET name='$name', acc='$primary_acc', species='$species', keyid='$keyid', keyword='$tmpkeyword'") unless switch('debug');
	}
	
	
	# Update uniprot
	# print "INSERT INTO uniprot SET name='$name', species='$species'\n";
	# Query("INSERT INTO uniprot SET name='$name', species='$species', fullname='".esc($fullname)."', symbol='".esc($gene)."', orderedlocus='".esc($orderedlocus)."', orf='".esc($orf)."'") unless switch('debug');
	$q = "INSERT INTO uniprot SET name='$name', acc='$primary_acc', species='$species', fullname='".esc($fullname)."', symbol='".esc($symbol)."', orderedlocus='".esc($orderedlocus)."', orf='".esc($orf)."', evidence='$existence', completeproteome='$completeproteome', referenceproteome='$referenceproteome', trembl='$trembl'";
	$q =~ s/=''/=NULL/g;
	Query($q) unless switch('debug');
	
	# Sequence
  	# <sequence checksum="BA860E2C779BDE52" length="1069" mass="121814" modified="2009-06-16" version="2">MSVIFNASVNTKSAVEYQTISSTQSHSAEEQSERLHKWISKDQLEKLYSSFLNTPERHVGIDELRIILEELDITFNDSMYTRLFLKINQNR</sequence>
	# Complete, but precursor:
	#   <sequence checksum="7D66503ACD865152" length="362" mass="40505" modified="1993-10-01" precursor="true" version="1">
	# Fragment:
	#   <sequence checksum="E999A5A3079D6FCA" fragment="single" length="61" mass="6760" modified="1992-12-01" version="1">
    die("Error: Duplicates for UniProt sequence for '$name'") if ($x -> children_count('sequence') != 1);
	my $seq = $x -> children_text('sequence');
	# Get status (fragments / precursor)
	$fragments = '';
	if (defined($x -> first_child('sequence') -> {'att'} -> {'fragment'}))
	{
		$fragments = $x -> first_child('sequence') -> {'att'} -> {'fragment'};
	}
	$precursor = '';
	if (defined($x -> first_child('sequence') -> {'att'} -> {'precursor'}))
	{
		$precursor = $x -> first_child('sequence') -> {'att'} -> {'precursor'};
	}
	# print "INSERT INTO uniseq SET name='$name', species='$species', type='UniProt', seq='$seq'\n";
	# $q = "INSERT INTO uniseq SET name='$name', species='$species', type='UniProt', fragments='$fragments', precursor='$precursor', trembl='$trembl', seq='$seq'";
	$q = "INSERT INTO uniseq SET name='$name', acc='$primary_acc', canon='$primary_acc', species='$species', type='UniProt', fragments='$fragments', precursor='$precursor', trembl='$trembl', seq='$seq'";
	$q =~ s/=''/=NULL/g;
	Query($q) unless switch('debug');
	
	
	# Get primary isoform's accession (can sometimes have a suffix, e.g. -2 or -3)
	# Needed for uniens_uniprot mapping below (in those dbReference entries, every accession has a suffix, including -1 for the canonical one)
	my $primary_isoform_acc = $primary_acc;
    foreach $comment ($x -> children('comment'))
	{
		if ($comment->{'att'}->{'type'} eq 'alternative products')
		{
			$primary_isoform_acc = $comment -> first_child('isoform') -> first_child('id') -> text;
		}
	}
	if ($primary_isoform_acc ne $primary_acc)
	{
		$tmp_primary_isoform_acc = $primary_isoform_acc;
		# Remove -1 at end
		$tmp_primary_isoform_acc =~ s/-1$//;
		if ($tmp_primary_isoform_acc eq $primary_isoform_acc)
		{
			addme("primary isoform acc differs from primary acc because of -1 for primary acc (kept)", $primary_acc);
			addme("primary isoform acc differs from primary acc because of -1 for primary isoform acc (kept)", $primary_isoform_acc);
			addme("primary isoform acc differs from primary acc because of -1 for primary acc|primary isoform acc (kept)", "$primary_acc|$primary_isoform_acc");
		}
		else
		{
			addme("primary isoform acc differs from primary acc because of something other than -1 for primary acc (kept)", $primary_acc);
			addme("primary isoform acc differs from primary acc because of something other than -1 for primary isoform acc (kept)", $primary_isoform_acc);
			addme("primary isoform acc differs from primary acc because of something other than -1 for primary acc|primary isoform acc (kept)", "$primary_acc|$primary_isoform_acc");
		}
	}
	
	
	
	
	
	# dbReferences
	# 
	# Includes:
	# ENSEMBL: ENSG/ENST/ENSP IDs (for uniens_uniprot)
	# PDB (for unipdb)
	# InterPro (for uniinterpro)
	# GO terms (for unigo)
	# HGNC, MGI, RefSeq (for uniid)
	# Proteome (for uniproteome)
	# 
    foreach $dbreference ($x -> children('dbReference'))
	{
		# print "\ntype=['".$dbreference->{'att'}->{'type'}."']\n\n";
		addme("dbReference types", $dbreference->{'att'}->{'type'});

		# Ensembl (ENSG/ENST/ENSP)
		# Example:
		# Old:
		# # <dbReference id="ENST00000412057" key="35" type="Ensembl">
		# # <property type="protein sequence ID" value="ENSP00000411730"/>
		# # <property type="gene ID" value="ENSG00000231834"/>
		# # </dbReference>
		# New:
		# 	    <dbReference type="Ensembl" id="ENST00000571732">
		# 	      <molecule id="P62258-2"/>
		# 	      <property type="protein sequence ID" value="ENSP00000461762"/>
		# 	      <property type="gene ID" value="ENSG00000108953"/>
		# 	    </dbReference>
		if (contains($dbreference->{'att'}->{'type'}, ('Ensembl', 'EnsemblFungi', 'EnsemblMetazoa')))
		{
			my $enstv = $dbreference->{'att'}->{'id'};
			my $ensgv = '';
			my $enspv = '';
			
			# print "\n >> ENSTV $enstv\n";
		    die("Error: ENSPV/ENSGV duplicates for '$name', ENSTV '$enstv'") if ($dbreference->children_count('property') != 2);
			
			foreach $property ($dbreference->children('property'))
			{
				if ($property->{'att'}->{'type'} eq 'protein sequence ID')
				{
					$enspv = $property->{'att'}->{'value'};
				}
				if ($property->{'att'}->{'type'} eq 'gene ID')
				{
					$ensgv = $property->{'att'}->{'value'};
				}
				# print "   >> Property Type ['".$property->{'att'}->{'type'}."']\n";
				# print "     >> Property Value ['".$property->{'att'}->{'value'}."']\n";
			}
			
			# Get accession (isoform)
			my $ensacc = '';
			if ($dbreference->children_count('molecule') > 0)
			{
			    die("Error: Multiple accs for '$name', ENSTV '$enstv'") if ($dbreference->children_count('molecule') > 1);
				$ensacc = $dbreference->first_child('molecule')->{'att'}->{'id'};
				
				# At least in UniProt 2022_01, Ensembl dbReferences always have an acc suffix (-1 etc.). -1 needs to be removed so that the uniens_uniprot table gets populated correctly.
				# Remove -1 from acc for canonical sequence
				# Canonical sequence isn't always -1
				# $ensacc =~ s/-1$//;
				# If this is the canonical isoform:
				if ($ensacc eq $primary_isoform_acc)
				{
					# Remove -1 (otherwise, keep suffix)
					$ensacc =~ s/-1$//;
				}
			}
			# Fall back on primary accession if there is no "molecule id" for this isoform
			if ($ensacc eq '')
			{
				$ensacc = $primary_acc;
			}
			
			
			# Get unversioned ENSG/ENST/ENSP
			$ensg = $ensgv;
			$ensg =~ s/\.\d+$//;
			$enst = $enstv;
			$enst =~ s/\.\d+$//;
			$ensp = $enspv;
			$ensp =~ s/\.\d+$//;
			
			
			# print "\n >> $name >> $ensg >> $enst >> $ensp\n\n";
			
			# Update uniens_uniprot
			$q = "INSERT INTO uniens_uniprot SET name='$name', acc='$ensacc', canon='$primary_acc', species='$species', ensgv='$ensgv', ensg='$ensg', enstv='$enstv', enst='$enst', enspv='$enspv', ensp='$ensp'";
			$q =~ s/=''/=NULL/g;
			Query($q) unless switch('debug2');
			
			# Sometimes, there's no "molecule id" (acc). I verified that this also happens for proteins that do have multiple isoforms (i.e. I can't just use the canonical primary_acc where acc is NULL):
			# SELECT e.*, COUNT(DISTINCT i.acc) AS isos FROM uniens_uniprot e LEFT OUTER JOIN uniiso i ON i.name=e.name GROUP BY e.id ORDER BY isos DESC;
			# e.g. TAF1_HUMAN
		}
		
		# PDB
		# Example:
		# <dbReference id="1LR5" key="43" type="PDB">
		#   <property type="method" value="X-ray"/>
		#   <property type="resolution" value="1.90 A"/>
		#   <property type="chains" value="A/B/C/D=39-198"/>
		# </dbReference>
		if ($dbreference->{'att'}->{'type'} eq 'PDB')
		{
			my $pdb = lc($dbreference->{'att'}->{'id'});
			my $chains = '';
			my $start = '';
			my $stop = '';
			
			my @chains = ();
			foreach $property ($dbreference->children('property'))
			{
				if ($property->{'att'}->{'type'} eq 'chains')
				{
					push(@chains, $property->{'att'}->{'value'});
				}
			}
			
			if (scalar(@chains) == 0)
			{
				# # die("Error: No chains for '$name' dbReference to PDB ('$pdb')");
				# warn("Warning: No chains for '$name' dbReference to PDB ('$pdb')");

				# No chain annotation, skip this mapping
				# e.g. PGFRA_HUMAN:
			    # <dbReference type="PDB" id="1GQ5">
			    #   <property type="method" value="X-ray"/>
			    #   <property type="resolution" value="2.20 A"/>
			    # </dbReference>
				# Should be:
			    # <dbReference type="PDB" id="6A32">
			    #   <property type="method" value="X-ray"/>
			    #   <property type="resolution" value="1.87 A"/>
			    #   <property type="chains" value="A=550-973"/>
			    # </dbReference>
				addme("no chain information for PDB mapping for name (skipped unipdb entry)", $name);
				addme("no chain information for PDB mapping for name|pdb (skipped unipdb entry)", "$name|$pdb");
				addme("no chain information for PDB mapping for pdb (skipped unipdb entry)", $pdb);
				next;
			}
			die("Error: '$name' PDB chain strings > 1 (pdb '$pdb')") if (scalar(@chains) > 1);
			
			@b = split(/, /, $chains[0]);
			foreach $s (@b)
			{
				# e.g. A/B=-
				if ($s =~ /=-$/)
				{
					# No start/stop coordinates, skip this mapping
					addme("no start/stop coordinates for PDB mapping for name (skipped unipdb entry)", $name);
					addme("no start/stop coordinates for PDB mapping for name|pdb (skipped unipdb entry)", "$name|$pdb");
					addme("no start/stop coordinates for PDB mapping for pdb (skipped unipdb entry)", $pdb);
					next;
				}
				if ($s eq '-')
				{
					# No start/stop coordinates, skip this mapping
					# e.g. RASF1_HUMAN:
				    # <dbReference type="PDB" id="2KZU">
				    #   <property type="method" value="NMR"/>
				    #   <property type="chains" value="B=-"/>
				    # </dbReference>
					addme("no start/stop coordinates for PDB mapping for name (skipped unipdb entry)", $name);
					addme("no start/stop coordinates for PDB mapping for name|pdb (skipped unipdb entry)", "$name|$pdb");
					addme("no start/stop coordinates for PDB mapping for pdb (skipped unipdb entry)", $pdb);
					next;
				}
				my @a = split(/=/, $s);
				die("Error: Chains look weird in '$chains[0]'") if (scalar(@a) != 2);
				$chains = $a[0];
				$chains =~ s/\///g;
				@a = split(/-/, $a[1]);
				if (scalar(@a) != 2)
				{
					# die("Error: Parsing error for '$name' dbReference to PDB ('$pdb'): chains[0] is '".$chains[0]."'");
					warn("Warning: Parsing error for '$name' dbReference to PDB ('$pdb'): chains[0] is '".$chains[0]."' (skipped)");
					next;
				}
				$start = $a[0];
				$stop = $a[1];

				if (($start < 1) or ($stop < 1))
				{
					die("Error: Start '$start', Stop '$stop' for '$name' PDB '$pdb' (chains $chains)");
					# next;
				}
				
				# No isoform information for dbReferences to PDB.
				# zcat input/uniprot_sprot_human.xml.gz | g -A3 'dbReference type="PDB"'|g -A3 "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-\d+"

				# Update unipdb
				# print "INSERT INTO unipdb SET name='$name', species='$species', pdb='$pdb', chains='$chains', start=$start, stop=$stop\n";
				Query("INSERT INTO unipdb SET name='$name', acc='$primary_acc', species='$species', pdb='$pdb', chains='$chains', start='$start', stop='$stop'") unless switch('debug');
			}
		}
		
		# PFAM
		# Example:
		# <dbReference id="PF00400" key="58" type="Pfam">
		#   <property type="entry name" value="WD40"/>
		#   <property type="match status" value="5"/>
		# </dbReference>
		if ($dbreference->{'att'}->{'type'} eq 'Pfam')
		{
			my $pfamid = $dbreference->{'att'}->{'id'};
			my $domain = '';
			my $occ = '';
			
		    die("Error: Properties not 2 for PFAM for '$name', PFAMID '$pfamid'") if ($dbreference->children_count('property') != 2);
		
			foreach $property ($dbreference->children('property'))
			{
				if ($property->{'att'}->{'type'} eq 'entry name')
				{
					$domain = $property->{'att'}->{'value'};
				}
				if ($property->{'att'}->{'type'} eq 'match status')
				{
					$occ = $property->{'att'}->{'value'};
				}
			}
		
			# No isoform information for dbReferences to Pfam.
			# zcat input/uniprot_sprot_human.xml.gz | g -A3 'dbReference.+Pfam'|g -A3 "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-\d+"

			# Update unipfam
			# print "INSERT INTO unipfam SET name='$name', species='$species', domain='$domain', pfamid='$pfamid', occ='$occ'\n";
			Query("INSERT INTO unipfam SET name='$name', acc='$primary_acc', species='$species', domain='$domain', pfamid='$pfamid', occ='$occ'") unless switch('debug');
		}

		# SUPFAM
		# Example:
		# <dbReference type="SUPFAM" id="SSF47370">
		#   <property type="entry name" value="Bromodomain"/>
		#   <property type="match status" value="1"/>
		# </dbReference>
		if ($dbreference->{'att'}->{'type'} eq 'SUPFAM')
		{
			my $supfamid = $dbreference->{'att'}->{'id'};
			my $domain = '';
			my $occ = '';
			
		    die("Error: Properties not 2 for SUPFAM for '$name', SUPFAMID '$supfamid'") if ($dbreference->children_count('property') != 2);
		
			foreach $property ($dbreference->children('property'))
			{
				if ($property->{'att'}->{'type'} eq 'entry name')
				{
					$domain = $property->{'att'}->{'value'};
				}
				if ($property->{'att'}->{'type'} eq 'match status')
				{
					$occ = $property->{'att'}->{'value'};
				}
			}
		
			# No isoform information for dbReferences to SUPFAM.
			# zcat input/uniprot_sprot_human.xml.gz | g -A3 'dbReference.+Pfam'|g -A3 "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-\d+"

			# Update unisupfam
			# print "INSERT INTO unisupfam SET name='$name', species='$species', domain='$domain', supfamid='$supfamid', occ='$occ'\n";
			Query("INSERT INTO unisupfam SET name='$name', acc='$primary_acc', species='$species', domain='".esc($domain)."', supfamid='$supfamid', occ='$occ'") unless switch('debug');
		}

		# InterPro
		# Example:
		# <dbReference id="IPR008816" key="8" type="InterPro"> 
		#   <property type="entry name" value="Rick_17kDa_Anti"/>
		# </dbReference> 
		if ($dbreference->{'att'}->{'type'} eq 'InterPro')
		{
			my $interproid = $dbreference->{'att'}->{'id'};
			my $domain = '';
			
		    die("Error: Properties not 1 for InterPro for '$name', InterPro ID '$interproid'") if ($dbreference->children_count('property') != 1);

			foreach $property ($dbreference->children('property'))
			{
				if ($property->{'att'}->{'type'} eq 'entry name')
				{
					$domain = $property->{'att'}->{'value'};
				}
			}

			# Update uniinterpro
			# print "INSERT INTO uniinterpro SET name='$name', species='$species', domain='$domain', interproid='$interproid'\n";
			Query("INSERT INTO uniinterpro SET name='$name', acc='$primary_acc', species='$species', domain='".esc($domain)."', interproid='$interproid'") unless switch('debug');
		}

		# GO terms
		# Example:
		# OLD:
		# <dbReference id="GO:0033644" key="9" type="Go"> 
		#   <property type="term" value="C:host cell membrane"/> 
		#   <property type="evidence" value="IEA:UniProtKB-SubCell"/>
		# </dbReference>
		# New:
        # <dbReference id="GO:0033644" type="GO">
        #   <property value="C:host cell membrane" type="term"/>
        #   <property value="ECO:0000501" type="evidence"/>
        #   <property value="UniProtKB-SubCell" type="project"/>
        # </dbReference>
		if ($dbreference->{'att'}->{'type'} eq 'GO')
		{
			my $goid = $dbreference->{'att'}->{'id'};
			my $gotype = '';
			my $goterm = '';
			my $evidence = '';
			my $evidencetype = '';
			my $project = '';
			
		    # die("Error: Properties not 3 for GO term for '$name', GO ID '$goid'") if ($dbreference->children_count('property') != 3);
		    die("Error: Properties not 3 or 2 for GO term for '$name', GO ID '$goid'") if (($dbreference->children_count('property') != 3) and ($dbreference->children_count('property') != 2));
			
			# Evidence code is sometimes not given (I guess in manual annotation cases by UniProt, e.g. GO:0042335 for APKC_DROME in UniProt 2017_12)
			$evidence = '';
			foreach $property ($dbreference->children('property'))
			{
				if ($property->{'att'}->{'type'} eq 'term')
				{
					$goterm = $property->{'att'}->{'value'};
				}
                if ($property->{'att'}->{'type'} eq 'evidence')
                {
                    $evidence = $property->{'att'}->{'value'};
                }
                if ($property->{'att'}->{'type'} eq 'project')
                {
                    $project = $property->{'att'}->{'value'};
                }
			}
			
			# Process GO term
			$goterm =~ /^(\w):(.+)/ or die("Error: Couldn't process GO term '$goterm' for name '$name'");
			($gotype, $goterm) = ($1, $2);
			if ($gotype eq 'C') { $gotype = 'CC'; }
			elsif ($gotype eq 'F') { $gotype = 'MF'; }
			elsif ($gotype eq 'P') { $gotype = 'BP'; }
			else { die("Error: Unknown GO type category '$gotype'"); }
			# Process evidence
			if ($evidence ne '')
			{
				$evidence =~ /^(\w+):(.+)/ or die("Error: Couldn't process GO evidence '$evidence' for name '$name'");
				($evidencetype, $evidence) = ($1, $2);
			}
			else
			{
				$evidence = '';
				$evidencetype = '';
			}
			
			# Update unigo
			$q = "INSERT INTO unigo SET name='$name', acc='$primary_acc', species='$species', goid='$goid', gotype='$gotype', goterm='".esc($goterm)."', evidencetype='$evidencetype', evidence='".esc($evidence)."', project='$project'";
			$q =~ s/=''/=NULL/g;
			# print "$q\n";
			Query($q) unless switch('debug');
		}

		# HGNC, MGI, RefSeq, IPI for uniid
		# Examples:
		#
		# <dbReference id="HGNC:4931" key="44" type="HGNC">
		#   <property type="gene designation" value="HLA-A"/>
		# </dbReference>
		#
		# <dbReference id="MGI:105305" key="22" type="MGI">
		#   <property type="gene designation" value="Slc1a5"/>
		# </dbReference>
		#
	    # <dbReference type="RefSeq" id="NP_001137430.1">
	    #   <molecule id="Q9NV96-2"/>
	    #   <property type="nucleotide sequence ID" value="NM_001143958.1"/>
	    # </dbReference>
		#
		# <dbReference id="IPI00144289" type="IPI"/>
		
		# HGNC/MGI
		if (contains($dbreference->{'att'}->{'type'}, ('HGNC', 'MGI')))
		{
			my $type = $dbreference->{'att'}->{'type'};
			# my $internal = $dbreference->{'att'}->{'id'};
			my $value = '';

		    die("Error: Properties not 1 for '$type' for '$name', internal ID '$internal'") if ($dbreference->children_count('property') != 1);
			# die("Error: Type not 'gene designation' for '$type' for '$name', internal ID '$internal'") if ($dbreference->first_child('property')->{'att'}->{'type'} ne 'gene designation');
			die("Error: Type not 'gene designation' for '$type' for '$name'") if ($dbreference->first_child('property')->{'att'}->{'type'} ne 'gene designation');

			$value = $dbreference->first_child('property')->{'att'}->{'value'};
			
			# Get isoform if given
			$isoacc = $primary_acc;
			if ($dbreference -> children_count('molecule') > 1)
			{
				die("Error: Multiple molecules for HGNC/MGI dbReference for '$name'");
			}
			if ($dbreference -> children_count('molecule') == 1)
			{
				$isoacc = $dbreference->first_child('molecule')->{'att'}->{'id'};
			}
			
			# Update uniid
			# print "INSERT INTO uniid SET name='$name', species='$species', type='$type', internal='$internal', value='$value'\n";
			Query("INSERT INTO uniid SET name='$name', acc='$isoacc', canon='$primary_acc', species='$species', type='$type', value='$value'") unless switch('debug');
		}
		
		# RefSeq
		if ($dbreference->{'att'}->{'type'} eq 'RefSeq')
		{
			my $type = $dbreference->{'att'}->{'type'};
			# my $internal = $dbreference->{'att'}->{'id'};
			my $value = '';

		    die("Error: Properties not 1 for '$type' for '$name', internal ID '$internal'") if ($dbreference->children_count('property') != 1);
			# die("Error: Type not 'nucleotide sequence ID' for '$type' for '$name', internal ID '$internal'") if ($dbreference->first_child('property')->{'att'}->{'type'} ne 'nucleotide sequence ID');
			die("Error: Type not 'nucleotide sequence ID' for '$type' for '$name'") if ($dbreference->first_child('property')->{'att'}->{'type'} ne 'nucleotide sequence ID');

			$value = $dbreference->first_child('property')->{'att'}->{'value'};

			# Get isoform if given
			$isoacc = $primary_acc;
			if ($dbreference -> children_count('molecule') > 1)
			{
				die("Error: Multiple molecules for RefSeq dbReference for '$name'");
			}
			if ($dbreference -> children_count('molecule') == 1)
			{
				$isoacc = $dbreference->first_child('molecule')->{'att'}->{'id'};
			}
			
			# Add unversioned value
			$valueshort = $value;
			$valueshort =~ s/\.\d+$//;

			# Update uniid
			# print "INSERT INTO uniid SET name='$name', species='$species', type='$type', internal='$internal', value='$value'\n";
			Query("INSERT INTO uniid SET name='$name', acc='$isoacc', canon='$primary_acc', species='$species', type='$type', value='$value', valueshort='$valueshort'") unless switch('debug');
		}

		# if ($dbreference->{'att'}->{'type'} eq 'IPI')
		# {
		# 	my $type = $dbreference->{'att'}->{'type'};
		# 	my $internal = $dbreference->{'att'}->{'id'};
		# 	my $value = '';
		#
		#     die("Error: Properties not 0 for '$type' for '$name', internal ID '$internal'") if ($dbreference->children_count('property') != 0);
		#
		# 	# Update uniid
		# 	# print "INSERT INTO uniid SET name='$name', species='$species', type='$type', internal='$internal', value='$value'\n";
		# 	Query("INSERT INTO uniid SET name='$name', species='$species', type='$type', internal='$internal', value='$value'") unless switch('debug');
		# }
		
		# Proteomes
		# Example:
	  	# <dbReference id="UP000005640" type="Proteomes"/>
		if ($dbreference->{'att'}->{'type'} eq 'Proteomes')
		{
			my $proteomeid = $dbreference->{'att'}->{'id'};

			# Update uniproteome
			# print "INSERT INTO uniproteome SET name='$name', species='$species', proteome='$proteomeid'\n";
			Query("INSERT INTO uniproteome SET name='$name', acc='$primary_acc', species='$species', proteome='$proteomeid'") unless switch('debug');
		}

		# NCBI Taxonomy (NCBI Taxon IDs) # NOW A SEPARATE ENTRY: ORGANISM.
		# Example:
		# 
		# Old:
		# <dbReference id="9606" key="1" type="NCBI Taxonomy"/>
		# 
		# New:
		  # <organism>
		  #   <name type="scientific">Pasteurella multocida (strain Pm70)</name>
		  #   <dbReference id="272843" type="NCBI Taxonomy"/>
		  #   <lineage>
		  #     <taxon>Bacteria</taxon>
		  #     <taxon>Proteobacteria</taxon>
		  #     <taxon>Gammaproteobacteria</taxon>
		  #     <taxon>Pasteurellales</taxon>
		  #     <taxon>Pasteurellaceae</taxon>
		  #     <taxon>Pasteurella</taxon>
		  #   </lineage>
		  # </organism>
		# 	    
		# 
		# if ($dbreference->{'att'}->{'type'} eq 'NCBI Taxonomy')
		# {
		# 	my $type = $dbreference->{'att'}->{'type'};
		# 	my $internal = $dbreference->{'att'}->{'id'};
		# 	my $value = '';
		# 	
		# 	# Update unitax, if this taxon ID isn't in there yet
		# 	$query = Query("SELECT species FROM unitax WHERE tax='$internal'");
		# 	if (Numrows($query) == 0)
		# 	{
		# 		Query("INSERT INTO unitax SET tax='$internal', species='$species'") unless switch('debug');
		# 	}
		# 	else
		# 	{
		# 		($thisspec) = FetchOne($query);
		# 		die("Error: Taxon ID '$internal' already exists in unitax, but it's linked to '$thisspec' instead of '$species'") if ($thisspec ne $species);
		# 	}
		# }
	}
	
	# Comments
	# 
	# Includes:
	# Isoform comments
	# Subcellular location
	# Function
	# 
	# Example:
    # <comment type="subcellular location">
    #   <molecule>Isoform 2</molecule>
    #   <subcellularLocation>
    #     <location evidence="7">Golgi apparatus</location>
    #   </subcellularLocation>
    # </comment>
	# 
	# <evidence key="1" type="ECO:0000006">
	#   <source>
	#     <dbReference id="16766513" type="PubMed"/>
	#   </source>
	# </evidence>
	#
	# Isoform comment (not yet parsed):
    # <comment type="miscellaneous">
    #   <molecule>Isoform 2</molecule>
    #   <text evidence="6">Based on a readthrough transcript which may produce a PTGES3L-AARSD1 fusion protein.</text>
    # </comment>
	
	%isoform_accs = ();
	# %isoform_names = ();

    foreach $comment ($x -> children('comment'))
	{
		# print "\ntype=['".$comment->{'att'}->{'type'}."']\n\n";
		addme("Comment types", $comment->{'att'}->{'type'});

		# "Alternative products" (Isoforms)
		# 
		# Example:
		# 
	    # <comment type="alternative products">
	    #   <event type="alternative splicing"/>
	    #   <isoform>
	    #     <id>P63034-1</id>		# Listed first at https://www.uniprot.org/uniprotkb/P63034/entry#sequences
	    #     <id>P97695-1</id>		# Listed second (P63034-1, P97695-1). P97695 is actually an obsolete entry ("demerged" back in 2004! https://www.uniprot.org/uniprotkb/P97695/history), and P97695-1 returns no hits on UniProt. I'll simply ignore such alternative accessions.
	    #     <name>1</name>		# Listed as the "Name" of the isoform
	    #     <name>CLM2-A</name>	# Listed under "Synonyms" (I have a separate column for synonyms separated by |)
	    #     <sequence type="displayed"/>
	    #   </isoform>
	    #   <isoform>
	    #     <id>P63034-2</id>
	    #     <id>P97695-2</id>
	    #     <name>2</name>
	    #     <sequence type="described" ref="VSP_006038"/>
	    #   </isoform>
	    #   <isoform>
	    #     <id>P63034-3</id>
	    #     <id>P97695-3</id>
	    #     <name>3</name>
	    #     <name>CLM2-B</name>
	    #     <sequence type="described" ref="VSP_006037"/>
	    #   </isoform>
	    # </comment>
		
		if ($comment->{'att'}->{'type'} eq 'alternative products')
		{
			# die("Error: Event type is '".$comment -> first_child('event') -> {'att'} -> {'type'}."' instead of 'alternative splicing' for alternative products in '$name'") if ($comment -> first_child('event') -> {'att'} -> {'type'} ne 'alternative splicing');
			my $isoform_event = '';
		    if ($comment -> children_count('event') == 0)
			{
				die("Error: No isoform events in '$name'");
			}
		    elsif ($comment -> children_count('event') == 1)
			{
				$isoform_event = $comment -> first_child('event') -> {'att'} -> {'type'};
			}
		    elsif ($comment -> children_count('event') > 1)
			{
				# die("Error: Multiple isoform events in '$name'");
				addme("multiple isoform events for name (kept)", $name);
				
				# Combine into one string
				$isoform_event = '';
				
				foreach $event ($comment -> children('event'))
				{
					$isoform_event .= $event -> {'att'} -> {'type'}."|";
				}
				$isoform_event =~ s/\|$//;
			}
			addme("total isoform events", $isoform_event);
			addme("isoform event is '$isoform_event' for name", $name);
			
			$first_isoform = 1;
			foreach $isoform ($comment -> children('isoform'))
			{
			    # die("Error: Multiple accessions ('id') for isoforms in '$name'") if ($isoform -> children_count('id') != 1);
			    if ($isoform -> children_count('id') != 1)
				{
					addme("multiple isoform accessions for same isoform entry for name (kept)", $name);
				}
			
				# foreach $isoform_acc ($isoform -> children('id'))
				# {
				# 	$isoform_acc = $isoform_acc -> text;
					
					# Only use the first isoform accession (see P63034 example above)
					$isoform_acc = $isoform -> first_child('id') -> text;
					
					addme("isoform event is '$isoform_event' for acc", $isoform_acc);
				
					my $isoform_canon = $isoform_acc;
					$isoform_canon =~ s/-\d+$//;
				
					if ($isoform_canon ne $primary_acc)
					{
						# die("Error: Isoform accession is '$isoform_acc' for canonical accession '$primary_acc' in '$name'");
						addme("isoform accession base is different from canonical accession for name (kept)", $name);
					}
				
					# # Remove -1 from acc for canonical sequence
					# See below
					# $isoform_acc =~ s/-1$//;
				
					# die("Error: Multiple names ('name') for isoforms in '$name'") if ($isoform  ->  children_count('name') != 1);
					# my $isoform_name = $isoform -> {'att'} -> {'name'};

				    die("Error: Multiple sequences for isoforms in '$name'") if ($isoform -> children_count('sequence') != 1);
					die("Error: No sequence type for isoform in '$name'") if (!$isoform -> first_child('sequence') -> att_exists('type'));
					my $isoform_sequencetype = $isoform -> first_child('sequence') -> {'att'} -> {'type'};

					addme("total isoform sequence type", $isoform_sequencetype);
					addme("isoform sequence type is '$isoform_sequencetype' for name", $name);
					addme("isoform sequence type is '$isoform_sequencetype' for acc", $isoform_acc);

					if (($first_isoform == 1) and ($isoform_sequencetype ne 'displayed'))
					{
						die("Error: First isoform sequence type is not 'displayed' in '$name'");
					}
					elsif (($first_isoform == 0) and ($isoform_sequencetype eq 'displayed'))
					{
						die("Error: Later isoform sequence type is 'displayed' in '$name'");
					}
				
					# Isoform names (previously 'text')
					# Storing these in an "isoform_accs" hash so that "molecule names" ("Isoform 4") can be translated to isoform accessions ("Q.....-3")
					# die("Error: Multiple texts for isoform in '$name'") if ($isoform -> children_count('text') > 1);
					my $isoform_text = '';
					my $isoform_text_aliases = '';
					# if ($isoform -> children_count('text') == 1)
					# {
					# 	$isoform_text = $isoform -> first_child('text') -> text;
					# }
					# This is normal:
					# if ($isoform -> children_count('name') > 1)
					# {
					# 	die("Error: Multiple names for isoform for '$name'")
					# }
					if ($isoform -> children_count('name') > 0)
					{
						# $isoform_text = join("|", $isoform -> children_text('name'));
						# $isoform_name = $isoform -> first_child('name')->text;
						foreach $isoform_name ($isoform -> children_text('name'))
						{
							# First isoform name
							if ($isoform_text eq '')
							{
								$isoform_text = $isoform_name;
							}
							else
							{
								# Aliases
								if ($isoform_text_aliases eq '')
								{
									$isoform_text_aliases = $isoform_name;
								}
								else
								{
									$isoform_text_aliases .= "|$isoform_name";
									# More readable: //
									# $isoform_text_aliases .= " // $isoform_name";
								}
							}
						
							# Not using this hash anyway
							# # There can be multiple names for one isoform acc:
							# 					    # <isoform>
							# 					    #   <id>A0MD28-1</id>
							# 					    #   <name>Replicase polyprotein 1ab</name>
							# 					    #   <name>pp1ab</name>
							# 					    #   <sequence type="external"/>
							# 					    # </isoform>
							# if (exists($isoform_names{$isoform_acc}))
							# {
							# 	$isoform_names{$isoform_acc} .= "|$isoform_name";
							# }
							# else
							# {
							# 	$isoform_names{$isoform_acc} = $isoform_name;
							# }
							if (exists($isoform_accs{$isoform_name}))
							{
								# die("Error: Isoform name '$isoform_name' is mapped to multiple isoform accs ('$isoform_accs{$isoform_name}', '$isoform_acc')")
								$isoform_accs{$isoform_name} .= "|$isoform_acc";
							}
							else
							{
								$isoform_accs{$isoform_name} = $isoform_acc;
							}
						}
					}
				
					$isoform_alternative = 1;
					if ($isoform_sequencetype eq 'displayed')
					{
						$isoform_alternative = 0;
					}
				
					if (($isoform_canon ne $primary_acc) and ($isoform_alternative == 0))
					{
						die("Error: Canonical isoform accession is '$isoform_acc' for canonical accession '$primary_acc' in '$name'");
					}
				
					# This is very rare in 2021_01 (just 7 or so):
					$isoform_confirmed = 1;
					if ($isoform_text =~ /No experimental confirmation available\./)
					{
						$isoform_confirmed = 0;
					}
				
					# # Skip "-1" (canonical) isoforms. If I don't, they'd still only be added for proteins with isoforms, so it's not a complete set of sequences then either.
					# Actually, proteins like https://www.uniprot.org/uniprot/Q6ZT62#sequences can have multiple isoforms with a -1 suffix (some are in external entries under their own accs). Only the "first_isoform", which is listed first in the XML, is then canonical. Can't use the -1 suffix to determine which is first.
					if (($first_isoform == 1) xor ($isoform_acc =~ /-1$/))
					{
						$isoform_acc =~ /(-\d+)$/ or die("Error: Couldn't parse isoform suffix from '$isoform_acc'");
						my $isoform_suffix = $1;
						if ($isoform_suffix ne '-1')
						{
							# die("first_isoform is '$first_isoform' for isoform_acc '$isoform_acc'") if (($first_isoform == 1) xor ($isoform_acc =~ /-1$/));
							addme("first_isoform is '$first_isoform' but isoform_acc doesn't end in '-1' for acc (kept)", $primary_acc);
							addme("first_isoform is '$first_isoform' but isoform_acc doesn't end in '-1' for isoform acc (kept)", $isoform_acc);
						}
						else
						{
							# die("first_isoform is '$first_isoform' for isoform_acc '$isoform_acc'") if (($first_isoform == 1) xor ($isoform_acc =~ /-1$/));
							addme("first_isoform is '$first_isoform' but isoform_acc ends in '-1' for acc (kept)", $primary_acc);
							addme("first_isoform is '$first_isoform' but isoform_acc ends in '-1' for isoform acc (kept)", $isoform_acc);
						}
					}
				
					
					# Handle first isoform (previously, these were skipped for table 'uniiso', but then I can't record their names etc.)
					my $isoform_seq = '';
					my $isoform_trembl = '';
					if ($first_isoform == 1)
					{
						# First isoform handled
						$first_isoform = 0;

						# This canonical isoform acc shouldn't be in the varsplic file
						if (exists($varsplic{$isoform_acc}))
						{
							die("Error: '$isoform_acc' is the first isoform, but it has a sequence in the varsplic file");
						}
						# Use canonical sequence for it
						$isoform_seq = $seq;
						$isoform_trembl = $trembl;
					}
					else
					{
						# This canonical isoform acc should be in the varsplic file
						if (!exists($varsplic{$isoform_acc}))
						{
							if ($isoform_sequencetype eq 'not described')
							{
								# Undefined sequence (gets skipped below)
								$isoform_seq = '';
								$isoform_trembl = $trembl;
							}
							else
							{
								$isoform_acc =~ /(-\d+)$/ or die("Error: Couldn't parse isoform suffix from '$isoform_acc'");
								my $isoform_suffix = $1;

								# die("Error: '$isoform_acc' isn't the first isoform, but it does not have a sequence in the varsplic file")
								addme("isoform isn't the first isoform, but it does not have a sequence in the varsplic file for acc (skipped)", $primary_acc);
								addme("isoform isn't the first isoform, but it does not have a sequence in the varsplic file for isoform acc (skipped)", $isoform_acc);
								addme("isoform isn't the first isoform, but it does not have a sequence in the varsplic file for isoform suffix (skipped)", $isoform_suffix);

								# Undefined sequence (gets skipped below)
								$isoform_seq = '';
								$isoform_trembl = $trembl;
							}
						}
						else
						{
							# Use varsplic alternative isoform sequence for it
							$isoform_seq = $varsplic{$isoform_acc};
							$isoform_trembl = $varsplic_trembl{$isoform_acc};
						}
					}
				
					# Both should always agree
					if ($isoform_trembl ne $trembl)
					{
						die("Error: Isoform '$isoform_acc' has trembl '$isoform_trembl', but primary accession '$primary_acc' has trembl '$trembl' (shouldn't happen)");
					}
					
						# Examples where -2, -3 etc. can be the default isoform:
						# https://www.uniprot.org/uniprot/Q8N423#Q8N423-2
						# https://www.uniprot.org/uniprot/Q16635#Q16635-3
						# https://www.uniprot.org/uniprot/Q6ZMS7#Q6ZMS7-2
						# https://www.uniprot.org/uniprot/Q9Y2M2#Q9Y2M2-2
						# Should therefore keep information on all listed isoforms.
						
						# No point inserting isoforms without sequences
						if ($isoform_seq ne '')
						{
							# Update uniiso (Isoform accession -> isoform annotation and existence evidence table)
							# For now, without sequence (which will be added by uniiso_seq.pl)
							$q = "INSERT INTO uniiso SET name='$name', acc='$isoform_acc', canon='$primary_acc', species='$species', event='$isoform_event', alternative='$isoform_alternative', sequencetype='$isoform_sequencetype', confirmed='$isoform_confirmed', isoform='".esc($isoform_text)."', isoform_aliases='".esc($isoform_text_aliases)."', seq='$isoform_seq'";
							$q =~ s/=''/=NULL/g;
							Query($q) unless switch('debug');
					
							# Also insert isoform into uniseq (with a new type, 'UniIso', so it doesn't interfere with any existing scripts).
							# # For now, without sequence (which will be added by uniiso_seq.pl)
							# Don't insert primary isoform into uniseq (alternative=0)
							if ($isoform_alternative != 0)
							{
								$q = "INSERT INTO uniseq SET name='$name', acc='$isoform_acc', canon='$primary_acc', species='$species', type='UniIso', trembl='$isoform_trembl', seq='$isoform_seq'";
								$q =~ s/=''/=NULL/g;
								Query($q) unless switch('debug');
							}
						}
						
					
					# }
					# else
					# {
					# 	# Caveat about isoform names ('isoform'): The numeric ones sometimes don't match the acc's -2, -3, -4 suffix, so I do always need to use the "molecule's" ID, and not convert the numeric isoform name into a suffix.
					# 	# Examples: SELECT * FROM uniiso WHERE isoform REGEXP '^[0-9]+$' AND acc NOT LIKE CONCAT('%-', isoform);
					#
					#
					# 	# For later isoforms:
					# 	# Update uniiso (Isoform accession -> isoform annotation and existence evidence table)
					# 	# For now, without sequence (which will be added by uniiso_seq.pl)
					# 	$q = "INSERT INTO uniiso SET name='$name', acc='$isoform_acc', canon='$primary_acc', species='$species', event='$isoform_event', alternative='$isoform_alternative', sequencetype='$isoform_sequencetype', confirmed='$isoform_confirmed', isoform='".esc($isoform_text)."', isoform_aliases='".esc($isoform_text_aliases)."'";
					# 	$q =~ s/=''/=NULL/g;
					# 	Query($q) unless switch('debug');
					#
					# 	# Also insert isoform into uniseq (with a new type, 'UniIso', so it doesn't interfere with any existing scripts).
					# 	# For now, without sequence (which will be added by uniiso_seq.pl)
					# 	$q = "INSERT INTO uniseq SET name='$name', acc='$isoform_acc', canon='$primary_acc', species='$species', type='UniIso', trembl='$trembl'";
					# 	$q =~ s/=''/=NULL/g;
					# 	Query($q) unless switch('debug');
					# }
				# }
				
				# # Once all accessions (<ids>) for this <isoform> entry are handled, consider the first isoform handled
				# Actually: Only use the first accession (<id>) for each <isoform> entry.
				# $first_isoform = 0;
			}
		}
	}

    foreach $comment ($x -> children('comment'))
	{
		# Subcellular location
		# 
		# Example:
		# 
	    # <comment type="subcellular location">
	    #   <molecule>Isoform 2</molecule>
	    #   <subcellularLocation>
	    #     <location evidence="7">Golgi apparatus</location>
	    #   </subcellularLocation>
	    # </comment>

		if ($comment->{'att'}->{'type'} eq 'subcellular location')
		{
			# Check if there is 'evidence' information
			my $evid = '';
			my $sources = '';
			foreach $subcellularlocation ($comment -> children('subcellularLocation'))
			{
				foreach $location ($subcellularlocation->children('location'))
				{
					if ($location->att_exists('evidence'))
					{
						my $evidence = $location->{'att'}->{'evidence'};
						my @keys = split(/ /, $evidence);
				
						my @sources = ();
						my @evid = ();
						foreach $evidence ($x -> children('evidence'))
						{
							if ($evidence -> children_count('source') == 0)
							{
								# die("Error: No PubMed source for subcellular location comment evidence for '$name'")
								next;
							}
							if ($evidence -> children_count('source') > 1)
							{
								die("Error: Multiple sources for evidence for '$name'");
							}
					
							if ((contains($evidence->{'att'}->{'key'}, @keys)))
							{
								push(@evid, $evidence->{'att'}->{'type'});
						
								foreach $dbref ($evidence -> first_child('source') -> children('dbReference'))
								{
									if ($dbref->{'att'}->{'type'} eq 'PubMed')
									{
										push(@sources, $dbref->{'att'}->{'id'});
									}
									else
									{
										# die("Error: Evidence is not PubMed for '$name' PTM (evidence key '".$evidence->{'att'}->{'key'}."')");
									}
								}
							}
						}
				
						# if (scalar(unique(@evid)) == 1) { @evid = unique(@evid); }
						# $evid = join("|", @evid);
						# $sources = join("|", @sources);
						$evid = join("|", unique(@evid));
						$sources = join("|", unique(@sources));
					}
				}
			}

		    # die("Error: Subcellular location duplicates for '$name'") if ($comment->children_count('subcellularLocation') != 1);
			foreach $subcellularlocation ($comment -> children('subcellularLocation'))
			{
				foreach $location ($subcellularlocation->children('location'))
				{
					my $loc = $location->text;
						
					# Get isoform information
				    # <comment type="subcellular location">
				    #   <molecule>Isoform 2</molecule>
				    #   <subcellularLocation>
				    #     <location evidence="7">Golgi apparatus</location>
				    #   </subcellularLocation>
				    # </comment>
					my $isoaccs = $primary_acc;
					if ($subcellularlocation -> children_count('molecule') > 1)
					{
						die("Error: Multiple molecules for subcellularLocation comment for '$name'");
					}
					elsif ($subcellularlocation -> children_count('molecule') == 1)
					{
						my $moleculename = $subcellularlocation->first_child('molecule')->text;
						
						if ($moleculename =~ /^Isoform /)
						{
						    # <isoform>
						    #   <id>Q2V2M9-1</id>
						    #   <name>1</name>
						    #   <sequence type="displayed"/>
						    # </isoform>
							
							# Strip "Isoform " (turning e.g. "Isoform 4" into "4")
							$moleculename =~ s/^Isoform //;
						}
						
						if (!exists($isoform_accs{$moleculename}))
						{
							addme("no isoform acc in isoform_accs hash for name|moleculename (used primary_acc instead)", "$name|$moleculename");
							addme("no isoform acc in isoform_accs hash for '$species' protein for name (used primary_acc instead)", $name);
						}
						else
						{
							$isoaccs = $isoform_accs{$moleculename};
						}
					}
					
					# There can be multiple isoform accs for a given moleculename
					foreach $isoacc (split(/\|/, $isoaccs))
					{
						# Update uniloc (Name -> Subcellular location table)
						$q = "INSERT INTO uniloc SET name='$name', acc='$isoacc', canon='$primary_acc', species='$species', loc='$loc', evidence='$evid', sources='$sources'";
						$q =~ s/=''/=NULL/g;
						Query($q) unless switch('debug');
					}
				}
			}
		}
	}

    foreach $comment ($x -> children('comment'))
	{
		# Function
		# 
		# Example:
		# 
		# SELECT * FROM unifunc WHERE name='FHOD3_HUMAN';
	    # <comment type="function">
	    #   <text evidence="1 9">Actin-organizing protein that may cause stress fiber formation together with cell elongation (By similarity). Isoform 4 may play a role in actin filament polymerization in cardiomyocytes.</text>
	    # </comment>
	    # <comment type="function">
	    #   <molecule>Isoform 2</molecule>
	    #   <text evidence="8">Metallocarboxypeptidase that mediates tubulin deglutamylation.</text>
	    # </comment>
		#
		# Example with two "isoforms" (XRCC4_HUMAN), though the latter isn't actually defined as an isoform ("C-terminus"). Its function instead gets added under the primary accession.
		#
		# SELECT * FROM unifunc WHERE name='XRCC4_HUMAN';
	    # <comment type="function">
	    #   <molecule>DNA repair protein XRCC4</molecule>
	    #   <text evidence="3 4 8 13 14 17 18 24 26 27 28 30 31 32 36 39 42 44 47 48 49 50 57 58">DNA non-homologous end joining (NHEJ) core factor, required for double-strand break repair and V(D)J recombination (PubMed:10757784, PubMed:10854421, PubMed:17124166, PubMed:16412978, PubMed:8548796, PubMed:25742519, PubMed:12517771, PubMed:17290226, PubMed:22228831, PubMed:25597996, PubMed:25934149, PubMed:26100018, PubMed:26774286). Acts as a scaffold protein that regulates recruitment of other proteins to DNA double-strand breaks (DSBs) (PubMed:15385968, PubMed:20852255, PubMed:26774286, PubMed:27437582). Associates with NHEJ1/XLF to form alternating helical filaments that bridge DNA and act like a bandage, holding together the broken DNA until it is repaired (PubMed:26100018, PubMed:27437582, PubMed:28500754, PubMed:21775435, PubMed:22287571, PubMed:21768349). The XRCC4-NHEJ1/XLF subcomplex binds to the DNA fragments of a DSB in a highly diffusive manner and robustly bridges two independent DNA molecules, holding the broken DNA fragments in close proximity to one other (PubMed:27437582). The mobility of the bridges ensures that the ends remain accessible for further processing by other repair factors (PubMed:27437582). Plays a key role in the NHEJ ligation step of the broken DNA during DSB repair via direct interaction with DNA ligase IV (LIG4): the LIG4-XRCC4 subcomplex reseals the DNA breaks after the gap filling is completed (PubMed:9242410, PubMed:10757784, PubMed:10854421, PubMed:12517771, PubMed:17290226, PubMed:19837014). XRCC4 stabilizes LIG4, regulates its subcellular localization and enhances LIG4's joining activity (PubMed:9242410, PubMed:10757784, PubMed:10854421, PubMed:12517771, PubMed:17290226, PubMed:21982441, PubMed:22228831). Binding of the LIG4-XRCC4 subcomplex to DNA ends is dependent on the assembly of the DNA-dependent protein kinase complex DNA-PK to these DNA ends (PubMed:10757784, PubMed:10854421). Promotes displacement of PNKP from processed strand break termini (PubMed:20852255, PubMed:28453785).</text>
	    # </comment>
	    # <comment type="function">
	    #   <molecule>Protein XRCC4, C-terminus</molecule>
	    #   <text evidence="55">Acts as an activator of the phospholipid scramblase activity of XKR4 (PubMed:33725486). This form, which is generated upon caspase-3 (CASP3) cleavage, translocates into the cytoplasm and interacts with XKR4, thereby promoting phosphatidylserine scramblase activity of XKR4 and leading to phosphatidylserine exposure on apoptotic cell surface (PubMed:33725486).</text>
	    # </comment>
		#
		# Example where the isoform annotation works nicely:
		# 
		# SELECT * FROM unifunc WHERE name='ACOX1_HUMAN';
	    # <comment type="function">
	    #   <text evidence="6 8 9 13 14 15">Involved in the initial and rate-limiting step of peroxisomal beta-oxidation of straight-chain saturated and unsaturated very-long-chain fatty acids (PubMed:7876265, PubMed:15060085, PubMed:17458872, PubMed:17603022, PubMed:32169171, PubMed:33234382). Catalyzes the desaturation of fatty acyl-CoAs such as palmitoyl-CoA (hexadecanoyl-CoA) to 2-trans-enoyl-CoAs ((2E)-enoyl-CoAs) such as (2E)-hexadecenoyl-CoA, and donates electrons directly to molecular oxygen (O(2)), thereby producing hydrogen peroxide (H(2)O(2)) (PubMed:7876265, PubMed:17458872, PubMed:17603022).</text>
	    # </comment>
	    # <comment type="function">
	    #   <molecule>Isoform 1</molecule>
	    #   <text evidence="9">Shows highest activity against medium-chain fatty acyl-CoAs. Shows optimum activity with a chain length of 10 carbons (decanoyl-CoA) in vitro.</text>
	    # </comment>
	    # <comment type="function">
	    #   <molecule>Isoform 2</molecule>
	    #   <text evidence="9">Is active against a much broader range of substrates and shows activity towards long-chain fatty acyl-CoAs.</text>
	    # </comment>
	    # <evidence type="ECO:0000269" key="6">
	    #   <source>
	    #     <dbReference type="PubMed" id="15060085"/>
	    #   </source>
	    # </evidence>
	    # <evidence type="ECO:0000269" key="8">
	    #   <source>
	    #     <dbReference type="PubMed" id="17458872"/>
	    #   </source>
	    # </evidence>
	    # <evidence type="ECO:0000269" key="9">
	    #   <source>
	    #     <dbReference type="PubMed" id="17603022"/>
	    #   </source>
	    # </evidence>
	    # <evidence type="ECO:0000269" key="13">
	    #   <source>
	    #     <dbReference type="PubMed" id="32169171"/>
	    #   </source>
	    # </evidence>
	    # <evidence type="ECO:0000269" key="14">
	    #   <source>
	    #     <dbReference type="PubMed" id="33234382"/>
	    #   </source>
	    # </evidence>
	    # <evidence type="ECO:0000269" key="15">
	    #   <source>
	    #     <dbReference type="PubMed" id="7876265"/>
	    #   </source>
	    # </evidence>
		
		if ($comment->{'att'}->{'type'} eq 'function')
		{
			if ($comment->children_count('text') != 1)
			{
				die("Error: More than one text tag under function comment for name '$name'");
			}
			
			# Get function comment text
			my $func = $comment->first_child('text')->text;
			
			# Get isoform if given
			my $isoaccs = $primary_acc;
			if ($comment -> children_count('molecule') > 1)
			{
				die("Error: Multiple molecules for function comment for '$name'");
			}
			if ($comment -> children_count('molecule') == 1)
			{
				my $moleculename = $comment->first_child('molecule')->text;
				
				if ($moleculename =~ /^Isoform /)
				{
				    # <isoform>
				    #   <id>Q2V2M9-1</id>
				    #   <name>1</name>
				    #   <sequence type="displayed"/>
				    # </isoform>
					
					# Strip "Isoform " (turning e.g. "Isoform 4" into "4")
					$moleculename =~ s/^Isoform //;
				}
				
				if (!exists($isoform_accs{$moleculename}))
				{
					addme("no isoform acc in isoform_accs hash for name|moleculename (used primary_acc instead)", "$name|$moleculename");
				}
				else
				{
					$isoaccs = $isoform_accs{$moleculename};
				}
			}

			# Check if there is 'evidence' information
			my $evid = '';
			my $sources = '';
			
			# if (!$comment->first_child('text')->att_exists('evidence'))
			# {
			# 	die("Error: No evidence attribute for '$name' in 'function' Comment with text '$func'");
			# }
			if ($comment->first_child('text')->att_exists('evidence'))
			{
				my $evidence = $comment->first_child('text')->{'att'}->{'evidence'};
				my @keys = split(/ /, $evidence);
				
				my @sources = ();
				my @evid = ();
				foreach $evidence ($x -> children('evidence'))
				{
					# Sometimes there is no PubMed ID, e.g. in:
				    # <comment type="function">
				    #   <text evidence="1">Plays a role in virus cell tropism, and may be required for efficient virus replication in macrophages.</text>
				    # </comment>
				    # <evidence type="ECO:0000250" key="1"/>	# sequence similarity evidence used in manual assertion
					
					# if ($evidence -> children_count('source') == 0)
					# {
					# 	die("Error: No PubMed source for Comment evidence for '$name'")
					# 	# next;
					# }
					if ($evidence -> children_count('source') > 1)
					{
						die("Error: Multiple sources for evidence for '$name'");
					}
					
					if ((contains($evidence->{'att'}->{'key'}, @keys)))
					{
						push(@evid, $evidence->{'att'}->{'type'});
						
						if ($evidence -> children_count('source') > 0)
						{
							foreach $dbref ($evidence -> first_child('source') -> children('dbReference'))
							{
								if ($dbref->{'att'}->{'type'} eq 'PubMed')
								{
									push(@sources, $dbref->{'att'}->{'id'});
								}
								else
								{
									# die("Error: Evidence is not PubMed for '$name' PTM (evidence key '".$evidence->{'att'}->{'key'}."')");
								}
							}
						}
					}
				}
				
				# if (scalar(unique(@evid)) == 1) { @evid = unique(@evid); }
				# $evid = join("|", @evid);
				# $sources = join("|", @sources);
				$evid = join("|", unique(@evid));
				$sources = join("|", unique(@sources));
			}

			# There can be multiple isoform accs for a given moleculename
			foreach $isoacc (split(/\|/, $isoaccs))
			{
				# Update unifunc (Name -> Protein Function table)
				$q = "INSERT INTO unifunc SET name='$name', acc='$isoacc', canon='$primary_acc', species='$species', func='".esc($func)."', evidence='$evid', sources='$sources'";
				$q =~ s/=''/=NULL/g;
				Query($q) unless switch('debug');
			}
		}
	}

	# Features
	# 
	# Includes:
	# PTM - modified residue
	# PTM - glycosylation site
	# PTM - lipid moiety-binding region (called a region, but it's always one residue)
	# PTM - crosslink (ubiquitin)
	# Variation (SNPs)
	# Misc features (unifeat)
	#  
    foreach $feature ($x -> children('feature'))
	{
		# print "\ntype=['".$feature->{'att'}->{'type'}."']\n\n";
		addme("Feature types", $feature->{'att'}->{'type'});
		# print "   >> ".$feature->{'att'}->{'type'}."\n" if (switch('debug'));

		# Modified residue (most PTMs) + glycosylation sites + lipid moiety-binding region + cross-link (ubiquitin) (all in one here)
		# 
		# Example:
		# <feature type="modified residue" status="by similarity" description="N-acetylmethionine; in 14-3-3 protein beta/alpha; alternate">
		#   <location>
		#     <position position="1"/>
		#   </location>
		# </feature>
		# 
		# <feature description="N-linked (GlcNAc...)" status="by similarity" type="glycosylation site">
		# <location>
		# <position position="110"/>
		# </location>
		# </feature>
		# 
		# <feature type="cross-link" status="by similarity" description="Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)">
		#   <location>
		#     <position position="334"/>
		#   </location>
		# </feature>
		#
		# New, isoform-specific:
	    # <feature type="modified residue" description="Phosphothreonine" evidence="9">
	    #   <location sequence="Q2V2M9-4">
	    #     <position position="1474"/>
	    #   </location>
	    # </feature>
	    # <feature type="modified residue" description="Phosphothreonine" evidence="9">
	    #   <location sequence="Q2V2M9-4">
	    #     <position position="1476"/>
	    #   </location>
	    # </feature>
		
		if (contains($feature->{'att'}->{'type'}, ('modified residue', 'glycosylation site', 'lipid moiety-binding region', 'cross-link')))
		{
			my $type = $feature->{'att'}->{'type'};
			
			

			# Check if there is 'evidence' information
			my $evid = '';
			my $sources = '';
			if ($feature->att_exists('evidence'))
			{
				my $evidence = $feature->{'att'}->{'evidence'};
				my @keys = split(/ /, $evidence);
				
				my @sources = ();
				my @evid = ();
				# Old routine
				# foreach $evidence ($x -> children('evidence'))
				# {
				# 	if ((contains($evidence->{'att'}->{'key'}, @keys)))
				# 	{
				# 		if ($evidence->{'att'}->{'type'} ne 'Literature')
				# 		{
				# 			# die("Error: Evidence is not PubMed for '$name' PTM (evidence key '".$evidence->{'att'}->{'key'}."')");
				# 			next;
				# 		}
				# 		
				# 		my $source = $evidence->{'att'}->{'attribute'};
				# 		if ($source =~ /PubMed=(\d+)/)
				# 		{
				# 			push(@sources, $1);
				# 		}
				# 		else
				# 		{
				# 			# die("Error: Source is '$source' instead of format PubMed=\\d+ for '$name'");
				# 		}
				# 	}
				# }
				# New routine
				foreach $evidence ($x -> children('evidence'))
				{
					# if ($evidence -> children_count('source') == 0)
					# {
					# 	die("Error: No PubMed source for Feature '$type' evidence for '$name'")
					# 	# next;
					# }
					
					# This seems to be stored differently now:
				    # <feature type="modified residue" description="N-acetylmethionine; in 14-3-3 protein beta/alpha; alternate" evidence="4 36">
				    #   <location>
				    #     <position position="1"/>
				    #   </location>
				    # </feature>

					# Evidence: 4 and 36 are evidence "keys":
				    # <evidence type="ECO:0000244" key="4">
				    #   <source ref="4378"/>
				    # </evidence>
				    # <evidence type="ECO:0000269" key="36">
				    #   <source ref="7"/>
				    # </evidence>

					# References:
				    # <reference key="7">
				    #   <citation type="submission" date="2008-12" db="UniProtKB">
				    #     <authorList>
				    #       <person name="Bienvenut W.V."/>
				    #       <person name="Zebisch A."/>
				    #       <person name="Kolch W."/>
				    #     </authorList>
				    #   </citation>
				    #   <scope>PROTEIN SEQUENCE OF 1-11; 14-57; 63-70; 106-117; 130-169 AND 215-246</scope>
				    #   <scope>CLEAVAGE OF INITIATOR METHIONINE</scope>
				    #   <scope>ACETYLATION AT MET-1 AND THR-2</scope>
				    #   <scope>IDENTIFICATION BY MASS SPECTROMETRY</scope>
				    #   <source>
				    #     <tissue>Colon carcinoma</tissue>
				    #   </source>
				    # </reference>

					# >> This isn't a PubMed paper, though.
					
  					
					if ($evidence -> children_count('source') == 0)
					{
						# die("Error: No PubMed source for Feature '$type' evidence for '$name'");
						addme("no pubmed source for feature '$type' evidence for name", $name);
						next;
					}
					if ($evidence -> children_count('source') > 1)
					{
						die("Error: Multiple sources for evidence for '$name'");
					}
					
					if ((contains($evidence->{'att'}->{'key'}, @keys)))
					{
						push(@evid, $evidence->{'att'}->{'type'});
						
						foreach $dbref ($evidence -> first_child('source') -> children('dbReference'))
						{
							if ($dbref->{'att'}->{'type'} eq 'PubMed')
							{
								push(@sources, $dbref->{'att'}->{'id'});
							}
							else
							{
								# die("Error: Evidence is not PubMed for '$name' PTM (evidence key '".$evidence->{'att'}->{'key'}."')");
							}
						}
					}
				}
				
				# if (scalar(unique(@evid)) == 1) { @evid = unique(@evid); }
				# $evid = join("|", @evid);
				# $sources = join("|", @sources);
				$evid = join("|", unique(@evid));
				$sources = join("|", unique(@sources));
			}
			
		    die("Error: Location duplicates for PTM of '$name'") if ($feature->children_count('location') != 1);
		    if ($feature->first_child('location')->children_count('position') != 1)
			{
				# this happens e.g. for internal cross-links (e.g. 970 -> 973 Cys-Gln crosslink in A2ML1_HUMAN)
				# die("Error: Location not a single amino acid for PTM '$type' of '$name'");
				next;
			}
			
			# my $description = $feature->{'att'}->{'type'}." ".$feature->{'att'}->{'description'};
			my $description = $feature->{'att'}->{'description'};
			my $site = $feature->first_child('location')->first_child('position')->{'att'}->{'position'};

			# Check if there is isoform information (in the location tag):
		    # <feature type="modified residue" description="Phosphothreonine" evidence="9">
		    #   <location sequence="Q2V2M9-4">
		    #     <position position="1476"/>
		    #   </location>
		    # </feature>
			my $isoacc = $primary_acc;
			if ($feature->children_count('location') > 1)
			{
				die("Error: Multiple locations for modified residue feature for '$name'")
			}
			if ($feature->children_count('location') == 1)
			{
				if ($feature->first_child('location')->att_exists('sequence'))
				{
					# Isoform acc (e.g. Q2V2M9-4)
					$isoacc = $feature->first_child('location')->{'att'}->{'sequence'};
					if (!defined($isoacc))
					{
						die("Error: Couldn't parse isoform acc for modified residue feature for '$name'")
					}
				}
			}
			
			# Check if there is 'status' information (e.g. "By similarity")
			my $status = '';
			if ($feature->att_exists('status'))
			{
				$status = $feature->{'att'}->{'status'};
			}
			
			# Scale
			# <reference key="10">
			# 		    <citation date="2009" first="834" last="840" name="Science" type="journal article" volume="325">
			# 		      <title>Lysine acetylation targets protein complexes and co-regulates major cellular functions.</title>
			# 		      <authorList>
			# 		        <person name="Choudhary C."/>
			# 		        <person name="Kumar C."/>
			# 		        <person name="Gnad F."/>
			# 		        <person name="Nielsen M.L."/>
			# 		        <person name="Rehman M."/>
			# 		        <person name="Walther T."/>
			# 		        <person name="Olsen J.V."/>
			# 		        <person name="Mann M."/>
			# 		      </authorList>
			# 		      <dbReference id="10.1126/science.1175371" key="11" type="DOI"/>
			# 		      <dbReference id="19608861" key="12" type="PubMed"/>
			# 		    </citation>
			# 		    <scope>ACETYLATION [LARGE SCALE ANALYSIS] AT LYS-248</scope>
			# 		    <scope>MASS SPECTROMETRY</scope>
			# 		  </reference>
		  	# @sources >> PMID >> reference >> scope >> scale large/small
			my $scale = '';
			my $studies = 0;
			my $studies_small = 0;
			my $studies_large = 0;
			foreach $source (split(/\|/, $sources))
			{
				foreach $reference ($x -> children('reference'))
				{
					foreach $citation ($reference -> children('citation'))
					{
						foreach $dbref ($citation -> children('dbReference'))
						{
							if ($dbref->{'att'}->{'type'} eq 'PubMed')
							{
								if ($dbref->{'att'}->{'id'} eq $source)
								{
									foreach $scope ($reference -> children_text('scope'))
									{
										# Canonical: Use sequence from primary accession
										my $tmpseq = $seq;
										if ($isoacc ne $primary_acc)
										{
											# Isoform: Use sequence from isoform accession
											$tmpseq = $varsplic{$isoacc};
										}
										if ($site > length($tmpseq))
										{
											# This happens if the annotation is for an isoform longer than the canonical one
											warn("\n\nWarning: Site is outside of sequence (it must be for an isoform longer than the canonical one) (skipped adding to unimod) (kept rest) for:\n\nname='$name', isoacc='$isoacc', species='$species', site='$site', ptm='', type='$type', status='$status', description='".esc($description)."', scale='$scale', evidence='$evid', sources='$sources'\n\n");
											next;
										}
										# print "INSERT INTO unimod SET name='$name', species='$species', site='$site', status='$status', description='".esc($description)."', evidence='$evid', sources='$sources'\n";
										# print "SCOPE='$scope'\n";
										$aa = substr($tmpseq, $site - 1, 1);
										# print "AA #$site in sequence = $aa";
										$aa = onetothree($aa);
										# print " >> $aa-$site\n";
										if ($scope =~ /$aa-$site/)
										{
											if ($scope =~ /\[LARGE SCALE ANALYSIS\]/)
											{
												if ($scale eq 'small')
												{
													$scale = 'both';
												}
												else
												{
													$scale = 'large';
												}
												$studies++;
												$studies_large++;
											}
											else
											{
												if ($scale eq 'large')
												{
													$scale = 'both';
												}
												else
												{
													$scale = 'small';
												}
												$studies++;
												$studies_small++;
											}
										}
										# print "SCALE='$scale'\n";
									}
								}
							}
						}
					}
				}
			}
			
			# Update unimod (Name -> PTMs table)
			$q = "INSERT INTO unimod SET name='$name', acc='$isoacc', canon='$primary_acc', species='$species', site='$site', ptm='', type='$type', description='".esc($description)."', source='UniProt', subset='', scale='$scale', evid='$evid', studies='$studies', studies_small='$studies_small', studies_large='$studies_large', pmids='$sources'";
			$q =~ s/=''/=NULL/g;
			Query($q) unless switch('debug');
		}
		
		# Example:
		# <feature description="(in allele A*29:04)" id="VAR_016347" type="sequence variant">
		# <original>N</original>
		# <variation>H</variation>
		# <location>
		# <position position="90"/>
		# </location>
		# <location>
		# <begin position="121"/>
		# <end position="122"/>
		# </location>
		# </feature>
		elsif ($feature->{'att'}->{'type'} eq 'sequence variant')
		{
			# Check if there is 'evidence' information
			my $evid = '';
			my $sources = '';
			if ($feature->att_exists('evidence'))
			{
				my $evidence = $feature->{'att'}->{'evidence'};
				my @keys = split(/ /, $evidence);
				
				my @sources = ();
				my @evid = ();
				# Old routine
				# foreach $evidence ($x -> children('evidence'))
				# {
				# 	if ((contains($evidence->{'att'}->{'key'}, @keys)))
				# 	{
				# 		die("Error: Evidence is not PubMed for '$name' PTM (evidence key '".$evidence->{'att'}->{'key'}."')") if ($evidence->{'att'}->{'type'} ne 'Literature');
				# 		
				# 		my $source = $evidence->{'att'}->{'attribute'};
				# 		if ($source =~ /PubMed=(\d+)/)
				# 		{
				# 			push(@sources, $1);
				# 		}
				# 		else
				# 		{
				# 			# die("Error: Source is '$source' instead of format PubMed=\\d+ for '$name'");
				# 		}
				# 	}
				# }
				# New routine
				foreach $evidence ($x -> children('evidence'))
				{
					if ($evidence -> children_count('source') == 0)
					{
						# die("Error: No PubMed source for 'sequence variant' Feature evidence for '$name'")
						# # I think I'll need to improve this parsing so it goes feature -> evidence -> reference to get PubMed IDs and to make "scale small/large" work again
						# Don't need study scale for sequence variants, though
						next;
					}
					if ($evidence -> children_count('source') > 1)
					{
						die("Error: Multiple sources for evidence for '$name'");
					}
					
					if ((contains($evidence->{'att'}->{'key'}, @keys)))
					{
						push(@evid, $evidence->{'att'}->{'type'});
						
						foreach $dbref ($evidence -> first_child('source') -> children('dbReference'))
						{
							if ($dbref->{'att'}->{'type'} eq 'PubMed')
							{
								push(@sources, $dbref->{'att'}->{'id'});
							}
							else
							{
								# die("Error: Evidence is not PubMed for '$name' PTM (evidence key '".$evidence->{'att'}->{'key'}."')");
							}
						}
					}
				}
				
				# if (scalar(unique(@evid)) == 1) { @evid = unique(@evid); }
				# $evid = join("|", @evid);
				# $sources = join("|", @sources);
				$evid = join("|", unique(@evid));
				$sources = join("|", unique(@sources));
			}

			if (($feature->children_count('original') == 0) or ($feature->children_count('variation') == 0))
			{
				# skip if SNP has no sequence information
				next;
			}
			
			die("Error: Location duplicates for SNP of '$name' (description '$description')") if ($feature->children_count('location') != 1);
			die("Error: Original duplicates for SNP of '$name' (description '$description')") if ($feature->children_count('original') != 1);
			die("Error: Variation duplicates for SNP of '$name' (description '$description')") if ($feature->children_count('variation') != 1);
			# die("Error: Location not a single amino acid for SNP of '$name'") if ($feature->first_child('location')->children_count('position') != 1);

			my $description = '';
			my $original = $feature->first_child_text('original');
			my $variant = $feature->first_child_text('variation');
			my $start = '';
			my $stop = '';
			
			# Check if there is a description (i.e. the "note" column in unisnp)
			if ($feature->att_exists('description'))
			{
				my $description = $feature->{'att'}->{'description'};
			}
			
			if ($feature->first_child('location')->children_count('position') == 1)
			{
				# SNP at only one residue
				$start = $feature->first_child('location')->first_child('position')->{'att'}->{'position'};
				$stop = $start;
				die("Error: Length mismatch for '$name' (original '$original', variant '$variant', position '$start')") if (length($original) != 1);
			}
			else
			{
				die("Error: Counts wrong") if (($feature->first_child('location')->children_count('begin') != 1) or ($feature->first_child('location')->children_count('end') != 1));
				
				# SNP over multiple residues
				$start = $feature->first_child('location')->first_child('begin')->{'att'}->{'position'};
				$stop = $feature->first_child('location')->first_child('end')->{'att'}->{'position'};
				die("Error: Length mismatch for '$name' (original '$original', variant '$variant', start '$start', stop '$stop')") if (length($original) != ($stop - $start + 1));
			}
			
			# Check if there is isoform information (in the location tag):
			my $isoacc = $primary_acc;
			if ($feature->children_count('location') > 1)
			{
				die("Error: Multiple locations for misc feature for '$name'")
			}
			if ($feature->children_count('location') == 1)
			{
				if ($feature->first_child('location')->att_exists('sequence'))
				{
					# Isoform acc (e.g. Q2V2M9-4)
					$isoacc = $feature->first_child('location')->{'att'}->{'sequence'};
					if (!defined($isoacc))
					{
						die("Error: Couldn't parse isoform acc for sequence variant feature for '$name'")
					}
				}
			}

			# Check if there is 'status' information (e.g. "By similarity")
			my $status = '';
			if ($feature->att_exists('status'))
			{
				$status = $feature->{'att'}->{'status'};
			}
			
			# Update unisnp (Name -> SNPs table)
			# print "INSERT INTO unisnp SET name='$name', species='$species', start='$start', stop='$stop', original='$original', variant='$variant', note='".esc($description)."', evidence='$evid', sources='$sources'\n";
			# Query("INSERT INTO unisnp SET name='$name', species='$species', start='$start', stop='$stop', original='$original', variant='$variant', note='".esc($description)."', evidence='$evid', sources='$sources'") unless switch('debug');
			$q = "INSERT INTO unisnp SET name='$name', acc='$isoacc', canon='$primary_acc', species='$species', start='$start', stop='$stop', original='$original', variant='$variant', note='".esc($description)."', evidence='$evid', sources='$sources'";
			$q =~ s/=''/=NULL/g;
			Query($q) unless switch('debug');
		}
		
		# Misc features (unifeat) - domain, short sequence motif, signal peptide, transit peptide, binding site, ... (all included now)
		# 
		# <feature description="Internalization signal" status="by similarity" type="short sequence motif">     
		#   <location>
		#     <begin position="332"/>   
		#     <end position="335"/>
		#   </location>
		# </feature>
		# 
		# <feature description="Nuclear localization signal" evidence="EC1" type="short sequence motif">
		#   <location>
		#     <begin position="100"/>
		#     <end position="114"/>
		#   </location>
		# </feature>
		# 
		# elsif (contains($feature->{'att'}->{'type'}, ('short sequence motif', 'signal peptide', 'transit peptide', 'binding site', 'initiator methionine', 'propeptide', 'peptide', 'chain')))
		# elsif (contains($feature->{'att'}->{'type'}, ('short sequence motif', 'signal peptide', 'transit peptide', 'binding site', 'initiator methionine', 'propeptide', 'peptide', 'chain')))
		else
		{
			my $type = $feature->{'att'}->{'type'};

			# Check if there is 'evidence' information
			my $evid = '';
			my $sources = '';
			if ($feature->att_exists('evidence'))
			{
				my $evidence = $feature->{'att'}->{'evidence'};
				my @keys = split(/ /, $evidence);
				
				my @sources = ();
				my @evid = ();
				# Old routine
				# foreach $evidence ($x -> children('evidence'))
				# {
				# 	if ((contains($evidence->{'att'}->{'key'}, @keys)))
				# 	{
				# 		if ($evidence->{'att'}->{'type'} ne 'Literature')
				# 		{
				# 			# die("Error: Evidence is not PubMed for '$name' PTM (evidence key '".$evidence->{'att'}->{'key'}."')");
				# 			next;
				# 		}
				# 		
				# 		my $source = $evidence->{'att'}->{'attribute'};
				# 		if ($source =~ /PubMed=(\d+)/)
				# 		{
				# 			push(@sources, $1);
				# 		}
				# 		else
				# 		{
				# 			# die("Error: Source is '$source' instead of format PubMed=\\d+ for '$name'");
				# 		}
				# 	}
				# }
				# New routine
				foreach $evidence ($x -> children('evidence'))
				{
					if ($evidence -> children_count('source') == 0)
					{
						# die("Error: No PubMed source for evidence for '$name'")
						next;
					}
					if ($evidence -> children_count('source') > 1)
					{
						die("Error: Multiple sources for evidence for '$name'");
					}
					
					if ((contains($evidence->{'att'}->{'key'}, @keys)))
					{
						push(@evid, $evidence->{'att'}->{'type'});
						
						foreach $dbref ($evidence -> first_child('source') -> children('dbReference'))
						{
							if ($dbref->{'att'}->{'type'} eq 'PubMed')
							{
								push(@sources, $dbref->{'att'}->{'id'});
							}
							else
							{
								# die("Error: Evidence is not PubMed for '$name' PTM (evidence key '".$evidence->{'att'}->{'key'}."')");
							}
						}
					}
				}
				
				# if (scalar(unique(@evid)) == 1) { @evid = unique(@evid); }
				# $evid = join("|", @evid);
				# $sources = join("|", @sources);
				$evid = join("|", unique(@evid));
				$sources = join("|", unique(@sources));
			}
			
			my $start = '';
			my $stop = '';
			
		    die("Error: Location duplicates for Feature '$type' of '$name'") if ($feature->children_count('location') != 1);
			if ($feature->first_child('location')->children_count('position') == 1)
			{
				# Feature at only one residue
				$start = $feature->first_child('location')->first_child('position')->{'att'}->{'position'};
				$stop = $start;
			}
			else
			{
				die("Error: Counts wrong") if (($feature->first_child('location')->children_count('begin') != 1) or ($feature->first_child('location')->children_count('end') != 1));
				
				# Feature over multiple residues
				$start = $feature->first_child('location')->first_child('begin')->{'att'}->{'position'};
				$stop = $feature->first_child('location')->first_child('end')->{'att'}->{'position'};
			}
			
			if (!defined($start) or !defined($stop) or ($start eq '') or ($stop eq ''))
			{
				# "position -> position" for e.g. variants (single-residue features)
				if (defined($feature->first_child('location')) and 
				defined($feature->first_child('location')->first_child('position')) and
				defined($feature->first_child('location')->first_child('position')->{'att'}->{'position'}))
				{
					$start = $feature->first_child('location')->first_child('position')->{'att'}->{'position'};
					$stop = $start;
				}
				
				# location information incomplete
				# die("Error: Location information incomplete for Feature '$type' of '$name'");
				# next;
				
				if (!defined($start))
				{
					$start = '';
				}
				if (!defined($stop))
				{
					$stop = '';
				}
			}
			
			# my $description = $feature->{'att'}->{'type'}." ".$feature->{'att'}->{'description'};
			my $description = '';
			if ($feature->{'att'}->{'description'})
			{
				$description = $feature->{'att'}->{'description'};
			}
			
			# Check if there is isoform information (in the location tag):
			my $isoacc = $primary_acc;
			if ($feature->children_count('location') > 1)
			{
				die("Error: Multiple locations for misc feature for '$name'")
			}
			if ($feature->children_count('location') == 1)
			{
				if ($feature->first_child('location')->att_exists('sequence'))
				{
					# Isoform acc (e.g. Q2V2M9-4)
					$isoacc = $feature->first_child('location')->{'att'}->{'sequence'};
					if (!defined($isoacc))
					{
						die("Error: Couldn't parse isoform acc for modified residue feature for '$name'")
					}
				}
			}

			# Check if there is 'status' information (e.g. "By similarity")
			my $status = '';
			if ($feature->att_exists('status'))
			{
				$status = $feature->{'att'}->{'status'};
			}
			
			# No scale information really needed for miscellaneous features
			# # Scale
			# # <reference key="10">
			# # 		    <citation date="2009" first="834" last="840" name="Science" type="journal article" volume="325">
			# # 		      <title>Lysine acetylation targets protein complexes and co-regulates major cellular functions.</title>
			# # 		      <authorList>
			# # 		        <person name="Choudhary C."/>
			# # 		        <person name="Kumar C."/>
			# # 		        <person name="Gnad F."/>
			# # 		        <person name="Nielsen M.L."/>
			# # 		        <person name="Rehman M."/>
			# # 		        <person name="Walther T."/>
			# # 		        <person name="Olsen J.V."/>
			# # 		        <person name="Mann M."/>
			# # 		      </authorList>
			# # 		      <dbReference id="10.1126/science.1175371" key="11" type="DOI"/>
			# # 		      <dbReference id="19608861" key="12" type="PubMed"/>
			# # 		    </citation>
			# # 		    <scope>ACETYLATION [LARGE SCALE ANALYSIS] AT LYS-248</scope>
			# # 		    <scope>MASS SPECTROMETRY</scope>
			# # 		  </reference>
			# 		  	# @sources >> PMID >> reference >> scope >> scale large/small
			my $scale = '';
			# if ($feature->{'att'}->{'type'}, ('short sequence motif', 'signal peptide', 'transit peptide', 'binding site', 'initiator methionine', 'propeptide', 'peptide', 'chain')))
			# foreach $source (split(/\|/, $sources))
			# {
			# 	foreach $reference ($x -> children('reference'))
			# 	{
			# 		foreach $citation ($reference -> children('citation'))
			# 		{
			# 			foreach $dbref ($citation -> children('dbReference'))
			# 			{
			# 				if ($dbref->{'att'}->{'type'} eq 'PubMed')
			# 				{
			# 					if ($dbref->{'att'}->{'id'} eq $source)
			# 					{
			# 						foreach $scope ($reference -> children_text('scope'))
			# 						{
			# 							if (defined $start)
			# 							if (($start > length($seq)) or ($stop > length($seq)))
			# 							{
			# 								# This happens if the annotation is for an isoform longer than the canonical one
			# 								warn("\n\nWarning: Start or Stop is outside of sequence (it must be for an isoform longer than the canonical one) (skipped adding to unifeat) (kept rest) for:\n\nname='$name', species='$species', start='$start', stop='$stop', type='$type', status='$status', description='".esc($description)."', scale='$scale', evidence='$evid', sources='$sources'\n\n");
			# 								next;
			# 							}
			# 							# print "INSERT INTO unimod SET name='$name', species='$species', start='$start', stop='$stop', status='$status', description='".esc($description)."', sources='$sources'\n";
			# 							# print "SCOPE='$scope'\n";
			# 							$aa = substr($seq, $start - 1, 1);
			# 							# print "AA #$start in sequence = $aa";
			# 							$aa = onetothree($aa);
			# 							# print " >> $aa-$start\n";
			# 							if ($scope =~ /$aa-$start/)
			# 							{
			# 								if ($scope =~ /\[LARGE SCALE ANALYSIS\]/)
			# 								{
			# 									if ($scale eq 'small')
			# 									{
			# 										$scale = 'both';
			# 									}
			# 									else
			# 									{
			# 										$scale = 'large';
			# 									}
			# 								}
			# 								else
			# 								{
			# 									if ($scale eq 'large')
			# 									{
			# 										$scale = 'both';
			# 									}
			# 									else
			# 									{
			# 										$scale = 'small';
			# 									}
			# 								}
			# 							}
			# 							# print "SCALE='$scale'\n";
			# 						}
			# 					}
			# 				}
			# 			}
			# 		}
			# 	}
			# }
			
			# Update unifeat (Name -> misc features table)
			$q = "INSERT INTO unifeat SET name='$name', acc='$isoacc', canon='$primary_acc', species='$species', start='$start', stop='$stop', type='$type', description='".esc($description)."', scale='$scale', evidence='$evid', sources='$sources'";
			$q =~ s/=''/=NULL/g;
			Query($q) unless switch('debug');
		}
	}
	
	# Output entry (for debugging)
	# print "\n\n";
	
	# $x -> flush;            # outputs the section and frees memory
	$x -> purge;              # frees memory

	# print "\n\n";
	# 
	# print "ACC=['".join("|", @acc)."']\n";
	# print "NAME=['$name']\n";
	# print "SPECIES=['$species']\n";
	# 
	# print "\n";

	# exit;
	# if (getme() > 20)
	# {
	# 	showme("dbReference types");
	# 	showme("Comment types");
	# 	showme("Feature types");
	# 	exit;
	# }

	stepme(100);
}

stopme();
stoptime();

$xml -> purge;

# showme("dbReference types");
# showme("Comment types");
# showme("Feature types");

# showmeallsorted(1);
showmeallsorted();

done();
