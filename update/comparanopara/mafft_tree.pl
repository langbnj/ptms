#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$scale = '';            # Disable newick2mafft branch length scaling
# $scale = '0.1 ';        # Enable newick2mafft branch length scaling (for the topology tree with all branches=1)
# $scale = '0.001 ';    # Enable newick2mafft branch length scaling (for MYA branch lengths from timetree, max is going to be 1.226)
# $scale = '0.000815660685154976 ';      # Enable newick2mafft branch length scaling (for MYA branch lengths from timetree, scaled to max 1)

our $usage = "$0 [species] [MAFFT mode] [-1para] [-direct]\n\n Species: Reference species.\n MAFFT mode: linsi, einsi or ginsi (for L-INS-i, E-INS-i or G-INS-i)\n -1para: Keeps the best-matching outparalog (by various metrics) for one2many and many2many cases.\n -direct: Run MAFFT right here in this process instead of submitting jobs to the cluster.\n\nExample: $0 human linsi";
($unispec, $mafftmode) = args(2);

# $treefile = "../compara/input/species_tree.nw";
# $treefile = "../compara/input/species_tree.branch_len.reordered_in_dendroscope.nw";
# $treefile = "../compara/input/species_tree.branch_len.reordered_in_dendroscope.yeastrooted.taxids.internal_nodes_labelled.nw";
$treefile = "../compara/input/species_tree_final.nw";

# Handling for 1para
$label = '';
# $labeldot = '';
# $labelslash = '';
if (switch('1parahc'))
{
	$label = '1parahc';
}
elsif (switch('1para'))
{
	$label = '1para';
}
# $labeldot .= '.' if ($label ne '');
# $labelslash .= '/' if ($label ne '');
# if ($label ne '')
# {
# 	run("Create tmp directory", "mkdir -p tmp/$label", 1);
# 	run("Create output directory", "mkdir -p output/$label", 1);
# }



# start

# Get Compara species tree

# Read from file
startme("Reading Compara species tree from '$treefile'");
open(TREE, $treefile) or die("Error: Couldn't open '$treefile'");
$speciestree = '';
while (<TREE>)
{
	$speciestree .= $_;
	stepme(10);
}
stopme();
close(TREE);

state("Species guide tree is '$speciestree'", 1) if (switch('debug'));

# Feeding in the processed tree just for checking:
# $speciestree = "(homo_sapiens:0.006267,((pan_paniscus:0.003330,pan_troglodytes:0.002210):0.004314,(gorilla_gorilla:0.008565,(pongo_abelii:0.016070,(nomascus_leucogenys:0.019679,(((rhinopithecus_roxellana:0.002093,rhinopithecus_bieti:0.003017):0.014320,(chlorocebus_sabaeus:0.011709,((macaca_nemestrina:0.004313,(macaca_mulatta:0.002549,macaca_fascicularis:0.003121):0.000632):0.004404,(papio_anubis:0.006543,(cercocebus_atys:0.005684,mandrillus_leucophaeus:0.006766):0.000577):0.001190):0.003970):0.004770):0.013696,((aotus_nancymaae:0.025342,(callithrix_jacchus:0.029787,(saimiri_boliviensis_boliviensis:0.026290,cebus_capucinus:0.026040):0.002968):0.000100):0.023606,(carlito_syrichta:0.079756,((otolemur_garnettii:0.074915,(microcebus_murinus:0.041368,(prolemur_simus:0.037167,propithecus_coquereli:0.038183):0.005197):0.027371):0.007201,((tupaia_belangeri:0.100279,((oryctolagus_cuniculus:0.077519,ochotona_princeps:0.091231):0.025426,(((octodon_degus:0.076196,(heterocephalus_glaber_female:0.072806,(chinchilla_lanigera:0.069808,cavia_porcellus:0.073542):0.003089):0.002173):0.021884,(sciurus_vulgaris:0.059890,(marmota_marmota_marmota:0.018723,(ictidomys_tridecemlineatus:0.013839,urocitellus_parryii:0.012741):0.006042):0.039685):0.030911):0.007277,(dipodomys_ordii:0.100796,(jaculus_jaculus:0.095059,(nannospalax_galili:0.088497,((microtus_ochrogaster:0.070130,(peromyscus_maniculatus_bairdii:0.064003,(cricetulus_griseus_chok1gshd:0.048927,mesocricetus_auratus:0.052393):0.016237):0.003158):0.010924,(rattus_norvegicus:0.062135,(mus_pahari:0.036579,(mus_caroli:0.020455,((mus_spretus:0.008996,mus_spicilegus:0.008514):0.001096,mus_musculus:0.010292):0.010886):0.016872):0.023114):0.016909):0.012163):0.009640):0.007501):0.000100):0.003594):0.001114):0.000100,(((erinaceus_europaeus:0.100896,sorex_araneus:0.107814):0.000100,(((felis_catus:0.011671,(panthera_tigris_altaica:0.004625,(panthera_leo:0.002369,panthera_pardus:0.002761):0.001185):0.007534):0.053232,((vulpes_vulpes:0.012223,(canis_lupus_familiaris:0.001075,canis_lupus_dingo:0.001295):0.010237):0.047311,((mustela_putorius_furo:0.016083,neovison_vison:0.016027):0.038935,(ailuropoda_melanoleuca:0.017262,(ursus_americanus:0.003588,ursus_maritimus:0.003302):0.014023):0.037337):0.006516):0.006228):0.016495,(((equus_caballus:0.006315,equus_asinus_asinus:0.006725):0.074616,(pteropus_vampyrus:0.076676,(rhinolophus_ferrumequinum:0.077482,myotis_lucifugus:0.081578):0.000100):0.004167):0.004614,((vicugna_pacos:0.020045,camelus_dromedarius:0.016175):0.060790,(((balaenoptera_musculus:0.018897,(physeter_catodon:0.020171,(tursiops_truncatus:0.013185,(phocoena_sinus:0.008454,(delphinapterus_leucas:0.003326,monodon_monoceros:0.002784):0.004886):0.004069):0.009181):0.000036):0.042769,(cervus_hanglu_yarkandensis:0.032005,(moschus_moschiferus:0.032769,((capra_hircus:0.011454,ovis_aries_rambouillet:0.010946):0.016875,(bison_bison_bison:0.004893,((bos_indicus_hybrid:0.002089,bos_taurus:0.002301):0.002516,(bos_mutus:0.001643,bos_grunniens:0.001797):0.003876):0.000100):0.022635):0.003043):0.002208):0.038168):0.010258,(catagonus_wagneri:0.047438,sus_scrofa:0.046350):0.027979):0.000100):0.009003):0.002140):0.016799):0.000100,((dasypus_novemcinctus:0.077661,choloepus_hoffmanni:0.078289):0.026667,((echinops_telfairi:0.096657,(loxodonta_africana:0.062164,procavia_capensis:0.076396):0.016773):0.023873,((monodelphis_domestica:0.074711,(sarcophilus_harrisii:0.066306,(notamacropus_eugenii:0.063073,(phascolarctos_cinereus:0.037496,vombatus_ursinus:0.035414):0.021487):0.012155):0.002102):0.041817,(ornithorhynchus_anatinus:0.126661,(((sphenodon_punctatus:0.128683,((podarcis_muralis:0.104185,salvator_merianae:0.104005):0.013026,(anolis_carolinensis:0.118737,(naja_naja:0.032020,(laticauda_laticaudata:0.026340,(notechis_scutatus:0.013853,pseudonaja_textilis:0.014906):0.008460):0.004870):0.083651):0.000100):0.012129):0.000665,((pelodiscus_sinensis:0.078041,((gopherus_evgoodei:0.022575,chelonoidis_abingdonii:0.022615):0.012439,(terrapene_carolina_triunguis:0.012521,chrysemys_picta_bellii:0.016049):0.017953):0.041302):0.040480,(crocodylus_porosus:0.112610,(struthio_camelus_australis:0.086911,(((anser_brachyrhynchus:0.031580,anas_platyrhynchos_platyrhynchos:0.028350):0.040115,(coturnix_japonica:0.044271,(gallus_gallus:0.038421,meleagris_gallopavo:0.040999):0.005829):0.031485):0.006220,((strigops_habroptila:0.058819,aquila_chrysaetos_chrysaetos:0.054101):0.017982,(ficedula_albicollis:0.045381,(parus_major:0.045098,(geospiza_fortis:0.037144,(serinus_canaria:0.040671,taeniopygia_guttata:0.041109):0.000100):0.011726):0.000986):0.026323):0.013547):0.001014):0.022625):0.001615):0.018651):0.006069,((xenopus_tropicalis:0.134210,leptobrachium_leishanense:0.119760):0.010564,(latimeria_chalumnae:0.119616,((erpetoichthys_calabaricus:0.131615,(lepisosteus_oculatus:0.157564,((paramormyrops_kingsleyae:0.119026,scleropages_formosus:0.114334):0.015702,(((denticeps_clupeoides:0.120175,clupea_harengus:0.101885):0.009364,((ictalurus_punctatus:0.100599,(electrophorus_electricus:0.107845,(astyanax_mexicanus:0.083273,pygocentrus_nattereri:0.089737):0.015239):0.004423):0.009471,(danio_rerio:0.078225,(carassius_auratus:0.041732,(sinocyclocheilus_grahami:0.038219,cyprinus_carpio_carpio:0.033091):0.008693):0.035893):0.032398):0.000399):0.000838,((esox_lucius:0.100082,(hucho_hucho:0.033785,((salmo_salar:0.010539,salmo_trutta:0.012411):0.017115,(oncorhynchus_mykiss:0.016195,(oncorhynchus_kisutch:0.014576,oncorhynchus_tshawytscha:0.014744):0.003275):0.012553):0.001002):0.055176):0.032039,(gadus_morhua:0.112842,(myripristis_murdjan:0.086403,(hippocampus_comes:0.114595,(((tetraodon_nigroviridis:0.069559,takifugu_rubripes:0.072351):0.037901,(labrus_bergylta:0.083828,((sparus_aurata:0.068108,(larimichthys_crocea:0.060969,dicentrarchus_labrax:0.062831):0.006507):0.006150,((sander_lucioperca:0.065758,cottoperca_gobio:0.064842):0.009602,(gasterosteus_aculeatus:0.086144,cyclopterus_lumpus:0.066556):0.006873):0.000498):0.007561):0.009147):0.000389,(((mastacembelus_armatus:0.080065,(betta_splendens:0.084898,anabas_testudineus:0.072182):0.008720):0.007410,((cynoglossus_semilaevis:0.096693,scophthalmus_maximus:0.098807):0.000100,(lates_calcarifer:0.060256,(seriola_lalandi_dorsalis:0.012786,seriola_dumerili:0.012774):0.043059):0.025646):0.000100):0.012554,((((oryzias_melastigma:0.044789,oryzias_javanicus:0.037601):0.017143,(oryzias_sinensis:0.020479,oryzias_latipes:0.018281):0.044962):0.049551,((nothobranchius_furzeri:0.092923,kryptolebias_marmoratus:0.088457):0.008070,((fundulus_heteroclitus:0.080179,cyprinodon_variegatus:0.083961):0.000046,(xiphophorus_maculatus:0.040289,(poecilia_reticulata:0.028708,(poecilia_latipinna:0.004392,poecilia_formosa:0.005228):0.018772):0.008664):0.044750):0.017139):0.009730):0.000100,((stegastes_partitus:0.051930,(acanthochromis_polyacanthus:0.039896,(amphiprion_ocellaris:0.006772,amphiprion_percula:0.005508):0.027059):0.021052):0.026953,(amphilophus_citrinellus:0.055699,(oreochromis_niloticus:0.022404,(neolamprologus_brichardi:0.015663,(haplochromis_burtoni:0.007454,(pundamilia_nyererei:0.008829,(astatotilapia_calliptera:0.002402,maylandia_zebra:0.002198):0.000566):0.000673):0.008121):0.011566):0.037173):0.021495):0.017832):0.000100):0.000100):0.018663):0.000100):0.013949):0.005874):0.004667):0.000100):0.006005):0.006235):0.000100,(callorhinchus_milii:0.145499,((eptatretus_burgeri:0.128661,petromyzon_marinus:0.168899):0.000100,((ciona_intestinalis:0.151626,ciona_savignyi:0.154164):0.032643,((drosophila_melanogaster:0.154675,caenorhabditis_elegans:0.188645):0.000100,saccharomyces_cerevisiae:0.241236):0.009979):0.042857):0.004591):0.000097):0.017450):0.000100):0.012491):0.002583):0.006190):0.012739):0.004473):0.000100):0.000100):0.021428):0.000041):0.027719):0.017662):0.010949):0.002677):0.008475):0.001864):0.000330);";

# # Old: reliant on manually processed tree
# # # Process tree (no longer necessary now, already processed) (see functions.inc.pl for what comparaprune did)
# # # state($speciestree);
# # # showbranchlengths($speciestree);
# # # $speciestree = comparaprune($speciestree, 1);
# # # state($speciestree);
# # # showbranchlengths($speciestree);
# # @species = returnterminals($speciestree);

# # Show species tree pre-pruning
# state("PRE Species tree, pre-pruning:");
# state($speciestree);
# state("PRE Species tree branch lengths, pre-pruning:");
# showbranchlengths($speciestree);

# Prune species tree
# $speciestree = comparaprune($speciestree);
$speciestree = comparaprune($speciestree, 1);

# # Show species tree post-pruning
# state("POST Species tree, post-pruning:");
# state($speciestree);
# state("POST Species tree branch lengths, post-pruning:");
# showbranchlengths($speciestree);
@species = returnterminals($speciestree);
# state("POST Species 'terminals':");
# show(@species);

# exit;


# Start aligning

# Translate UniProt species mnemonic to full Ensembl species name
$query = Query("SELECT species FROM comparaspecies WHERE unispec='$unispec'");
($species) = FetchOne($query);

if (!switch('1para') and !switch('1parahc'))
{
	# normal
	# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e, unimod m WHERE c.species='$species' AND e.ensp=c.ensp AND m.name=e.name AND m.ptm!='' GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
	# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e WHERE c.species='$species' AND e.ensp=c.ensp GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln LIMIT 1000");
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
startme("Aligning one2one orthologs for '$unispec' ($species) using MAFFT, while giving a species tree", 0, Numrows($mainquery));
starttime();

$free = 0;
while (($ensp) = Fetch($mainquery))
{
	stepme(1);

	print(" >> $ensp\n") if (switch('debug'));

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
	
	# if input FASTA file doesn't exist or is zero
	if (!-s $infile)
	{
		addme("couldn't find FASTA file for ensp (probably because it doesn't have any one-to-one orthologs (only one-to-many) in the orthologue mart files) (skipped)", $ensp);
		next;
	}
	
	# Set MAFFT output file names
	# $outfile = "output/mafft.$ensp.$mafftmode.tree.txt";
	if (switch('1para') or switch('1parahc'))
	{
		$outfile = "output/mafft.$ensp.$mafftmode.tree.$label.txt";
	}
	else
	{
		$outfile = "output/mafft.$ensp.$mafftmode.tree.txt";
	}

    if (-s $outfile)
    {
		# Skip this ensp if output file already exists and is non-zero
        addme("alignment file already exists for ensp (skipped)", $ensp);
        next;
		# # Crash instead
		# die("Error: Alignment file '$outfile' already exists for ensp '$ensp' (should run clean.pl first)");
    }

	# Build input tree from Compara species tree
	# Need to do this for each alignment individually since MAFFT wants the tree leaves to have numbers as names, in the same order as the input sequences
	# See also http://mafft.cbrc.jp/alignment/software/treein.html, http://mafft.cbrc.jp/alignment/software/newick2mafft.rb
	$tree = $speciestree;

	# # Get ENSPs in this alignment from table 'comparahomology'
	# @titles = ();
	# push(@titles, $ensp);
	# if (switch('1para') or switch('1parahc'))
	# {
	# 	# 1para
	# 	$query = Query("SELECT DISTINCT ensp2 FROM comparahomology WHERE species1='$unispec' AND ensp1='$ensp' AND ensp2 IS NOT NULL ORDER BY ensp2");
	# }
	# else
	# {
	# 	# normal
	# 	$query = Query("SELECT DISTINCT ensp2 FROM comparahomology WHERE species1='$unispec' AND ensp1='$ensp' AND homology='ortholog_one2one' AND hc=1 AND ensp2 IS NOT NULL ORDER BY ensp2");
	# }
	# while (($homoensp) = Fetch($query))
	# {
	# 	push(@titles, $homoensp);
	# }

	state("Reading alignment file '$infile'", 1) if (switch('debug'));
	open(IN, "$infile") or die("Error: Couldn't open '$infile'");
	fastabreak();
	@titles = ();
	while (<IN>)
	{
		($title, $seq) = getfasta();
		$seq = undef;	# Just to suppress unused warning
		
		push(@titles, $title);
	}
	close(IN);

	# Build substitution hash (species >> number)
	state("Building species >> number substitution hash from species guide tree '$speciestree'", 1) if (switch('debug'));
	%sub = ();
	$i = 0;
	foreach $ensp (@titles)
	{
		$i++;
		$query = Query("SELECT species FROM comparaenspspecies WHERE ensp='$ensp'");
		($thisspecies) = FetchOne($query);
		# print " >> $thisspecies >> $ensp >> $i\n";
		
		if (exists($sub{$thisspecies}))
		{
			die("Error: Species '$thisspecies' already being substituted to $sub{$thisspecies}") ;
		}
		if (!contains($thisspecies, @species))
		{
		    show(@species);
            die("Error: Species '$thisspecies' not found in tree");
            # addme("species not found in ensembl guide tree for species (skipped)", $thisspecies);
            # addme("species not found in ensembl guide tree for ensp (skipped)", $ensp);
            # next;
		}
		$sub{$thisspecies} = $i;
    	state("Will substitute '$thisspecies' >> '$i'", 1) if (switch('debug'));
	}
	
	# Make substitutions in tree & remove unnecessary nodes
	state("Make substitutions in species guide tree with Bio::Phylo:\n\n$speciestree\n\n", 1) if (switch('debug'));
	$t = Bio::Phylo::IO->parse(-string => $tree, -format => 'newick')->first;
	foreach $node (@{$t->get_root->get_terminals()})
	{
		if (exists($sub{$node->get_name}))
		{
			# Substitute number for node name
			# addme("substituted from", $node->get_name);
			$node->set_name($sub{$node->get_name});
			# addme("substituted to", $node->get_name);
		}
		else
		{
			# addme("removed", $node->get_name);
			# mafftprune removes it below
		}
	}
	# showmeall();
	state("Guide tree after substitutions:\n\n$tree\n\n", 1) if (switch('debug'));

	# Process tree (see functions.inc.pl)
	$tree = $t->to_newick;
	state("Guide tree in Newick format:\n\n$tree\n\n", 1) if (switch('debug'));
	# state($tree);
	# showbranchlengths($tree);
	# showterminaldistances($tree);
	$tree = mafftprune($tree);
	state("Guide tree in Newick format, pruned:\n\n$tree\n\n", 1) if (switch('debug'));
	# state($tree);
	# showbranchlengths($tree);
	# showterminaldistances($tree);

	# Write NH tree
	if (switch('1para') or switch('1parahc'))
	{
		# 1para
		$tmpfile = "tmp/tmp-tree-$ensp-$mafftmode-$label-raw.txt";
	}
	else
	{
		# normal
		$tmpfile = "tmp/tmp-tree-$ensp-$mafftmode-raw.txt";
	}
	state("Writing tree to '$tmpfile'", 1) if (switch('debug'));
	open(TMP, ">$tmpfile") or die("Error: Couldn't open '$tmpfile'");
	print TMP $tree."\n";
	close(TMP);
	
	# Convert to MAFFT internal tree format
	if (switch('1para') or switch('1parahc'))
	{
		# 1para
		$treeinfile = "tmp/tmp-tree-$ensp-$mafftmode-$label.mafft";
	}
	else
	{
		# normal
		$treeinfile = "tmp/tmp-tree-$ensp-$mafftmode.mafft";
	}
    # system("ruby ./bin/newick2mafft.rb $tmpfile > $treeinfile");
	state("Converting this Newick tree:\n\n".`cat $tmpfile`."\n\n", 1) if (switch('debug'));
	state("ruby ./bin/newick2mafft.rb $scale$tmpfile > $treeinfile") if (switch('debug'));
	system("ruby ./bin/newick2mafft.rb $scale$tmpfile > $treeinfile");
	state("..to MAFFT format:\n\n".`cat $treeinfile`."\n\n", 1) if (switch('debug'));
	# run("", "ruby ./bin/newick2mafft.rb $tmpfile > $treeinfile");

	# Run MAFFT!
	chdir("tmp"); # make sure log files end up in ./tmp/
	
    # Wait for grid nodes to free up
    if (!switch('direct'))
    {
        while ($free <= 0)
        {
          $free = freenodes();
        }
    }

	# Run MAFFT
	# run("MAFFT L-INS-i", "~/\.sh mafft --localpair --maxiterate 16 --reorder ../$infile \\> ../$outfile", 1);
	# run("MAFFT $mafftmode", "~/scripts/qsub.sh mafft --localpair --maxiterate 1000 --anysymbol ../$infile \\> ../$outfile", 1);

	state("Running MAFFT:", 1) if (switch('debug'));
	addme("ran MAFFT for ensp", $ensp);
    if (!switch('direct'))
    {
        state("~/scripts/qsub.sh mafft-$mafftmode --quiet --anysymbol --treein ../$treeinfile ../$infile \\> ../$outfile") if (switch('debug'));
        run("MAFFT $mafftmode", "~/scripts/qsub.sh mafft-$mafftmode --quiet --anysymbol --treein ../$treeinfile ../$infile \\> ../$outfile", 1);
        $free--;
    }
    else
    {
        state("mafft-$mafftmode --quiet --anysymbol --treein ../$treeinfile ../$infile > ../$outfile") if (switch('debug'));
        if (switch('debug'))
        {
            run("MAFFT $mafftmode", "mafft-$mafftmode --anysymbol --treein ../$treeinfile ../$infile > ../$outfile");
        }
        else
        {
            run("MAFFT $mafftmode", "mafft-$mafftmode --quiet --anysymbol --treein ../$treeinfile ../$infile > ../$outfile", 1);
        }
        $free--;
    }
	
    # showmeall();
    # stopme();
    # stoptime();
    # die("DEBUG");
	
	# Done!
	chdir("..");
}
stopme();
stoptime();

# showmeall();
showmesome(50);

done();
