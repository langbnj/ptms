#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species] [MAFFT mode] [-tree]\n\n Species: Reference species.\n MAFFT mode: linsi, einsi or ginsi (for L-INS-i, E-INS-i or G-INS-i)\n".
" -tree: Use MAFFT results with species tree.\n -1para: Keeps the best-matching outparalog (by various metrics) for one2many and many2many cases.\n\n\nExample: $0 human linsi -tree";
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



# $speciestreefile = "../compara/input/species_tree.taxids.nw";

# To do taxid substitutions:
# SELECT CONCAT('s/',species,':/',tax,':/; ') FROM comparaspecies;
# cat species_tree.branch_len.reordered_in_dendroscope.nw | perl -npe "s/homo_sapiens:/9606:/; s/pan_paniscus:/9597:/; s/pan_troglodytes:/9598:/; s/gorilla_gorilla:/9595:/; s/pongo_abelii:/9601:/; s/nomascus_leucogenys:/61853:/; s/rhinopithecus_roxellana:/61622:/; s/rhinopithecus_bieti:/61621:/; s/chlorocebus_sabaeus:/60711:/; s/macaca_nemestrina:/9545:/; s/macaca_mulatta:/9544:/; s/macaca_fascicularis:/9541:/; s/papio_anubis:/9555:/; s/cercocebus_atys:/9531:/; s/mandrillus_leucophaeus:/9568:/; s/aotus_nancymaae:/37293:/; s/callithrix_jacchus:/9483:/; s/saimiri_boliviensis_boliviensis:/39432:/; s/cebus_capucinus:/2715852:/; s/carlito_syrichta:/1868482:/; s/otolemur_garnettii:/30611:/; s/microcebus_murinus:/30608:/; s/prolemur_simus:/1328070:/; s/propithecus_coquereli:/379532:/; s/tupaia_belangeri:/37347:/; s/oryctolagus_cuniculus:/9986:/; s/ochotona_princeps:/9978:/; s/octodon_degus:/10160:/; s/heterocephalus_glaber_female:/10181:/; s/chinchilla_lanigera:/34839:/; s/cavia_porcellus:/10141:/; s/sciurus_vulgaris:/55149:/; s/marmota_marmota_marmota:/9994:/; s/ictidomys_tridecemlineatus:/43179:/; s/urocitellus_parryii:/9999:/; s/dipodomys_ordii:/10020:/; s/jaculus_jaculus:/51337:/; s/nannospalax_galili:/1026970:/; s/microtus_ochrogaster:/79684:/; s/peromyscus_maniculatus_bairdii:/230844:/; s/cricetulus_griseus_chok1gshd:/10029:/; s/mesocricetus_auratus:/10036:/; s/rattus_norvegicus:/10116:/; s/mus_pahari:/10093:/; s/mus_caroli:/10089:/; s/mus_spretus:/10096:/; s/mus_spicilegus:/10103:/; s/mus_musculus:/10090:/; s/erinaceus_europaeus:/9365:/; s/sorex_araneus:/42254:/; s/felis_catus:/9685:/; s/panthera_tigris_altaica:/74533:/; s/panthera_leo:/9689:/; s/panthera_pardus:/9691:/; s/vulpes_vulpes:/9627:/; s/canis_lupus_familiaris:/9615:/; s/canis_lupus_dingo:/286419:/; s/mustela_putorius_furo:/9669:/; s/neovison_vison:/452646:/; s/ailuropoda_melanoleuca:/9646:/; s/ursus_americanus:/9643:/; s/ursus_maritimus:/29073:/; s/equus_caballus:/9796:/; s/equus_asinus_asinus:/83772:/; s/pteropus_vampyrus:/132908:/; s/rhinolophus_ferrumequinum:/59479:/; s/myotis_lucifugus:/59463:/; s/vicugna_pacos:/30538:/; s/camelus_dromedarius:/9838:/; s/balaenoptera_musculus:/9771:/; s/physeter_catodon:/9755:/; s/tursiops_truncatus:/9739:/; s/phocoena_sinus:/42100:/; s/delphinapterus_leucas:/9749:/; s/monodon_monoceros:/40151:/; s/cervus_hanglu_yarkandensis:/84702:/; s/moschus_moschiferus:/68415:/; s/capra_hircus:/9925:/; s/ovis_aries_rambouillet:/9940:/; s/bison_bison_bison:/43346:/; s/bos_indicus_hybrid:/30522:/; s/bos_taurus:/9913:/; s/bos_mutus:/72004:/; s/bos_grunniens:/30521:/; s/catagonus_wagneri:/51154:/; s/sus_scrofa:/9823:/; s/dasypus_novemcinctus:/9361:/; s/choloepus_hoffmanni:/9358:/; s/echinops_telfairi:/9371:/; s/loxodonta_africana:/9785:/; s/procavia_capensis:/9813:/; s/monodelphis_domestica:/13616:/; s/sarcophilus_harrisii:/9305:/; s/notamacropus_eugenii:/9315:/; s/phascolarctos_cinereus:/38626:/; s/vombatus_ursinus:/29139:/; s/ornithorhynchus_anatinus:/9258:/; s/sphenodon_punctatus:/8508:/; s/podarcis_muralis:/64176:/; s/salvator_merianae:/96440:/; s/anolis_carolinensis:/28377:/; s/naja_naja:/35670:/; s/laticauda_laticaudata:/8630:/; s/notechis_scutatus:/8663:/; s/pseudonaja_textilis:/8673:/; s/pelodiscus_sinensis:/13735:/; s/gopherus_evgoodei:/1825980:/; s/chelonoidis_abingdonii:/106734:/; s/terrapene_carolina_triunguis:/2587831:/; s/chrysemys_picta_bellii:/8478:/; s/crocodylus_porosus:/8502:/; s/struthio_camelus_australis:/441894:/; s/anser_brachyrhynchus:/132585:/; s/anas_platyrhynchos_platyrhynchos:/8840:/; s/coturnix_japonica:/93934:/; s/gallus_gallus:/9031:/; s/meleagris_gallopavo:/9103:/; s/strigops_habroptila:/2489341:/; s/aquila_chrysaetos_chrysaetos:/223781:/; s/ficedula_albicollis:/59894:/; s/parus_major:/9157:/; s/geospiza_fortis:/48883:/; s/serinus_canaria:/9135:/; s/taeniopygia_guttata:/59729:/; s/xenopus_tropicalis:/8364:/; s/leptobrachium_leishanense:/445787:/; s/latimeria_chalumnae:/7897:/; s/erpetoichthys_calabaricus:/27687:/; s/lepisosteus_oculatus:/7918:/; s/paramormyrops_kingsleyae:/1676925:/; s/scleropages_formosus:/113540:/; s/denticeps_clupeoides:/299321:/; s/clupea_harengus:/7950:/; s/ictalurus_punctatus:/7998:/; s/electrophorus_electricus:/8005:/; s/astyanax_mexicanus:/7994:/; s/pygocentrus_nattereri:/42514:/; s/danio_rerio:/7955:/; s/carassius_auratus:/7957:/; s/sinocyclocheilus_grahami:/75366:/; s/cyprinus_carpio_carpio:/630221:/; s/esox_lucius:/8010:/; s/hucho_hucho:/62062:/; s/salmo_salar:/8030:/; s/salmo_trutta:/8032:/; s/oncorhynchus_mykiss:/8022:/; s/oncorhynchus_kisutch:/8019:/; s/oncorhynchus_tshawytscha:/74940:/; s/gadus_morhua:/8049:/; s/myripristis_murdjan:/586833:/; s/hippocampus_comes:/109280:/; s/tetraodon_nigroviridis:/99883:/; s/takifugu_rubripes:/31033:/; s/labrus_bergylta:/56723:/; s/sparus_aurata:/8175:/; s/larimichthys_crocea:/215358:/; s/dicentrarchus_labrax:/13489:/; s/sander_lucioperca:/283035:/; s/cottoperca_gobio:/56716:/; s/gasterosteus_aculeatus:/69293:/; s/cyclopterus_lumpus:/8103:/; s/mastacembelus_armatus:/205130:/; s/betta_splendens:/158456:/; s/anabas_testudineus:/64144:/; s/cynoglossus_semilaevis:/244447:/; s/scophthalmus_maximus:/52904:/; s/lates_calcarifer:/8187:/; s/seriola_lalandi_dorsalis:/1841481:/; s/seriola_dumerili:/41447:/; s/oryzias_melastigma:/30732:/; s/oryzias_javanicus:/123683:/; s/oryzias_sinensis:/183150:/; s/oryzias_latipes:/8090:/; s/nothobranchius_furzeri:/105023:/; s/kryptolebias_marmoratus:/37003:/; s/fundulus_heteroclitus:/8078:/; s/cyprinodon_variegatus:/28743:/; s/xiphophorus_maculatus:/8083:/; s/poecilia_reticulata:/8081:/; s/poecilia_latipinna:/48699:/; s/poecilia_formosa:/48698:/; s/stegastes_partitus:/144197:/; s/acanthochromis_polyacanthus:/80966:/; s/amphiprion_ocellaris:/80972:/; s/amphiprion_percula:/161767:/; s/amphilophus_citrinellus:/61819:/; s/oreochromis_niloticus:/8128:/; s/neolamprologus_brichardi:/32507:/; s/haplochromis_burtoni:/8153:/; s/pundamilia_nyererei:/303518:/; s/astatotilapia_calliptera:/8154:/; s/maylandia_zebra:/106582:/; s/callorhinchus_milii:/7868:/; s/eptatretus_burgeri:/7764:/; s/petromyzon_marinus:/7757:/; s/ciona_intestinalis:/7719:/; s/ciona_savignyi:/51511:/; s/drosophila_melanogaster:/7227:/; s/caenorhabditis_elegans:/6239:/; s/saccharomyces_cerevisiae:/559292:/;" > species_tree.branch_len.reordered_in_dendroscope.taxids_backup_before_removing_extra_species.nh
# ...then remove extra mus_musculus and sus_scrofa species in Dendroscope (anything that is still text, not numeric)
# This led to 200 species here:
# $speciestreefile = "../compara/input/species_tree.branch_len.reordered_in_dendroscope.taxids.nw";
# $speciestreefile = "../compara/input/species_tree.branch_len.reordered_in_dendroscope.yeastrooted.taxids.internal_nodes_labelled.nw";
$speciestreefile = "../compara/input/species_tree_final_taxids.nw";


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

startme("Running TreeBeST for '$unispec' ($species) using MAFFT alignments of type '$mafftmode'", 0, Numrows($mainquery));
starttime();

$free = 0;
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
		addme("couldn't find FASTA file for ensp (skipped)", $ensp);
		next;
	}
	
	# Check MAFFT alignment file
	$mafftfile = "output/mafft.$ensp.$mafftmode.txt";
	if (!-s $mafftfile)
	{
		addme("couldn't find MAFFT file for ensp (skipped)", $ensp);
		next;
	}

	# Check output CDS alignment file
	$cdsfile = "output/cdstax.$ensp.$mafftmode.txt";
	if (!-s $cdsfile)
	{
		addme("couldn't find back-translated CDS file for ensp (skipped)", $ensp);
		next;
	}
	
	# Wait for cluster nodes to free up
	while ($free <= 0)
	{
		$free = freenodes();
	}

	# Run TreeBeST

	# Way to run Treebest to use method 2, same as Compara
	# -D is for debug, so it shows what it's doing on STDERR
	# treebest best -D -f in.spectree.nh -o out.tree.nhx in.align.mfa
	
	$treeoutfile = "output/tree.$ensp.$mafftmode.txt";

    if (-s $treeoutfile)
    {
	    # Skip this ensp if output file already exists and is non-zero
        addme("output tree already exists for ensp (skipped)", $ensp);
        next;
		# # Crash instead
		# die("Error: Output file '$treeoutfile' already exists for ensp '$ensp' (should run clean.pl first)");
    }

	# The input species tree should be in the New Hampshire format and can be multifurcated. Each internal node MUST have a taxon name. Here is a simple example:
	# 
	# (HUMAN*,(RAT*,MOUSE*)Murinae,(CANFA*-dog,PIG-comment)Laurasiatheria)Eutheria
	# In this example, HUMAN, RAT, MOUSE, CANFA and PIG are species names. Murinae, Laurasiatheria and Eutheria are taxon names. A star `*' indicates that the species is completely sequenced and any gene losses should be counted. A hyphen `-' marks the start of a comment which will not be parsed. Both star and hyphen will not be parsed as species name.
	
	chdir("tmp"); # make sure log files end up in ./tmp/
	run("TreeBeST $mafftmode", "~/scripts/qsub.sh treebest best -D -f ../$speciestreefile -o ../$treeoutfile ../$cdsfile", 1);
	chdir("..");
	
	$free--;
	
	stepme(1);
}
stopme();
stoptime();

# showmesome(0);
showmeall();

done();
