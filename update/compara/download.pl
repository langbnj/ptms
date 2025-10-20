#!/usr/bin/env perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');




# Download ENSEMBL
# $rel = 90;
# $rel = 104; 	# May 2021
# $rel = 105;				# December 2021 (used in UniProt 2022_01)
# $rel = 106;				# April 2022 (presumably will be used in the upcoming UniProt 2022_02)
# $rel = 107;
$rel = 108;

# This explains why Compara doesn't contain all Ensembl species anymore (it only contains 200, out of 311):
# (Performance reasons)
# https://www.ensembl.info/2021/01/20/important-changes-of-data-availability-in-ensembl-gene-trees-and-biomart/

# run("Download", "wget -O input/Compara.protein_trees.aa.fasta.gz 'ftp://ftp.ensembl.org/pub/release-$rel/emf/ensembl-compara/homologies/Compara.$rel.protein_default.aa.fasta.gz'");
# run("Download", "wget -O input/Compara.protein_trees.cds.fasta.gz 'ftp://ftp.ensembl.org/pub/release-$rel/emf/ensembl-compara/homologies/Compara.$rel.protein_default.cds.fasta.gz'");
# run("Download", "wget -O input/Compara.protein_trees.nh.emf.gz 'ftp://ftp.ensembl.org/pub/release-$rel/emf/ensembl-compara/homologies/Compara.$rel.protein_default.nh.emf.gz'");
# run("Download", "wget -O input/Compara.protein_trees.nhx.emf.gz 'ftp://ftp.ensembl.org/pub/release-$rel/emf/ensembl-compara/homologies/Compara.$rel.protein_default.nhx.emf.gz'");

run("Download", "wget -N -P input 'ftp://ftp.ensembl.org/pub/release-$rel/emf/ensembl-compara/homologies/Compara.$rel.protein_default.aa.fasta.gz'");
run("Download", "wget -N -P input 'ftp://ftp.ensembl.org/pub/release-$rel/emf/ensembl-compara/homologies/Compara.$rel.protein_default.cds.fasta.gz'");
run("Download", "wget -N -P input 'ftp://ftp.ensembl.org/pub/release-$rel/emf/ensembl-compara/homologies/Compara.$rel.protein_default.nh.emf.gz'");
run("Download", "wget -N -P input 'ftp://ftp.ensembl.org/pub/release-$rel/emf/ensembl-compara/homologies/Compara.$rel.protein_default.nhx.emf.gz'");



# Species trees

# run("Download latest complete species topology (this is the only complete one) (no branch lengths in here though)", "wget -O input/species_tree.eukaryotes.topology.nw 'http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/ensembl-compara/scripts/pipeline/species_tree.eukaryotes.topology.nw?root=ensembl&view=co'");
# run("Download latest complete species topology (this is the only complete one) (no branch lengths in here though)", "wget -O input/species_tree.ensembl.topology.nw 'https://raw.githubusercontent.com/Ensembl/ensembl-compara/release/$rel/scripts/pipeline/species_tree.ensembl.topology.nw'");



# "A tree with branch-lengths computed in-house "

# Species tree from Ensembl GitHub, with branch lengths:
# https://github.com/Ensembl/ensembl-compara/blob/release/106/conf/vertebrates/species_tree.branch_len.nw
# Contains 227 species rather than the expected 200 in Compara, but that should be okay. It seems to be things like multiple mouse strains that add to the total.

# List of 227 "allowed species":
# https://github.com/Ensembl/ensembl-compara/blob/release/106/conf/vertebrates/allowed_species.json

# However, Compara then actually only uses 200!

# This tree is officially linked to at https://useast.ensembl.org/info/genome/compara/species_trees.html, so it should really be fine to use.
# "A tree with branch-lengths computed in-house"
# This page also says that the protein trees pipeline actually uses the topology of the NCBI taxonomy, i.e. no branch lengths.

# This has the best-looking branch lengths (vs. cafe tree below).
run("Download latest species tree from GitHub (has 227 species, and nice-looking floating point branch lengths, and this is the only one that has correct species names)", "wget -O input/species_tree.branch_len.nw 'https://raw.githubusercontent.com/Ensembl/ensembl-compara/release/$rel/conf/vertebrates/species_tree.branch_len.nw'");

# Reordered so human is at the top:
# species_tree.branch_len.reordered_in_dendroscope.nw
# Reordering steps:
# - Click "Homo sapiens"
# - Click the "reroot" button in the top toolbar
# - Click the "ladderize right" button in the top toolbar
# - Click File > Export to Newick format (saving as .nw rather than .tree)
# Then, it gets a lot more tedious: MAFFT needs a version of this tree that has numeric taxon IDs instead of species names.
# species_tree.branch_len.reordered_in_dendroscope.taxids.nw
# - Use query "SELECT CONCAT('s/',species,':/',tax,'*:/; ') FROM comparaspecies;" to get regexp replacements and run them on species_tree.branch_len.reordered_in_dendroscope.nw (the asterisks are for TreeBeST to mark them as "fully assembled" species, according to https://treesoft.sourceforge.net/treebest.shtml).
# - Delete any still-text (unreplaced) species (which will be extra mus_musculus and sus_scrofa breeds) in Dendroscope by shift-selecting their labels and pressing Cmd-Delete.
# # - In a text editor, manually add 199 labels at every bifurcation (after every ")"). I chose two-character "aa", "ab" etc. labels.
# - Run:
#   label_internal_nodes_in_tree.pl input/species_tree.branch_len.reordered_in_dendroscope.yeastrooted.taxids.nw
#   to produce 
#   input/species_tree.branch_len.reordered_in_dendroscope.yeastrooted.taxids.internal_nodes_labelled.nw
# - Done! TreeBeST can run now.




# "the CAFE pipelines (Gene Gain/Loss trees), with branch lengths coming from the TimeTree database"

# Another option could be: default_protein-trees_cafe.nh
# This one has exactly 200 species (so none would need to be removed), and it has integer branch lengths that look fairly good and comparable to species_tree.branch_len.nw.
# Some branch lengths are very short, though, such as the Oryzias species.
# compara_branchlengths.R definitely shows that the cafe tree has more outliers and a more "powerful" power law distribution. species_tree.branch_len.nw would be preferable, if I can find a good way to remove the extra species (or if TreeBeST can tolerate them).

# According to https://useast.ensembl.org/info/genome/compara/species_trees.html, this cafe tree is for "the CAFE pipelines (Gene Gain/Loss trees), with branch lengths coming from the TimeTree database".
# This could actually be quite good ("TimeTree"). However, molecular evolutionary time probably makes more sense than million years.
# Might be worth trying both.
run("Download latest protein-trees CAFE species tree from Ensembl FTP (has 200 species, and integer branch lengths from TimeTree, but the species names don't quite match those of Ensembl (see mapping below))", "wget -O input/default_protein-trees_cafe.nh 'ftp://ftp.ensembl.org/pub/release-$rel/compara/species_trees/default_protein-trees_cafe.nh'");

# # vertebrates_species-tree_Ensembl.nh here looked promising, but the species names in it don't match well:
# # Mapping table (made by hand for release 106):
# Ailuropoda_melanoleuca_reference_isolate|ailuropoda_melanoleuca
# Anolis_carolinensis_reference_strain|anolis_carolinensis
# Bos_indicus_x_Bos_taurus|bos_indicus_hybrid
# Caenorhabditis_elegans_strain_N2|caenorhabditis_elegans
# Canis_lupus_familiaris_breed_Labrador_retriever|canis_lupus_familiaris
# Capra_hircus_reference_breed|capra_hircus
# Cebus_imitator|cebus_capucinus
# Cricetulus_griseus|cricetulus_griseus_chok1gshd
# Ficedula_albicollis_reference_strain|ficedula_albicollis
# Gorilla_gorilla_gorilla|gorilla_gorilla
# Heterocephalus_glaber|heterocephalus_glaber_female
# Meleagris_gallopavo_reference_strain|meleagris_gallopavo
# Mus_caroli_strain_CAROLI_EIJ|mus_caroli
# Mus_musculus_reference_CL57BL6_strain|mus_musculus
# Mus_pahari_strain_PAHARI_EIJ|mus_pahari
# Mus_spretus_strain_SPRET/EiJ|mus_spretus
# Ovis_aries_reference_breed|ovis_aries_rambouillet
# Rattus_norvegicus_strain_BN/NHsdMcwi|rattus_norvegicus
# Saccharomyces_cerevisiae_S288C|saccharomyces_cerevisiae
# Scophthalmus_maximus_reference_strain|scophthalmus_maximus





# To download all input/alternative_species_trees:
# wget -N -r -e robots=off -l 1 http://ftp.ensembl.org/pub/release-108/compara/species_trees/
# Some of these contain all the comparahomology clades:
# ~/update/comparanopara/test_alternative_input_format/species_trees >> mv ftp.ensembl.org/pub/release-108/compara/species_trees .
# ~/update/comparanopara/test_alternative_input_format/species_trees >> g Homininae * | g Hominidae | g Hominoidea | g Catarrhini | g Simiiformes | g Haplorrhini | g Primates | g Euarchontoglires | g Boreoeutheria | g Eutheria | g Theria | g Mammalia | g Amniota | g Tetrapoda | g Sarcopterygii | g Chordata | g Euteleostomi | g Gnathostomata | g Vertebrata | g Opisthokonta | g Bilateria
# >> default_ncRNA-trees_full_species_tree.nh
# >> default_protein-trees_cafe.nh
# >> default_protein-trees_default.nh
# >> vertebrates_species-tree_Ensembl.nh
# >> vertebrates_species-tree_NCBI Taxonomy.nh
# ...but none of them contain all 200 Compara species as well:
# SELECT species FROM comparaspecies; # These (200 species) are all written identically to those in ~/update/compara/input/species_tree_final.nh (227 species)
# ~/update/comparanopara/test_alternative_input_format/species_trees >> g Homininae * | g Hominidae | g Hominoidea | g Catarrhini | g Simiiformes | g Haplorrhini | g Primates | g Euarchontoglires | g Boreoeutheria | g Eutheria | g Theria | g Mammalia | g Amniota | g Tetrapoda | g Sarcopterygii | g Chordata | g Euteleostomi | g Gnathostomata | g Vertebrata | g Opisthokonta | g Bilateria| g acanthochromis_polyacanthus | g ailuropoda_melanoleuca | g amphilophus_citrinellus | g amphiprion_ocellaris | g amphiprion_percula | g anabas_testudineus | g anas_platyrhynchos_platyrhynchos | g anolis_carolinensis | g anser_brachyrhynchus | g aotus_nancymaae | g aquila_chrysaetos_chrysaetos | g astatotilapia_calliptera | g astyanax_mexicanus | g balaenoptera_musculus | g betta_splendens | g bison_bison_bison | g bos_grunniens | g bos_indicus_hybrid | g bos_mutus | g bos_taurus | g caenorhabditis_elegans | g callithrix_jacchus | g callorhinchus_milii | g camelus_dromedarius | g canis_lupus_dingo | g canis_lupus_familiaris | g capra_hircus | g carassius_auratus | g carlito_syrichta | g catagonus_wagneri | g cavia_porcellus | g cebus_capucinus | g cercocebus_atys | g cervus_hanglu_yarkandensis | g chelonoidis_abingdonii | g chinchilla_lanigera | g chlorocebus_sabaeus | g choloepus_hoffmanni | g chrysemys_picta_bellii | g ciona_intestinalis | g ciona_savignyi | g clupea_harengus | g cottoperca_gobio | g coturnix_japonica | g cricetulus_griseus_chok1gshd | g crocodylus_porosus | g cyclopterus_lumpus | g cynoglossus_semilaevis | g cyprinodon_variegatus | g cyprinus_carpio_carpio | g danio_rerio | g dasypus_novemcinctus | g delphinapterus_leucas | g denticeps_clupeoides | g dicentrarchus_labrax | g dipodomys_ordii | g drosophila_melanogaster | g echinops_telfairi | g electrophorus_electricus | g eptatretus_burgeri | g equus_asinus_asinus | g equus_caballus | g erinaceus_europaeus | g erpetoichthys_calabaricus | g esox_lucius | g felis_catus | g ficedula_albicollis | g fundulus_heteroclitus | g gadus_morhua | g gallus_gallus | g gasterosteus_aculeatus | g geospiza_fortis | g gopherus_evgoodei | g gorilla_gorilla | g haplochromis_burtoni | g heterocephalus_glaber_female | g hippocampus_comes | g homo_sapiens | g hucho_hucho | g ictalurus_punctatus | g ictidomys_tridecemlineatus | g jaculus_jaculus | g kryptolebias_marmoratus | g labrus_bergylta | g larimichthys_crocea | g lates_calcarifer | g laticauda_laticaudata | g latimeria_chalumnae | g lepisosteus_oculatus | g leptobrachium_leishanense | g loxodonta_africana | g macaca_fascicularis | g macaca_mulatta | g macaca_nemestrina | g mandrillus_leucophaeus | g marmota_marmota_marmota | g mastacembelus_armatus | g maylandia_zebra | g meleagris_gallopavo | g mesocricetus_auratus | g microcebus_murinus | g microtus_ochrogaster | g monodelphis_domestica | g monodon_monoceros | g moschus_moschiferus | g mustela_putorius_furo | g mus_caroli | g mus_musculus | g mus_pahari | g mus_spicilegus | g mus_spretus | g myotis_lucifugus | g myripristis_murdjan | g naja_naja | g nannospalax_galili | g neolamprologus_brichardi | g neovison_vison | g nomascus_leucogenys | g notamacropus_eugenii | g notechis_scutatus | g nothobranchius_furzeri | g ochotona_princeps | g octodon_degus | g oncorhynchus_kisutch | g oncorhynchus_mykiss | g oncorhynchus_tshawytscha | g oreochromis_niloticus | g ornithorhynchus_anatinus | g oryctolagus_cuniculus | g oryzias_javanicus | g oryzias_latipes | g oryzias_melastigma | g oryzias_sinensis | g otolemur_garnettii | g ovis_aries_rambouillet | g panthera_leo | g panthera_pardus | g panthera_tigris_altaica | g pan_paniscus | g pan_troglodytes | g papio_anubis | g paramormyrops_kingsleyae | g parus_major | g pelodiscus_sinensis | g peromyscus_maniculatus_bairdii | g petromyzon_marinus | g phascolarctos_cinereus | g phocoena_sinus | g physeter_catodon | g podarcis_muralis | g poecilia_formosa | g poecilia_latipinna | g poecilia_reticulata | g pongo_abelii | g procavia_capensis | g prolemur_simus | g propithecus_coquereli | g pseudonaja_textilis | g pteropus_vampyrus | g pundamilia_nyererei | g pygocentrus_nattereri | g rattus_norvegicus | g rhinolophus_ferrumequinum | g rhinopithecus_bieti | g rhinopithecus_roxellana | g saccharomyces_cerevisiae | g saimiri_boliviensis_boliviensis | g salmo_salar | g salmo_trutta | g salvator_merianae | g sander_lucioperca | g sarcophilus_harrisii | g sciurus_vulgaris | g scleropages_formosus | g scophthalmus_maximus | g serinus_canaria | g seriola_dumerili | g seriola_lalandi_dorsalis | g sinocyclocheilus_grahami | g sorex_araneus | g sparus_aurata | g sphenodon_punctatus | g stegastes_partitus | g strigops_habroptila | g struthio_camelus_australis | g sus_scrofa | g taeniopygia_guttata | g takifugu_rubripes | g terrapene_carolina_triunguis | g tetraodon_nigroviridis | g tupaia_belangeri | g tursiops_truncatus | g urocitellus_parryii | g ursus_americanus | g ursus_maritimus | g vicugna_pacos | g vombatus_ursinus | g vulpes_vulpes | g xenopus_tropicalis | g xiphophorus_maculatus
# >> None
# >> Can't use any of these alternative species trees to get clades.



# Vertebrates tree (not described at https://useast.ensembl.org/info/genome/compara/species_trees.html):
# The branch lengths and topology in this one and in the GitHub species_tree.branch_len.nw agree very well, as far as I can see in Dendroscope. Using species_tree.branch_len.nw since its names are correct
# # run("Download latest species tree (despite the name, this contains all Compara species, not just vertebrates) (has branch lengths)", "wget -O input/vertebrates_species-tree_Ensembl.nh 'http://ftp.ensembl.org/pub/release-$rel/compara/species_trees/vertebrates_species-tree_Ensembl.nh'");
# run("Download latest species tree with branch lengths for 39 mammals (this one has accurate branch lengths according to the Ensembl Help Desk, but is mammals only)", "wget -O input/species_tree.39mammals.branch_len.nw 'http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/ensembl-compara/scripts/pipeline/species_tree.39mammals.branch_len.nw?root=ensembl&view=co'");
 





# show directory
run("ls", "ls -lah input");
state("Downloaded ENSEMBL Compara release '$rel'");

done();
