#!/usr/bin/env perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');


# start

# Cite http://nar.oxfordjournals.org/content/43/D1/D512.full.pdf+html?sid=e6f09a39-3ea8-4f63-a63b-469430ae096a

# # # # # # # Note - PhosphoSitePlus is continuously updated (current version 110313, from first line of files, i.e. 2013-11-03)
# # # # # # Note - PhosphoSitePlus is continuously updated (current version 121417, from first line of files, i.e. 2017-12-14) (and 2017-12-20 from file timestamps)
# # # # # Note - PhosphoSitePlus is continuously updated (current version 022618, from first line of files, i.e. 2018-02-26) (and the same from the file timestamps)
# # # # Note - PhosphoSitePlus is continuously updated (current version 040122, from first line of files, i.e. 2022-04-01) (and the same from the .gz file timestamps, and also the PhosphoSitePlus download section)
# # # Note - PhosphoSitePlus is continuously updated (current version 060322, from first line of files, i.e. 2022-06-03) (and the same from the .gz file timestamps, and also the PhosphoSitePlus download section)
# # Note - PhosphoSitePlus is continuously updated (current version 080822, from first line of files, i.e. 2022-08-08) (and the same from the .gz file timestamps, and also the PhosphoSitePlus download section)
# Note - PhosphoSitePlus is continuously updated (current version 102022, from first line of files, i.e. 2022-10-20) (and the same from the .gz file timestamps, and also the PhosphoSitePlus download section: Thu Oct 20 15:20:26 EDT 2022)

warn("
- Download manually from https://www.phosphosite.org/staticDownloads (anything ending in _site_dataset (7), and Kinase_Substrate_Dataset.gz and Regulatory_sites.gz, although these last two aren't used anywhere yet) (9 files in total)
- mkdir -p input/phosphositeplus
- Copy into input/phosphositeplus
- Run
");

# - Download manually from https://www.phosphosite.org/staticDownloads (anything ending in _site_dataset, and Kinase_Substrate_Dataset and Regulatory_sites, although these last two aren't used anywhere yet) (9 files total)
# - mkdir input/phosphositeplus
# - Copy into input/phosphositeplus
# - Run:
run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/Acetylation_site_dataset.gz' 'input/unimod_phosphositeplus_K-ac.txt.gz'");
run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/Methylation_site_dataset.gz' 'input/unimod_phosphositeplus_KR-me.txt.gz'");
run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/O-GalNAc_site_dataset.gz' 'input/unimod_phosphositeplus_O-gal.txt.gz'");
run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/O-GlcNAc_site_dataset.gz' 'input/unimod_phosphositeplus_O-g.txt.gz'");
run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/Phosphorylation_site_dataset.gz' 'input/unimod_phosphositeplus_STY-p.txt.gz'");
run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/Sumoylation_site_dataset.gz' 'input/unimod_phosphositeplus_K-sumo.txt.gz'");
run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/Ubiquitination_site_dataset.gz' 'input/unimod_phosphositeplus_K-ub.txt.gz'");

run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/Kinase_Substrate_Dataset.gz' 'input/unimod_kinase_targets_phosphositeplus.txt.gz'");
run("Rename PhosphoSitePlus manual download", "mv -fv 'input/phosphositeplus/Regulatory_sites.gz' 'input/unimod_regulatory_sites_phosphositeplus.txt.gz'");

# Remove input/phosphositeplus
run("Remove temporary directory", "rmdir -v input/phosphositeplus");

# Unpack
run("Unpack", "gunzip -f input/unimod_phosphositeplus_K-ac.txt.gz");
run("Unpack", "gunzip -f input/unimod_phosphositeplus_KR-me.txt.gz");
run("Unpack", "gunzip -f input/unimod_phosphositeplus_O-gal.txt.gz");
run("Unpack", "gunzip -f input/unimod_phosphositeplus_O-g.txt.gz");
run("Unpack", "gunzip -f input/unimod_phosphositeplus_STY-p.txt.gz");
run("Unpack", "gunzip -f input/unimod_phosphositeplus_K-sumo.txt.gz");
run("Unpack", "gunzip -f input/unimod_phosphositeplus_K-ub.txt.gz");

run("Unpack", "gunzip -f input/unimod_kinase_targets_phosphositeplus.txt.gz");
run("Unpack", "gunzip -f input/unimod_regulatory_sites_phosphositeplus.txt.gz");

# show directory
run("ls", "ls -1 input");

done();
