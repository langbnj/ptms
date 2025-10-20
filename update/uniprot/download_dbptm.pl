#!/usr/bin/env perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');

# Text-mining-based PTM database
# There is a new release, now (2022):
# "The download of experimental PTM sites in dbPTM"
# "Due to the inaccessibility of database contents in several online PTM resources, a total 41 biological databases related to PTMs are integrated in dbPTM. To solve the heterogeneity among the data collected from different sources, the reported modification sites are mapped to the UniProtKB protein entries using sequence comparison. With the high-throughput of mass spectrometry-based methods in post-translational proteomics, this update also includes manually curated MS/MS-identified peptides associated with PTMs from research articles through a literature survey. First, a table list of PTM-related keywords is constructed by referring to the UniProtKB PTM list (http://www.uniprot.org/docs/ptmlist.txt) and the annotations of RESID. Then, all fields in the PubMed database are searched based on the keywords of the constructed table list. This is then followed by downloading the full text of the research articles. For the various experiments of proteomic identification, a text-mining system is developed to survey full-text literature that potentially describes the site-specific identification of modified sites.Furthermore, in order to determine the locations of PTMs on a full-length protein sequence, the experimentally verified MS/MS peptides are then mapped to UniProtKB protein entries based on its database identifier (ID) and sequence identity. In the process of data mapping, MS/MS peptides that cannot align exactly to a protein sequence are discarded. Finally, each mapped PTM site is attributed with a corresponding literature (PubMed ID). All types of PTM were categorized by the modified amino acid, with tab-delimited format. These datasets provide UniProt ID, modified position, PTM type, and the sequence with upstream 10 amino acids to downstream 10 amino acids. However, some types of PTM, which were occurred in N-terminal or C-terminal protein, were extracted the sequences with dashes ('-')."
# NAR https://academic.oup.com/nar/article/50/D1/D471/6426061 (PMID 34788852)
# https://awi.cuhk.edu.cn/dbPTM/


# Download

# # dbPTM version 3, last updated: 2012-09-24 (apparently, from the file timestamp)
# run("Download", "wget -nv -O input/dbPTM.tgz 'http://dbptm.mbc.nctu.edu.tw/download/dbPTM.tgz'");
# run("Unpack", "tar -xzvpf dbPTM.tgz");
# run("Clean up", "rm -f dbPTM.tgz");

run("Make dbPTM download directory", "mkdir -p input/dbptm");

run("Download Download page", "wget -nv -O input/unimod_dbptm_download.php 'https://awi.cuhk.edu.cn/dbPTM/download.php'");
$down = `cat input/unimod_dbptm_download.php`;
startme("Downloading individual PTM archives");
cd("input/dbptm");
while ($down =~ /download\/experiment\/([^"]+).gz/g)
{
    stepme(1);

    $ptm = $1;
    # URL-encode PTM name
    $tmpptm = $ptm;
    $tmpptm =~ s/ /%20/g;
    # # Escape spaces with underscores in PTM name for file name
    # $ptm =~ s/ /_/g;

    # print " >> $ptm >> $tmpptm\n";
    # print " >> $ptm >> $tmpptm >> wget -c -O 'input/dbptm/$ptm.gz' 'https://awi.cuhk.edu.cn/dbPTM/download/experiment/$tmpptm.gz'\n";
    # print " >> input/dbptm/$ptm.gz >> https://awi.cuhk.edu.cn/dbPTM/download/experiment/$tmpptm.gz\n";

    run("Download", "wget -nv -O '$ptm.tar.gz' 'https://awi.cuhk.edu.cn/dbPTM/download/experiment/$tmpptm.gz'");
    run("Unpack", "tar -xzvpf '$ptm.tar.gz'");
    run("Clean up", "rm -f '$ptm.tar.gz'");
}
cd("../..");
stopme();

# show directory
run("ls", "ls -1 input/dbptm");

# Can safely concatenate all these into one file
run("Concatenate dbPTM individual PTM files into a single file", "cat input/dbptm/* > input/unimod_dbptm.txt");
run("Clean up", "rm -f input/unimod_dbptm_download.php");
# run("Clean up", "rm -rf input/dbptm");

# show directory
run("ls", "ls -1 input");



done();
