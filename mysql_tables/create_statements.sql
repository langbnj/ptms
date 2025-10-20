CREATE TABLE `PTMs` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(45) NOT NULL DEFAULT '',
  `display` varchar(45) NOT NULL DEFAULT '',
  `class` varchar(45) NOT NULL DEFAULT '',
  `amino` varchar(15) NOT NULL DEFAULT '',
  PRIMARY KEY (`id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='List of PTMs of interest';

CREATE TABLE `comparacds` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(25) DEFAULT NULL,
  `species` varchar(50) DEFAULT NULL,
  `aln` int DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='ENSEMBL Compara GeneTrees CDS Alignments (FASTA format)';

CREATE TABLE `comparaenspspecies` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) DEFAULT NULL,
  `species` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='ENSEMBL Compara ENSP to Species mapping';

CREATE TABLE `comparafasta` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='ENSEMBL Compara GeneTrees Alignments (FASTA format)';

CREATE TABLE `comparafasta_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Alignments (FASTA format)';

CREATE TABLE `comparafasta_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Alignments (FASTA format)';

CREATE TABLE `comparafasta_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Alignments (FASTA format)';

CREATE TABLE `comparafasta_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Alignments (FASTA format)';

CREATE TABLE `comparafasta_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Alignments (FASTA format)';

CREATE TABLE `comparafasta_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Alignments (FASTA format)';

CREATE TABLE `comparafasta_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Alignments (FASTA format)';

CREATE TABLE `comparahomology` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg1` varchar(25) DEFAULT NULL,
  `enst1` varchar(25) DEFAULT NULL,
  `ensp1` varchar(25) DEFAULT NULL,
  `acc1` varchar(25) DEFAULT NULL,
  `name1` varchar(25) DEFAULT NULL,
  `symbol1` varchar(50) DEFAULT NULL,
  `species1` varchar(25) DEFAULT NULL,
  `fullspecies1` varchar(50) DEFAULT NULL,
  `ensg2` varchar(25) DEFAULT NULL,
  `enst2` varchar(25) DEFAULT NULL,
  `ensp2` varchar(25) DEFAULT NULL,
  `acc2` varchar(25) DEFAULT NULL,
  `name2` varchar(25) DEFAULT NULL,
  `symbol2` varchar(50) DEFAULT NULL,
  `species2` varchar(25) DEFAULT NULL,
  `fullspecies2` varchar(50) DEFAULT NULL,
  `homology` varchar(25) DEFAULT NULL,
  `hc` tinyint DEFAULT NULL,
  `clade` varchar(25) DEFAULT NULL,
  `target_percent_id` float DEFAULT NULL,
  `query_percent_id` float DEFAULT NULL,
  `gene_order_score` tinyint DEFAULT NULL,
  `whole_genome_align_score` float DEFAULT NULL,
  `target_chr` varchar(50) DEFAULT NULL,
  `target_chr_start` int DEFAULT NULL,
  `target_chr_stop` int DEFAULT NULL,
  `target_chr_strand` char(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ensg1` (`ensg1`),
  KEY `Enst1` (`enst1`),
  KEY `Ensp1` (`ensp1`),
  KEY `Acc1` (`acc1`),
  KEY `Name1` (`name1`),
  KEY `Symbol1` (`symbol1`),
  KEY `Species1` (`species1`),
  KEY `Fullspecies1` (`fullspecies1`),
  KEY `Ensg2` (`ensg2`),
  KEY `Enst2` (`enst2`),
  KEY `Ensp2` (`ensp2`),
  KEY `Acc2` (`acc2`),
  KEY `Name2` (`name2`),
  KEY `Symbol2` (`symbol2`),
  KEY `Species2` (`species2`),
  KEY `Fullspecies2` (`fullspecies2`),
  KEY `Homology` (`homology`),
  KEY `Hc` (`hc`),
  KEY `Clade` (`clade`),
  KEY `Target_percent_id` (`target_percent_id`),
  KEY `Query_percent_id` (`query_percent_id`),
  KEY `Gene_order_score` (`gene_order_score`),
  KEY `Whole_genome_align_score` (`whole_genome_align_score`),
  KEY `Target_chr` (`target_chr`),
  KEY `Target_chr_start` (`target_chr_start`),
  KEY `Target_chr_stop` (`target_chr_stop`),
  KEY `Target_chr_strand` (`target_chr_strand`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='ENSEMBL Compara homology relationships';

CREATE TABLE `comparanh` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='ENSEMBL Compara GeneTrees (Newick format, NH)';

CREATE TABLE `comparanh_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Trees (Newick format, NH)';

CREATE TABLE `comparanh_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Trees (Newick format, NH)';

CREATE TABLE `comparanh_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Trees (Newick format, NH)';

CREATE TABLE `comparanh_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Trees (Newick format, NH)';

CREATE TABLE `comparanh_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Trees (Newick format, NH)';

CREATE TABLE `comparanh_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Trees (Newick format, NH)';

CREATE TABLE `comparanh_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Trees (Newick format, NH)';

CREATE TABLE `comparanhx` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` mediumtext NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='ENSEMBL Compara GeneTrees (Newick format, NHX)';

CREATE TABLE `comparaspecies` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `species` varchar(50) DEFAULT NULL,
  `display` varchar(50) DEFAULT NULL,
  `shortdisplay` varchar(50) DEFAULT NULL,
  `unispec` varchar(5) DEFAULT NULL,
  `tax` mediumint DEFAULT NULL,
  `type` varchar(25) DEFAULT NULL,
  `homininae` tinyint DEFAULT NULL,
  `greatape` tinyint DEFAULT NULL,
  `primate` tinyint DEFAULT NULL,
  `mammal` tinyint DEFAULT NULL,
  `vertebrate` tinyint DEFAULT NULL,
  `metazoan` tinyint DEFAULT NULL,
  `eukaryote` tinyint DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Species` (`species`),
  KEY `Unispec` (`unispec`),
  KEY `Type` (`type`),
  KEY `Primate` (`primate`),
  KEY `Mammal` (`mammal`),
  KEY `Vertebrate` (`vertebrate`),
  KEY `Metazoan` (`metazoan`),
  KEY `Eukaryote` (`eukaryote`),
  KEY `Tax` (`tax`),
  KEY `Greatape` (`greatape`),
  KEY `Homininae` (`homininae`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Compara species table (in Compara species tree order)';

CREATE TABLE `ensembl` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `species` varchar(15) DEFAULT NULL,
  `fullspecies` varchar(45) DEFAULT NULL,
  `ensgv` varchar(25) DEFAULT NULL,
  `ensg` varchar(25) DEFAULT NULL,
  `enstv` varchar(25) DEFAULT NULL,
  `enst` varchar(25) DEFAULT NULL,
  `enspv` varchar(25) DEFAULT NULL,
  `ensp` varchar(25) DEFAULT NULL,
  `seq` text,
  `cds` mediumtext,
  `symbol` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Species` (`species`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Ensp` (`ensp`),
  KEY `Ensgv` (`ensgv`),
  KEY `Enstv` (`enstv`),
  KEY `Enspv` (`enspv`),
  KEY `Symbol` (`symbol`),
  KEY `Fullspecies` (`fullspecies`),
  KEY `Seq` (`seq`(200)),
  KEY `Cds` (`cds`(200))
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='ENSEMBL proteomes';

CREATE TABLE `ensembl_gff3_cds` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `enst` varchar(25) DEFAULT NULL,
  `ensp` varchar(25) DEFAULT NULL,
  `fullspecies` varchar(50) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `phase` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Enst` (`enst`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`),
  KEY `Phase` (`phase`),
  KEY `Chr` (`chr`),
  KEY `Fullspecies` (`fullspecies`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Ensembl CDS annotation from GFF3 files';

CREATE TABLE `ensembl_gff3_exon` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `enst` varchar(25) DEFAULT NULL,
  `biotype` varchar(50) DEFAULT NULL,
  `ense` varchar(25) DEFAULT NULL,
  `ensev` varchar(25) DEFAULT NULL,
  `fullspecies` varchar(50) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `constitutive` tinyint(1) DEFAULT NULL,
  `rank` mediumint DEFAULT NULL,
  `ensembl_phase` tinyint(1) DEFAULT NULL,
  `ensembl_end_phase` tinyint(1) DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `coding` tinyint(1) DEFAULT NULL,
  `cds_start` int unsigned DEFAULT NULL,
  `cds_stop` int unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `phase` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ense` (`ense`),
  KEY `Ensev` (`ensev`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Constitutive` (`constitutive`),
  KEY `Rank` (`rank`),
  KEY `Ensembl_phase` (`ensembl_phase`),
  KEY `Ensembl_end_phase` (`ensembl_end_phase`),
  KEY `Enst` (`enst`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`),
  KEY `Phase` (`phase`),
  KEY `Chr` (`chr`),
  KEY `Cds_start` (`cds_start`),
  KEY `Cds_stop` (`cds_stop`),
  KEY `Biotype` (`biotype`),
  KEY `Coding` (`coding`),
  KEY `Fullspecies` (`fullspecies`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Ensembl exon annotation from GFF3 files';

CREATE TABLE `ensembl_gff3_gene` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(50) DEFAULT NULL,
  `ensgv` varchar(50) DEFAULT NULL,
  `fullspecies` varchar(50) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `logic_name` varchar(250) DEFAULT NULL,
  `symbol` varchar(250) DEFAULT NULL,
  `description` varchar(2500) DEFAULT NULL,
  `biotype` varchar(50) DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ensg` (`ensg`),
  KEY `Ensgv` (`ensgv`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Logic_name` (`logic_name`),
  KEY `Symbol` (`symbol`),
  KEY `Biotype` (`biotype`),
  KEY `Chr` (`chr`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`),
  KEY `Fullspecies` (`fullspecies`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Ensembl gene annotation from GFF3 files';

CREATE TABLE `ensembl_gff3_transcript` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(50) DEFAULT NULL,
  `enst` varchar(50) DEFAULT NULL,
  `enstv` varchar(50) DEFAULT NULL,
  `fullspecies` varchar(50) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `symbol` varchar(250) DEFAULT NULL,
  `symbolv` varchar(250) DEFAULT NULL,
  `biotype` varchar(50) DEFAULT NULL,
  `tsl` tinyint(1) DEFAULT NULL,
  `basic` tinyint(1) DEFAULT NULL,
  `appris` varchar(4) DEFAULT NULL,
  `longest` tinyint(1) DEFAULT NULL,
  `longest_for_symbol` tinyint(1) DEFAULT NULL,
  `ccds` varchar(25) DEFAULT NULL,
  `ccdsv` varchar(25) DEFAULT NULL,
  `tags` varchar(1200) DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `seq` mediumtext,
  `aaseq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Enstv` (`enstv`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Symbol` (`symbol`),
  KEY `Symbolv` (`symbolv`),
  KEY `Biotype` (`biotype`),
  KEY `Tsl` (`tsl`),
  KEY `Basic` (`basic`),
  KEY `Ccds` (`ccds`),
  KEY `Ccdsv` (`ccdsv`),
  KEY `Chr` (`chr`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`),
  KEY `Appris` (`appris`),
  KEY `Seq` (`seq`(8)),
  KEY `Aaseq` (`aaseq`(8)),
  KEY `Longest` (`longest`),
  KEY `Longest_for_symbol` (`longest_for_symbol`),
  KEY `Fullspecies` (`fullspecies`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Ensembl transcript annotation from GFF3 files';

CREATE TABLE `ensembl_gtf` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(25) DEFAULT NULL,
  `enst` varchar(25) DEFAULT NULL,
  `ensp` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `symbol` varchar(250) DEFAULT NULL,
  `biotype` varchar(25) DEFAULT NULL,
  `tsl` tinyint(1) DEFAULT NULL,
  `basic` tinyint(1) DEFAULT NULL,
  `ccds` varchar(25) DEFAULT NULL,
  `ensembl_tags` varchar(1200) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Symbol` (`symbol`),
  KEY `Biotype` (`biotype`),
  KEY `Tsl` (`tsl`),
  KEY `Basic` (`basic`),
  KEY `Ccds` (`ccds`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Ensembl transcript annotation from GTF files';

CREATE TABLE `ensembl_species` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(250) DEFAULT NULL,
  `fullspecies` varchar(250) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `tax` int DEFAULT NULL,
  `division` varchar(250) DEFAULT NULL,
  `assembly_id` varchar(250) DEFAULT NULL,
  `assembly_name` varchar(250) DEFAULT NULL,
  `genebuild` varchar(250) DEFAULT NULL,
  `nproteincoding` int DEFAULT NULL,
  `nproteincodinguniprotkbswissprot` int DEFAULT NULL,
  `nproteincodinguniprotkbtrembl` int DEFAULT NULL,
  `uniprotcoverage` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Assembly_id` (`assembly_id`),
  KEY `Assembly_name` (`assembly_name`),
  KEY `Division` (`division`),
  KEY `Genebuild` (`genebuild`),
  KEY `Name` (`name`),
  KEY `Nproteincoding` (`nproteincoding`),
  KEY `Nproteincodinguniprotkbswissprot` (`nproteincodinguniprotkbswissprot`),
  KEY `Nproteincodinguniprotkbtrembl` (`nproteincodinguniprotkbtrembl`),
  KEY `Species` (`fullspecies`),
  KEY `Taxonomy_id` (`tax`),
  KEY `Unispec` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='input-uniprot_report_EnsemblVertebrates-108.txt';

CREATE TABLE `evorate_analysis` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `species` varchar(250) DEFAULT NULL,
  `predictor` varchar(250) DEFAULT NULL,
  `evorate` varchar(250) DEFAULT NULL,
  `mafftmode` varchar(250) DEFAULT NULL,
  `source` varchar(250) DEFAULT NULL,
  `mode` varchar(250) DEFAULT NULL,
  `resamples` int DEFAULT NULL,
  `ptm` varchar(250) DEFAULT NULL,
  `dis` varchar(250) DEFAULT NULL,
  `nptm` int DEFAULT NULL,
  `ncontrol` int DEFAULT NULL,
  `medianptm` double DEFAULT NULL,
  `mediancontrol` double DEFAULT NULL,
  `meanptm` double DEFAULT NULL,
  `meancontrol` double DEFAULT NULL,
  `pmw_nonresampled` double DEFAULT NULL,
  `pmw_nonresampled_averaged_median` double DEFAULT NULL,
  `pmw_nonresampled_averaged_mean` double DEFAULT NULL,
  `pmw_nonresampled_averaged_median_unpaired` double DEFAULT NULL,
  `pmw_nonresampled_averaged_mean_unpaired` double DEFAULT NULL,
  `ptt_nonresampled` double DEFAULT NULL,
  `pmw` double DEFAULT NULL,
  `ptt` double DEFAULT NULL,
  `pmedian` double DEFAULT NULL,
  `pmean` double DEFAULT NULL,
  `mediandif` double DEFAULT NULL,
  `meandif` double DEFAULT NULL,
  `wilcox_successes` int DEFAULT NULL,
  `ttest_successes` int DEFAULT NULL,
  `median_successes` int DEFAULT NULL,
  `mean_successes` int DEFAULT NULL,
  `pmw_nonresampled_corr` double DEFAULT NULL,
  `pmw_nonresampled_averaged_median_corr` double DEFAULT NULL,
  `pmw_nonresampled_averaged_mean_corr` double DEFAULT NULL,
  `pmw_nonresampled_averaged_median_unpaired_corr` double DEFAULT NULL,
  `pmw_nonresampled_averaged_mean_unpaired_corr` double DEFAULT NULL,
  `ptt_nonresampled_corr` double DEFAULT NULL,
  `pmwcorr` double DEFAULT NULL,
  `pttcorr` double DEFAULT NULL,
  `pmediancorr` double DEFAULT NULL,
  `pmeancorr` double DEFAULT NULL,
  `pmw_nonresampled_sig` double DEFAULT NULL,
  `pmw_nonresampled_averaged_median_sig` double DEFAULT NULL,
  `pmw_nonresampled_averaged_mean_sig` double DEFAULT NULL,
  `pmw_nonresampled_averaged_median_unpaired_sig` double DEFAULT NULL,
  `pmw_nonresampled_averaged_mean_unpaired_sig` double DEFAULT NULL,
  `ptt_nonresampled_sig` double DEFAULT NULL,
  `pmwsig` double DEFAULT NULL,
  `pttsig` double DEFAULT NULL,
  `pmediansig` double DEFAULT NULL,
  `pmeansig` double DEFAULT NULL,
  `sig` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Dis` (`dis`),
  KEY `Evorate` (`evorate`),
  KEY `Mafftmode` (`mafftmode`),
  KEY `Mean_successes` (`mean_successes`),
  KEY `Median_successes` (`median_successes`),
  KEY `Mode` (`mode`),
  KEY `Ncontrol` (`ncontrol`),
  KEY `Nptm` (`nptm`),
  KEY `Predictor` (`predictor`),
  KEY `Ptm` (`ptm`),
  KEY `Resamples` (`resamples`),
  KEY `Source` (`source`),
  KEY `Species` (`species`),
  KEY `Ttest_successes` (`ttest_successes`),
  KEY `Wilcox_successes` (`wilcox_successes`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='output-evorate_analysis-10000-human.tsv';

CREATE TABLE `evorate_capra0_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra0_einsi';

CREATE TABLE `evorate_capra0_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra0_einsi_tree';

CREATE TABLE `evorate_capra0_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra0_einsi_tree_1para';

CREATE TABLE `evorate_capra0_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra0_ginsi';

CREATE TABLE `evorate_capra0_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra0_ginsi_tree';

CREATE TABLE `evorate_capra0_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra0_linsi';

CREATE TABLE `evorate_capra0_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra0_linsi_tree';

CREATE TABLE `evorate_capra0_para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra0_para';

CREATE TABLE `evorate_capra1_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra1_einsi';

CREATE TABLE `evorate_capra1_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra1_einsi_tree';

CREATE TABLE `evorate_capra1_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra1_einsi_tree_1para';

CREATE TABLE `evorate_capra1_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra1_ginsi';

CREATE TABLE `evorate_capra1_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra1_ginsi_tree';

CREATE TABLE `evorate_capra1_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra1_linsi';

CREATE TABLE `evorate_capra1_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra1_linsi_tree';

CREATE TABLE `evorate_capra1_para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra1_para';

CREATE TABLE `evorate_capra3_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra3_einsi';

CREATE TABLE `evorate_capra3_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra3_einsi_tree';

CREATE TABLE `evorate_capra3_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra3_einsi_tree_1para';

CREATE TABLE `evorate_capra3_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra3_ginsi';

CREATE TABLE `evorate_capra3_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra3_ginsi_tree';

CREATE TABLE `evorate_capra3_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra3_linsi';

CREATE TABLE `evorate_capra3_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra3_linsi_tree';

CREATE TABLE `evorate_capra3_para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra3_para';

CREATE TABLE `evorate_capra5_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra5_einsi';

CREATE TABLE `evorate_capra5_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra5_einsi_tree';

CREATE TABLE `evorate_capra5_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra5_einsi_tree_1para';

CREATE TABLE `evorate_capra5_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra5_ginsi';

CREATE TABLE `evorate_capra5_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra5_ginsi_tree';

CREATE TABLE `evorate_capra5_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra5_linsi';

CREATE TABLE `evorate_capra5_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra5_linsi_tree';

CREATE TABLE `evorate_capra5_para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates capra5_para';

CREATE TABLE `evorate_lichtarge_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates lichtarge_einsi';

CREATE TABLE `evorate_lichtarge_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates lichtarge_einsi_tree';

CREATE TABLE `evorate_lichtarge_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates lichtarge_einsi_tree_1para';

CREATE TABLE `evorate_lichtarge_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates lichtarge_ginsi';

CREATE TABLE `evorate_lichtarge_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates lichtarge_ginsi_tree';

CREATE TABLE `evorate_lichtarge_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates lichtarge_linsi';

CREATE TABLE `evorate_lichtarge_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates lichtarge_linsi_tree';

CREATE TABLE `evorate_lichtarge_para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`name`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates lichtarge_para';

CREATE TABLE `evorate_norm_capra0_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra0_einsi';

CREATE TABLE `evorate_norm_capra0_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra0_einsi_tree';

CREATE TABLE `evorate_norm_capra0_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra0_einsi_tree_1para';

CREATE TABLE `evorate_norm_capra0_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra0_ginsi';

CREATE TABLE `evorate_norm_capra0_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra0_ginsi_tree';

CREATE TABLE `evorate_norm_capra0_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra0_linsi';

CREATE TABLE `evorate_norm_capra0_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra0_linsi_tree';

CREATE TABLE `evorate_norm_capra1_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra1_einsi';

CREATE TABLE `evorate_norm_capra1_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra1_einsi_tree';

CREATE TABLE `evorate_norm_capra1_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra1_einsi_tree_1para';

CREATE TABLE `evorate_norm_capra1_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra1_ginsi';

CREATE TABLE `evorate_norm_capra1_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra1_ginsi_tree';

CREATE TABLE `evorate_norm_capra1_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra1_linsi';

CREATE TABLE `evorate_norm_capra1_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra1_linsi_tree';

CREATE TABLE `evorate_norm_capra3_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra3_einsi';

CREATE TABLE `evorate_norm_capra3_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra3_einsi_tree';

CREATE TABLE `evorate_norm_capra3_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra3_einsi_tree_1para';

CREATE TABLE `evorate_norm_capra3_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra3_ginsi';

CREATE TABLE `evorate_norm_capra3_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra3_ginsi_tree';

CREATE TABLE `evorate_norm_capra3_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra3_linsi';

CREATE TABLE `evorate_norm_capra3_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra3_linsi_tree';

CREATE TABLE `evorate_norm_capra5_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra5_einsi';

CREATE TABLE `evorate_norm_capra5_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra5_einsi_tree';

CREATE TABLE `evorate_norm_capra5_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra5_einsi_tree_1para';

CREATE TABLE `evorate_norm_capra5_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra5_ginsi';

CREATE TABLE `evorate_norm_capra5_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra5_ginsi_tree';

CREATE TABLE `evorate_norm_capra5_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra5_linsi';

CREATE TABLE `evorate_norm_capra5_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_capra5_linsi_tree';

CREATE TABLE `evorate_norm_lichtarge_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_lichtarge_einsi';

CREATE TABLE `evorate_norm_lichtarge_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_lichtarge_einsi_tree';

CREATE TABLE `evorate_norm_lichtarge_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_lichtarge_einsi_tree_1para';

CREATE TABLE `evorate_norm_lichtarge_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_lichtarge_ginsi';

CREATE TABLE `evorate_norm_lichtarge_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_lichtarge_ginsi_tree';

CREATE TABLE `evorate_norm_lichtarge_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_lichtarge_linsi';

CREATE TABLE `evorate_norm_lichtarge_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_lichtarge_linsi_tree';

CREATE TABLE `evorate_norm_rate4site_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_rate4site_einsi';

CREATE TABLE `evorate_norm_rate4site_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_rate4site_einsi_tree';

CREATE TABLE `evorate_norm_rate4site_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_rate4site_einsi_tree_1para';

CREATE TABLE `evorate_norm_rate4site_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_rate4site_ginsi';

CREATE TABLE `evorate_norm_rate4site_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_rate4site_ginsi_tree';

CREATE TABLE `evorate_norm_rate4site_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_rate4site_linsi';

CREATE TABLE `evorate_norm_rate4site_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates norm_rate4site_linsi_tree';

CREATE TABLE `evorate_quantile_logos` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `species` varchar(250) DEFAULT NULL,
  `source` varchar(250) DEFAULT NULL,
  `evorate` varchar(250) DEFAULT NULL,
  `predictor` varchar(250) DEFAULT NULL,
  `bootstraps` varchar(250) DEFAULT NULL,
  `dis` varchar(250) DEFAULT NULL,
  `alpha` double DEFAULT NULL,
  `ptm` varchar(250) DEFAULT NULL,
  `quantiles` int DEFAULT NULL,
  `quant` int DEFAULT NULL,
  `position` int DEFAULT NULL,
  `aa` varchar(250) DEFAULT NULL,
  `freq` double DEFAULT NULL,
  `pvalue` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Aa` (`aa`),
  KEY `Bootstraps` (`bootstraps`),
  KEY `Dis` (`dis`),
  KEY `Evorate` (`evorate`),
  KEY `Position` (`position`),
  KEY `Predictor` (`predictor`),
  KEY `Ptm` (`ptm`),
  KEY `Quant` (`quant`),
  KEY `Quantiles` (`quantiles`),
  KEY `Source` (`source`),
  KEY `Species` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='output-evorate_quantile_logos-human.tsv';

CREATE TABLE `evorate_rate4site_einsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates rate4site_einsi';

CREATE TABLE `evorate_rate4site_einsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates rate4site_einsi_tree';

CREATE TABLE `evorate_rate4site_einsi_tree_1para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates rate4site_einsi_tree_1para';

CREATE TABLE `evorate_rate4site_ginsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates rate4site_ginsi';

CREATE TABLE `evorate_rate4site_ginsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates rate4site_ginsi_tree';

CREATE TABLE `evorate_rate4site_linsi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates rate4site_linsi';

CREATE TABLE `evorate_rate4site_linsi_tree` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates rate4site_linsi_tree';

CREATE TABLE `evorate_rate4site_para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` char(15) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `rate` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `NameSite` (`name`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates rate4site_para';

CREATE TABLE `evorate_strat` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `nosurf` varchar(250) DEFAULT NULL,
  `test` varchar(250) DEFAULT NULL,
  `alternative` varchar(250) DEFAULT NULL,
  `source` varchar(250) DEFAULT NULL,
  `evorate` varchar(250) DEFAULT NULL,
  `mafftmode` varchar(250) DEFAULT NULL,
  `pred` varchar(250) DEFAULT NULL,
  `minsize` int DEFAULT NULL,
  `disfilt` varchar(250) DEFAULT NULL,
  `strat` varchar(250) DEFAULT NULL,
  `pval` double DEFAULT NULL,
  `sig` varchar(250) DEFAULT NULL,
  `min` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Alternative` (`alternative`),
  KEY `Disfilt` (`disfilt`),
  KEY `Evorate` (`evorate`),
  KEY `Mafftmode` (`mafftmode`),
  KEY `Min` (`min`),
  KEY `Minsize` (`minsize`),
  KEY `Nosurf` (`nosurf`),
  KEY `Pred` (`pred`),
  KEY `Sig` (`sig`),
  KEY `Source` (`source`),
  KEY `Strat` (`strat`),
  KEY `Test` (`test`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='output-table.tsv';

CREATE TABLE `gencode39` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) DEFAULT NULL,
  `enst` varchar(250) DEFAULT NULL,
  `ensg` varchar(250) DEFAULT NULL,
  `otthumg` varchar(250) DEFAULT NULL,
  `otthumt` varchar(250) DEFAULT NULL,
  `orf` varchar(250) DEFAULT NULL,
  `symbol` varchar(250) DEFAULT NULL,
  `seqlen` int DEFAULT NULL,
  `seq` text,
  PRIMARY KEY (`id`),
  KEY `Ensg` (`ensg`),
  KEY `Ensp` (`ensp`),
  KEY `Enst` (`enst`),
  KEY `Orf` (`orf`),
  KEY `Otthumg` (`otthumg`),
  KEY `Otthumt` (`otthumt`),
  KEY `Seqlen` (`seqlen`),
  KEY `Symbol` (`symbol`),
  KEY `Seq` (`seq`(25))
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='input/translations_39.txt';

CREATE TABLE `gencode_gff3_cds` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(25) DEFAULT NULL,
  `ensgv` varchar(25) DEFAULT NULL,
  `enst` varchar(25) DEFAULT NULL,
  `enstv` varchar(25) DEFAULT NULL,
  `ensp` varchar(25) DEFAULT NULL,
  `enspv` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `phase` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Enst` (`enst`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`),
  KEY `Phase` (`phase`),
  KEY `Chr` (`chr`),
  KEY `Enstv` (`enstv`),
  KEY `Enspv` (`enspv`),
  KEY `Ensg` (`ensg`),
  KEY `Ensgv` (`ensgv`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='GENCODE CDS annotation from GFF3 files';

CREATE TABLE `gencode_gff3_exon` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(25) DEFAULT NULL,
  `ensgv` varchar(25) DEFAULT NULL,
  `enst` varchar(25) DEFAULT NULL,
  `enstv` varchar(25) DEFAULT NULL,
  `biotype` varchar(50) DEFAULT NULL,
  `ense` varchar(25) DEFAULT NULL,
  `ensev` varchar(25) DEFAULT NULL,
  `exonid` varchar(35) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `exon_number` mediumint DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `coding` tinyint(1) DEFAULT NULL,
  `cds_start` int unsigned DEFAULT NULL,
  `cds_stop` int unsigned DEFAULT NULL,
  `phase` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ense` (`ense`),
  KEY `Ensev` (`ensev`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Exon_number` (`exon_number`),
  KEY `Enst` (`enst`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`),
  KEY `Phase` (`phase`),
  KEY `Chr` (`chr`),
  KEY `Cds_start` (`cds_start`),
  KEY `Cds_stop` (`cds_stop`),
  KEY `Biotype` (`biotype`),
  KEY `Coding` (`coding`),
  KEY `Enstv` (`enstv`),
  KEY `Exonid` (`exonid`),
  KEY `Ensg` (`ensg`),
  KEY `Ensgv` (`ensgv`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='GENCODE exon annotation from GFF3 files';

CREATE TABLE `gencode_gff3_gene` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(25) DEFAULT NULL,
  `ensgv` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `symbol` varchar(250) DEFAULT NULL,
  `biotype` varchar(50) DEFAULT NULL,
  `level` tinyint(1) DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ensg` (`ensg`),
  KEY `Ensgv` (`ensgv`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Symbol` (`symbol`),
  KEY `Biotype` (`biotype`),
  KEY `Level` (`level`),
  KEY `Chr` (`chr`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='GENCODE gene annotation from GFF3 files';

CREATE TABLE `gencode_gff3_transcript` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(25) DEFAULT NULL,
  `ensgv` varchar(25) DEFAULT NULL,
  `enst` varchar(25) DEFAULT NULL,
  `enstv` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `symbol` varchar(250) DEFAULT NULL,
  `symbolv` varchar(250) DEFAULT NULL,
  `biotype` varchar(50) DEFAULT NULL,
  `tsl` tinyint(1) DEFAULT NULL,
  `basic` tinyint(1) DEFAULT NULL,
  `appris` varchar(4) DEFAULT NULL,
  `longest` tinyint(1) DEFAULT NULL,
  `longest_for_symbol` tinyint(1) DEFAULT NULL,
  `ncbiids` varchar(250) DEFAULT NULL,
  `ccds` varchar(25) DEFAULT NULL,
  `ccdsv` varchar(25) DEFAULT NULL,
  `tags` varchar(1200) DEFAULT NULL,
  `chr` varchar(50) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `seq` mediumtext,
  `aaseq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Enstv` (`enstv`),
  KEY `Species` (`species`),
  KEY `Source` (`source`),
  KEY `Symbol` (`symbol`),
  KEY `Symbolv` (`symbolv`),
  KEY `Biotype` (`biotype`),
  KEY `Tsl` (`tsl`),
  KEY `Basic` (`basic`),
  KEY `Ccds` (`ccds`),
  KEY `Ccdsv` (`ccdsv`),
  KEY `Ensgv` (`ensgv`),
  KEY `Chr` (`chr`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Strand` (`strand`),
  KEY `Seq` (`seq`(8)),
  KEY `Aaseq` (`aaseq`(8)),
  KEY `Ncbiid` (`ncbiids`),
  KEY `Appris` (`appris`),
  KEY `Longest_for_symbol` (`longest_for_symbol`),
  KEY `Longest` (`longest`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='GENCODE transcript annotation from GFF3 files';

CREATE TABLE `gencode_uniprot` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(25) DEFAULT NULL,
  `ensgv` varchar(25) DEFAULT NULL,
  `symbol` varchar(25) DEFAULT NULL,
  `symbolv` varchar(25) DEFAULT NULL,
  `enst` varchar(25) DEFAULT NULL,
  `enstv` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `primary_acc` varchar(25) DEFAULT NULL,
  `seq` mediumtext,
  `aaseq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensg` (`ensg`),
  KEY `Ensgv` (`ensgv`),
  KEY `Enst` (`enst`),
  KEY `Enstv` (`enstv`),
  KEY `Species` (`species`),
  KEY `Acc` (`acc`),
  KEY `Symbol` (`symbol`),
  KEY `Symbolv` (`symbolv`),
  KEY `Primary_acc` (`primary_acc`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='GENCODE UniProt mapping';

CREATE TABLE `hgnc` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `hgnc_id` varchar(250) DEFAULT NULL,
  `symbol` varchar(250) DEFAULT NULL,
  `name` varchar(250) DEFAULT NULL,
  `locus_group` varchar(250) DEFAULT NULL,
  `locus_type` varchar(250) DEFAULT NULL,
  `status` varchar(250) DEFAULT NULL,
  `location` varchar(250) DEFAULT NULL,
  `location_sortable` varchar(250) DEFAULT NULL,
  `alias_symbol` varchar(250) DEFAULT NULL,
  `alias_name` text,
  `prev_symbol` varchar(250) DEFAULT NULL,
  `prev_name` text,
  `gene_group` varchar(250) DEFAULT NULL,
  `gene_group_id` varchar(250) DEFAULT NULL,
  `date_approved_reserved` varchar(250) DEFAULT NULL,
  `date_symbol_changed` varchar(250) DEFAULT NULL,
  `date_name_changed` varchar(250) DEFAULT NULL,
  `date_modified` varchar(250) DEFAULT NULL,
  `entrez_id` int DEFAULT NULL,
  `ensembl_gene_id` varchar(250) DEFAULT NULL,
  `vega_id` varchar(250) DEFAULT NULL,
  `ucsc_id` varchar(250) DEFAULT NULL,
  `ena` varchar(250) DEFAULT NULL,
  `refseq_accession` varchar(250) DEFAULT NULL,
  `ccds_id` text,
  `uniprot_ids` varchar(250) DEFAULT NULL,
  `pubmed_id` varchar(250) DEFAULT NULL,
  `mgd_id` varchar(250) DEFAULT NULL,
  `rgd_id` varchar(250) DEFAULT NULL,
  `lsdb` text,
  `cosmic` varchar(250) DEFAULT NULL,
  `omim_id` varchar(250) DEFAULT NULL,
  `mirbase` varchar(250) DEFAULT NULL,
  `homeodb` int DEFAULT NULL,
  `snornabase` varchar(250) DEFAULT NULL,
  `bioparadigms_slc` varchar(250) DEFAULT NULL,
  `orphanet` int DEFAULT NULL,
  `pseudogene_org` varchar(250) DEFAULT NULL,
  `horde_id` varchar(250) DEFAULT NULL,
  `merops` varchar(250) DEFAULT NULL,
  `imgt` varchar(250) DEFAULT NULL,
  `iuphar` varchar(250) DEFAULT NULL,
  `kznf_gene_catalog` int DEFAULT NULL,
  `mamit_trnadb` int DEFAULT NULL,
  `cd` varchar(250) DEFAULT NULL,
  `lncrnadb` varchar(250) DEFAULT NULL,
  `enzyme_id` varchar(250) DEFAULT NULL,
  `intermediate_filament_db` int DEFAULT NULL,
  `rna_central_ids` int DEFAULT NULL,
  `lncipedia` varchar(250) DEFAULT NULL,
  `gtrnadb` varchar(250) DEFAULT NULL,
  `agr` varchar(250) DEFAULT NULL,
  `mane_select` varchar(250) DEFAULT NULL,
  `gencc` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Agr` (`agr`),
  KEY `Alias_symbol` (`alias_symbol`),
  KEY `Bioparadigms_slc` (`bioparadigms_slc`),
  KEY `Cd` (`cd`),
  KEY `Cosmic` (`cosmic`),
  KEY `Date_approved_reserved` (`date_approved_reserved`),
  KEY `Date_modified` (`date_modified`),
  KEY `Date_name_changed` (`date_name_changed`),
  KEY `Date_symbol_changed` (`date_symbol_changed`),
  KEY `Ena` (`ena`),
  KEY `Ensembl_gene_id` (`ensembl_gene_id`),
  KEY `Entrez_id` (`entrez_id`),
  KEY `Enzyme_id` (`enzyme_id`),
  KEY `Gencc` (`gencc`),
  KEY `Gene_group` (`gene_group`),
  KEY `Gene_group_id` (`gene_group_id`),
  KEY `Gtrnadb` (`gtrnadb`),
  KEY `Hgnc_id` (`hgnc_id`),
  KEY `Homeodb` (`homeodb`),
  KEY `Horde_id` (`horde_id`),
  KEY `Imgt` (`imgt`),
  KEY `Intermediate_filament_db` (`intermediate_filament_db`),
  KEY `Iuphar` (`iuphar`),
  KEY `Kznf_gene_catalog` (`kznf_gene_catalog`),
  KEY `Lncipedia` (`lncipedia`),
  KEY `Lncrnadb` (`lncrnadb`),
  KEY `Location` (`location`),
  KEY `Location_sortable` (`location_sortable`),
  KEY `Locus_group` (`locus_group`),
  KEY `Locus_type` (`locus_type`),
  KEY `Mamit_trnadb` (`mamit_trnadb`),
  KEY `Mane_select` (`mane_select`),
  KEY `Merops` (`merops`),
  KEY `Mgd_id` (`mgd_id`),
  KEY `Mirbase` (`mirbase`),
  KEY `Name` (`name`),
  KEY `Omim_id` (`omim_id`),
  KEY `Orphanet` (`orphanet`),
  KEY `Prev_symbol` (`prev_symbol`),
  KEY `Pseudogene_org` (`pseudogene_org`),
  KEY `Pubmed_id` (`pubmed_id`),
  KEY `Refseq_accession` (`refseq_accession`),
  KEY `Rgd_id` (`rgd_id`),
  KEY `Rna_central_ids` (`rna_central_ids`),
  KEY `Snornabase` (`snornabase`),
  KEY `Status` (`status`),
  KEY `Symbol` (`symbol`),
  KEY `Ucsc_id` (`ucsc_id`),
  KEY `Uniprot_ids` (`uniprot_ids`),
  KEY `Vega_id` (`vega_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='input/hgnc_complete_set_noquotes.txt';

CREATE TABLE `hgnc_aliases` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `hgncid` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `symbol` varchar(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `alias` varchar(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `ensg` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `accs` varchar(250) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Hgncid` (`hgncid`),
  KEY `Symbol` (`symbol`),
  KEY `Alias` (`alias`),
  KEY `Ensg` (`ensg`),
  KEY `Accs` (`accs`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci COMMENT='HGNC alternative gene symbols utility table';

CREATE TABLE `lists` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `list` varchar(255) NOT NULL DEFAULT '',
  `type` varchar(255) NOT NULL DEFAULT '',
  `category` varchar(255) NOT NULL DEFAULT '',
  `item` varchar(255) NOT NULL DEFAULT '',
  `map` varchar(255) NOT NULL DEFAULT '',
  `maptype` varchar(255) NOT NULL DEFAULT '',
  PRIMARY KEY (`id`),
  KEY `List` (`list`),
  KEY `Type` (`type`),
  KEY `Category` (`category`),
  KEY `Name` (`item`),
  KEY `Item` (`item`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='General purpose table for lists of e.g. proteins';

CREATE TABLE `m_csa` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(250) DEFAULT NULL,
  `acc` varchar(250) DEFAULT NULL,
  `species` varchar(250) DEFAULT NULL,
  `site` int DEFAULT NULL,
  `aa` varchar(250) DEFAULT NULL,
  `aa3` varchar(250) DEFAULT NULL,
  `ptm` varchar(250) DEFAULT NULL,
  `mcsa_id` int DEFAULT NULL,
  `pdb_id` varchar(250) DEFAULT NULL,
  `site_pdb` int DEFAULT NULL,
  `chain_name` varchar(250) DEFAULT NULL,
  `assembly_chain_name` varchar(250) DEFAULT NULL,
  `assembly` int DEFAULT NULL,
  `aa_pdb` varchar(250) DEFAULT NULL,
  `aa3_pdb` varchar(250) DEFAULT NULL,
  `domain_name` varchar(250) DEFAULT NULL,
  `domain_cath_id` varchar(250) DEFAULT NULL,
  `roles_summary` varchar(250) DEFAULT NULL,
  `main_annotation` text,
  PRIMARY KEY (`id`),
  KEY `Aa` (`aa`),
  KEY `Aa3` (`aa3`),
  KEY `Aa3_pdb` (`aa3_pdb`),
  KEY `Aa_pdb` (`aa_pdb`),
  KEY `Acc` (`acc`),
  KEY `Assembly` (`assembly`),
  KEY `Assembly_chain_name` (`assembly_chain_name`),
  KEY `Chain_name` (`chain_name`),
  KEY `Domain_cath_id` (`domain_cath_id`),
  KEY `Domain_name` (`domain_name`),
  KEY `Mcsa_id` (`mcsa_id`),
  KEY `Name` (`name`),
  KEY `Pdb_id` (`pdb_id`),
  KEY `Ptm` (`ptm`),
  KEY `Roles_summary` (`roles_summary`),
  KEY `Site` (`site`),
  KEY `Site_pdb` (`site_pdb`),
  KEY `Species` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='tmp-m_csa.txt';

CREATE TABLE `ncbi` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `species` varchar(15) DEFAULT NULL,
  `gi` varchar(25) DEFAULT NULL,
  `ncbiid` varchar(25) DEFAULT NULL,
  `ncbiidshort` varchar(25) DEFAULT NULL,
  `ncbiname` text,
  `sequence` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Species` (`species`),
  KEY `Gi` (`gi`),
  KEY `Ncbiid` (`ncbiid`),
  KEY `Ncbiidshort` (`ncbiidshort`),
  KEY `Sequence` (`sequence`(8))
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='NCBI RefSeq mRNA + Proteins. Only contains Refseq GI IDs';

CREATE TABLE `ncbi_mane` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `symbol` varchar(25) DEFAULT NULL,
  `species` varchar(15) DEFAULT NULL,
  `ensgv` varchar(25) DEFAULT NULL,
  `ensg` varchar(25) DEFAULT NULL,
  `enstv` varchar(25) DEFAULT NULL,
  `enst` varchar(25) DEFAULT NULL,
  `enspv` varchar(25) DEFAULT NULL,
  `ensp` varchar(25) DEFAULT NULL,
  `nmv` varchar(25) DEFAULT NULL,
  `nm` varchar(25) DEFAULT NULL,
  `npv` varchar(25) DEFAULT NULL,
  `np` varchar(25) DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Symbol` (`symbol`),
  KEY `Species` (`species`),
  KEY `Ensgv` (`ensgv`),
  KEY `Ensg` (`ensg`),
  KEY `Enstv` (`enstv`),
  KEY `Enst` (`enst`),
  KEY `Enspv` (`enspv`),
  KEY `Ensp` (`ensp`),
  KEY `Nmv` (`nmv`),
  KEY `nm` (`nm`),
  KEY `npv` (`npv`),
  KEY `np` (`np`),
  KEY `Seq` (`seq`(50))
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Mapping of NCBI RefSeq IDs to Ensembl IDs (from MANE)';

CREATE TABLE `ptm_contact_aas` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ptm` varchar(250) DEFAULT NULL,
  `acc` varchar(250) DEFAULT NULL,
  `alphacc` varchar(250) DEFAULT NULL,
  `afdb` varchar(250) DEFAULT NULL,
  `site` varchar(250) DEFAULT NULL,
  `source` varchar(250) DEFAULT NULL,
  `class` varchar(250) DEFAULT NULL,
  `type` varchar(250) DEFAULT NULL,
  `dis` varchar(250) DEFAULT NULL,
  `surf` varchar(250) DEFAULT NULL,
  `plddt` varchar(250) DEFAULT NULL,
  `plddt10` varchar(250) DEFAULT NULL,
  `aa1` varchar(250) DEFAULT NULL,
  `aa2` varchar(250) DEFAULT NULL,
  `contact` varchar(250) DEFAULT NULL,
  `atom1` varchar(250) DEFAULT NULL,
  `atom2` varchar(250) DEFAULT NULL,
  `dist` varchar(250) DEFAULT NULL,
  `pae` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Aa1` (`aa1`),
  KEY `Aa2` (`aa2`),
  KEY `Acc` (`acc`),
  KEY `Afdb` (`afdb`),
  KEY `Alphacc` (`alphacc`),
  KEY `Atom1` (`atom1`),
  KEY `Atom2` (`atom2`),
  KEY `Class` (`class`),
  KEY `Contact` (`contact`),
  KEY `Dis` (`dis`),
  KEY `Dist` (`dist`),
  KEY `Pae` (`pae`),
  KEY `Plddt` (`plddt`),
  KEY `Plddt10` (`plddt10`),
  KEY `Ptm` (`ptm`),
  KEY `Site` (`site`),
  KEY `Source` (`source`),
  KEY `Surf` (`surf`),
  KEY `Type` (`type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='output-contacts-human-1000.txt';

CREATE TABLE `snps_clinvar` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(20) DEFAULT NULL,
  `canon` varchar(20) DEFAULT NULL,
  `ensg` varchar(50) DEFAULT NULL,
  `enst` varchar(50) DEFAULT NULL,
  `species` varchar(10) DEFAULT NULL,
  `site` smallint unsigned DEFAULT NULL,
  `original` char(1) DEFAULT NULL,
  `variant` char(1) DEFAULT NULL,
  `clinsig` tinyint(1) DEFAULT NULL,
  `clin` varchar(100) DEFAULT NULL,
  `genetictest` tinyint unsigned DEFAULT NULL,
  `submitters` smallint unsigned DEFAULT NULL,
  `evaluated` date DEFAULT NULL,
  `rsid` varchar(50) DEFAULT NULL,
  `phenoids` text,
  `phenotype` text,
  `origin` varchar(100) DEFAULT NULL,
  `originsimple` varchar(50) DEFAULT NULL,
  `status` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Acc` (`acc`),
  KEY `Site` (`site`),
  KEY `Original` (`original`),
  KEY `Variant` (`variant`),
  KEY `Enst` (`enst`),
  KEY `Clinsig` (`clinsig`),
  KEY `Genetictest` (`genetictest`),
  KEY `Canon` (`canon`),
  KEY `Originsimple` (`originsimple`),
  KEY `Ensg` (`ensg`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='NCBI ClinVar disease-associated variants';

CREATE TABLE `snps_cosmic` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `chr` varchar(5) DEFAULT NULL,
  `pos` int DEFAULT NULL,
  `originalbase` varchar(10) DEFAULT NULL,
  `variantbase` varchar(10) DEFAULT NULL,
  `cosv` varchar(20) DEFAULT NULL,
  `source` varchar(60) DEFAULT NULL,
  `cell_line` tinyint DEFAULT NULL,
  `name` varchar(20) DEFAULT NULL,
  `acc` varchar(20) DEFAULT NULL,
  `canon` varchar(20) DEFAULT NULL,
  `ensg` varchar(20) DEFAULT NULL,
  `enst` varchar(20) DEFAULT NULL,
  `species` varchar(10) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `original` char(1) DEFAULT NULL,
  `variant` char(1) DEFAULT NULL,
  `ac` mediumint DEFAULT NULL,
  `an` mediumint DEFAULT NULL,
  `af` double DEFAULT NULL,
  `cgc` tinyint DEFAULT NULL,
  `loc` varchar(100) DEFAULT NULL,
  `subloc` varchar(100) DEFAULT NULL,
  `hist` varchar(200) DEFAULT NULL,
  `subhist` varchar(200) DEFAULT NULL,
  `pmid` varchar(20) DEFAULT NULL,
  `zyg` varchar(10) DEFAULT NULL,
  `loh` varchar(10) DEFAULT NULL,
  `som` varchar(50) DEFAULT NULL,
  `tum` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Chr` (`chr`),
  KEY `Pos` (`pos`),
  KEY `Originalbase` (`originalbase`),
  KEY `Variantbase` (`variantbase`),
  KEY `Cosv` (`cosv`),
  KEY `Source` (`source`),
  KEY `Cell_line` (`cell_line`),
  KEY `Name` (`name`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Species` (`species`),
  KEY `Site` (`site`),
  KEY `Original` (`original`),
  KEY `Variant` (`variant`),
  KEY `Ac` (`ac`),
  KEY `An` (`an`),
  KEY `Cgc` (`cgc`),
  KEY `Loc` (`loc`),
  KEY `Hist` (`hist`),
  KEY `Pmid` (`pmid`),
  KEY `Zyg` (`zyg`),
  KEY `Loh` (`loh`),
  KEY `Som` (`som`),
  KEY `Tum` (`tum`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Catalogue of Somatic Mutations in Cancer (COSMIC) human cancer missense mutations';

CREATE TABLE `snps_gnomad` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `chr` varchar(5) DEFAULT NULL,
  `pos` int DEFAULT NULL,
  `originalbase` char(1) DEFAULT NULL,
  `variantbase` char(1) DEFAULT NULL,
  `rsid` int DEFAULT NULL,
  `name` varchar(20) DEFAULT NULL,
  `acc` varchar(20) DEFAULT NULL,
  `canon` varchar(20) DEFAULT NULL,
  `ensg` varchar(20) DEFAULT NULL,
  `enst` varchar(20) DEFAULT NULL,
  `species` varchar(10) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `original` char(1) DEFAULT NULL,
  `variant` char(1) DEFAULT NULL,
  `source` varchar(20) DEFAULT NULL,
  `ac` mediumint DEFAULT NULL,
  `ac_xx` mediumint DEFAULT NULL,
  `ac_xy` mediumint DEFAULT NULL,
  `ac_grpmax` mediumint DEFAULT NULL,
  `an` mediumint DEFAULT NULL,
  `an_xx` mediumint DEFAULT NULL,
  `an_xy` mediumint DEFAULT NULL,
  `an_grpmax` mediumint DEFAULT NULL,
  `af` double DEFAULT NULL,
  `af_xx` double DEFAULT NULL,
  `af_xy` double DEFAULT NULL,
  `af_grpmax` double DEFAULT NULL,
  `grpmax` varchar(5) DEFAULT NULL,
  `n_alt_alleles` tinyint DEFAULT NULL,
  `allele_type` varchar(20) DEFAULT NULL,
  `variant_type` varchar(20) DEFAULT NULL,
  `was_mixed` tinyint DEFAULT NULL,
  `nhomalt` mediumint DEFAULT NULL,
  `nhomalt_xx` mediumint DEFAULT NULL,
  `nhomalt_xy` mediumint DEFAULT NULL,
  `nhomalt_grpmax` mediumint DEFAULT NULL,
  `inbreeding_coeff` double DEFAULT NULL,
  `af_controls_and_biobanks` double DEFAULT NULL,
  `af_non_cancer` double DEFAULT NULL,
  `af_non_neuro` double DEFAULT NULL,
  `af_non_topmed` double DEFAULT NULL,
  `af_non_v2` double DEFAULT NULL,
  `af_afr` double DEFAULT NULL,
  `af_ami` double DEFAULT NULL,
  `af_amr` double DEFAULT NULL,
  `af_asj` double DEFAULT NULL,
  `af_eas` double DEFAULT NULL,
  `af_fin` double DEFAULT NULL,
  `af_mid` double DEFAULT NULL,
  `af_nfe` double DEFAULT NULL,
  `af_sas` double DEFAULT NULL,
  `af_remaining` double DEFAULT NULL,
  `fafmax95` double DEFAULT NULL,
  `fafmax_faf95_max` double DEFAULT NULL,
  `faf95` double DEFAULT NULL,
  `faf95_afr` double DEFAULT NULL,
  `faf95_amr` double DEFAULT NULL,
  `faf95_eas` double DEFAULT NULL,
  `faf95_nfe` double DEFAULT NULL,
  `faf95_sas` double DEFAULT NULL,
  `as_fs` double DEFAULT NULL,
  `as_sor` double DEFAULT NULL,
  `as_pab_max` double DEFAULT NULL,
  `fs` double DEFAULT NULL,
  `sor` double DEFAULT NULL,
  `lcr` tinyint DEFAULT NULL,
  `non_par` tinyint DEFAULT NULL,
  `segdup` tinyint DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Species` (`species`),
  KEY `Site` (`site`),
  KEY `Original` (`original`),
  KEY `Variant` (`variant`),
  KEY `Subset` (`source`),
  KEY `Ac` (`ac`),
  KEY `Ac_xx` (`ac_xx`),
  KEY `Ac_xy` (`ac_xy`),
  KEY `Ac_grpmax` (`ac_grpmax`),
  KEY `An` (`an`),
  KEY `An_xx` (`an_xx`),
  KEY `An_xy` (`an_xy`),
  KEY `An_grpmax` (`an_grpmax`),
  KEY `N_alt_alleles` (`n_alt_alleles`),
  KEY `Allele_type` (`allele_type`),
  KEY `Variant_type` (`variant_type`),
  KEY `Was_mixed` (`was_mixed`),
  KEY `Nhomalt` (`nhomalt`),
  KEY `Nhomalt_xx` (`nhomalt_xx`),
  KEY `Nhomalt_xy` (`nhomalt_xy`),
  KEY `Nhomalt_grpmax` (`nhomalt_grpmax`),
  KEY `Lcr` (`lcr`),
  KEY `Non_par` (`non_par`),
  KEY `Segdup` (`segdup`),
  KEY `Chr` (`chr`),
  KEY `Pos` (`pos`),
  KEY `Originalbase` (`originalbase`),
  KEY `Variantbase` (`variantbase`),
  KEY `Rsid` (`rsid`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Natural variants from gnomAD (v2 and v4 exomes, v2, v3 and v4 genomes)';

CREATE TABLE `snps_human` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` char(11) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `original` char(1) DEFAULT NULL,
  `variant` char(1) DEFAULT NULL,
  `chromo` char(2) DEFAULT NULL,
  `dbsnp` char(12) DEFAULT NULL,
  `allelefreq` float DEFAULT NULL,
  `quality` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Site` (`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='1000 Genomes Project SNPs (first wave)';

CREATE TABLE `snps_mastermind` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `chr` varchar(5) DEFAULT NULL,
  `pos` int DEFAULT NULL,
  `originalbase` char(3) DEFAULT NULL,
  `variantbase` char(3) DEFAULT NULL,
  `name` varchar(20) DEFAULT NULL,
  `acc` varchar(20) DEFAULT NULL,
  `canon` varchar(20) DEFAULT NULL,
  `ensg` varchar(20) DEFAULT NULL,
  `enst` varchar(20) DEFAULT NULL,
  `species` varchar(10) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `original` char(1) DEFAULT NULL,
  `variant` char(1) DEFAULT NULL,
  `papers` smallint DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Species` (`species`),
  KEY `Site` (`site`),
  KEY `Original` (`original`),
  KEY `Variant` (`variant`),
  KEY `Papers` (`papers`),
  KEY `Chr` (`chr`),
  KEY `Pos` (`pos`),
  KEY `Originalbase` (`originalbase`),
  KEY `Variantbase` (`variantbase`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Genomenon Mastermind Cited Variants Reference (CVR) variants';

CREATE TABLE `snps_unisnp` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` char(11) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `original` char(1) DEFAULT NULL,
  `variant` char(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Site` (`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt SNPs';

CREATE TABLE `snps_univar` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) DEFAULT NULL,
  `acc` varchar(20) DEFAULT NULL,
  `canon` varchar(20) DEFAULT NULL,
  `ensg` varchar(20) DEFAULT NULL,
  `enst` varchar(20) DEFAULT NULL,
  `species` varchar(10) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `original` char(1) DEFAULT NULL,
  `variant` char(1) DEFAULT NULL,
  `sources` smallint DEFAULT NULL,
  `uniprot` tinyint DEFAULT NULL,
  `cosmic_curated` tinyint DEFAULT NULL,
  `dbsnp` tinyint DEFAULT NULL,
  `tcga` tinyint DEFAULT NULL,
  `thousand_genomes` tinyint DEFAULT NULL,
  `esp` tinyint DEFAULT NULL,
  `clinvar` tinyint DEFAULT NULL,
  `topmed` tinyint DEFAULT NULL,
  `exac` tinyint DEFAULT NULL,
  `gnomad` tinyint DEFAULT NULL,
  `clingen` tinyint DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Species` (`species`),
  KEY `Site` (`site`),
  KEY `Original` (`original`),
  KEY `Variant` (`variant`),
  KEY `Uniprot` (`uniprot`),
  KEY `Cosmic_curated` (`cosmic_curated`),
  KEY `Dbsnp` (`dbsnp`),
  KEY `Tcga` (`tcga`),
  KEY `Thousand_genomes` (`thousand_genomes`),
  KEY `Esp` (`esp`),
  KEY `Clinvar` (`clinvar`),
  KEY `Topmed` (`topmed`),
  KEY `Exac` (`exac`),
  KEY `Gnomad` (`gnomad`),
  KEY `Clingen` (`clingen`),
  KEY `Sources` (`sources`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt variants';

CREATE TABLE `snps_vervet` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensg` varchar(25) DEFAULT NULL,
  `enst` varchar(25) DEFAULT NULL,
  `enstv` varchar(25) DEFAULT NULL,
  `ensp` varchar(25) DEFAULT NULL,
  `enspv` varchar(25) DEFAULT NULL,
  `ensite` int DEFAULT NULL,
  `original` varchar(250) DEFAULT NULL,
  `human` varchar(250) DEFAULT NULL,
  `variant` varchar(250) DEFAULT NULL,
  `name` varchar(25) DEFAULT NULL,
  `site` int DEFAULT NULL,
  `chromo` varchar(250) DEFAULT NULL,
  `allelecount` int DEFAULT NULL,
  `allelefreq` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Ensp` (`ensp`),
  KEY `Name` (`name`),
  KEY `Site` (`site`),
  KEY `Variant` (`variant`),
  KEY `Allelecount` (`allelecount`),
  KEY `Enstv` (`enstv`),
  KEY `Enspv` (`enspv`),
  KEY `Ensite` (`ensite`),
  KEY `Original` (`original`),
  KEY `Human` (`human`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Vervet monkey SNPs (& mapped to human via comparafasta)';

CREATE TABLE `string` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp1` char(18) DEFAULT NULL,
  `ensp2` char(18) DEFAULT NULL,
  `symbol1` char(18) DEFAULT NULL,
  `symbol2` char(18) DEFAULT NULL,
  `species` char(5) DEFAULT NULL,
  `score` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Ensp1` (`ensp1`),
  KEY `Ensp2` (`ensp2`),
  KEY `species` (`species`),
  KEY `Symbol1` (`symbol1`),
  KEY `Symbol2` (`symbol2`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='STRING Protein-protein interactions';

CREATE TABLE `string_genes` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `symbol1` char(18) DEFAULT NULL,
  `symbol2` char(18) DEFAULT NULL,
  `species` char(5) DEFAULT NULL,
  `score` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Symbol1` (`symbol1`),
  KEY `Symbol2` (`symbol2`),
  KEY `species` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='STRING Protein-protein interactions at the gene level';

CREATE TABLE `tax` (
  `id` int NOT NULL AUTO_INCREMENT,
  `tax` mediumint DEFAULT NULL,
  `name` varchar(500) DEFAULT NULL,
  `type` varchar(20) DEFAULT NULL,
  `rank` varchar(20) DEFAULT NULL,
  `parents` varchar(1000) DEFAULT NULL,
  `speciestype` varchar(20) DEFAULT NULL,
  `homininae` tinyint DEFAULT NULL,
  `greatape` tinyint DEFAULT NULL,
  `primate` tinyint DEFAULT NULL,
  `mammal` tinyint DEFAULT NULL,
  `vertebrate` tinyint DEFAULT NULL,
  `metazoan` tinyint DEFAULT NULL,
  `eukaryote` tinyint DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Tax` (`tax`),
  KEY `Speciestype` (`speciestype`),
  KEY `Primate` (`primate`),
  KEY `Mammal` (`mammal`),
  KEY `Vertebrate` (`vertebrate`),
  KEY `Metazoan` (`metazoan`),
  KEY `Eukaryote` (`eukaryote`),
  KEY `Rank` (`rank`),
  KEY `Type` (`type`),
  KEY `Greatape` (`greatape`),
  KEY `Homininae` (`homininae`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='NCBI Taxon ID table';

CREATE TABLE `uniacc` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` char(20) DEFAULT NULL,
  `species` char(15) DEFAULT NULL,
  `acc` char(10) DEFAULT NULL,
  `canon` char(10) DEFAULT NULL,
  `trembl` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Index_1` (`name`),
  KEY `Index_2` (`species`),
  KEY `Index_3` (`acc`),
  KEY `Primary_acc` (`canon`),
  KEY `Trembl` (`trembl`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt accessions';

CREATE TABLE `uniasa` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `site` int DEFAULT NULL,
  `aa` varchar(1) DEFAULT NULL,
  `pdb` varchar(20) DEFAULT NULL,
  `chain` char(1) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `pdbsite` int DEFAULT NULL,
  `method` varchar(25) DEFAULT NULL,
  `resolution` float DEFAULT NULL,
  `asa` float DEFAULT NULL,
  `relasa` float DEFAULT NULL,
  `dis` varchar(25) DEFAULT NULL,
  `surf` varchar(250) DEFAULT NULL,
  `type` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Site` (`site`),
  KEY `Pdb` (`pdb`),
  KEY `Type` (`type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt Accessible Surface Area (new SIFTS version)';

CREATE TABLE `uniasa_blast` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `site` int DEFAULT NULL,
  `aa` varchar(1) DEFAULT NULL,
  `pdb` varchar(20) DEFAULT NULL,
  `chain` char(1) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `pdbsite` int DEFAULT NULL,
  `asa` float DEFAULT NULL,
  `relasa` float DEFAULT NULL,
  `dis` varchar(1) DEFAULT NULL,
  `surf` varchar(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Site` (`site`),
  KEY `Call` (`surf`),
  KEY `Pdb` (`pdb`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt Accessible Surface Area (obsolete BLAST version)';

CREATE TABLE `uniens` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `canonical` tinyint DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `ensgv` varchar(40) DEFAULT NULL,
  `ensg` varchar(40) DEFAULT NULL,
  `enstv` varchar(40) DEFAULT NULL,
  `enst` varchar(40) DEFAULT NULL,
  `enspv` varchar(40) DEFAULT NULL,
  `ensp` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Ensp` (`ensp`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Ensgv` (`ensgv`),
  KEY `Enstv` (`enstv`),
  KEY `Enspv` (`enspv`),
  KEY `Canonical` (`canonical`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt-ENSEMBL mappings (by sequence, mapped using sequence matching)';

CREATE TABLE `uniens_uniprot` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `ensgv` varchar(40) DEFAULT NULL,
  `ensg` varchar(40) DEFAULT NULL,
  `enstv` varchar(40) DEFAULT NULL,
  `enst` varchar(40) DEFAULT NULL,
  `enspv` varchar(40) DEFAULT NULL,
  `ensp` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Ensg` (`ensg`),
  KEY `Enst` (`enst`),
  KEY `Ensp` (`ensp`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Ensgv` (`ensgv`),
  KEY `Enstv` (`enstv`),
  KEY `Enspv` (`enspv`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt-ENSEMBL mappings (as provided by UniProt, and usually quite outdated)';

CREATE TABLE `unievo` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `acc` varchar(250) DEFAULT NULL,
  `species` varchar(250) DEFAULT NULL,
  `site` int DEFAULT NULL,
  `evorate` varchar(250) DEFAULT NULL,
  `rate` double DEFAULT NULL,
  `ptm` int DEFAULT NULL,
  `ptms` varchar(250) DEFAULT NULL,
  `hgmd` int DEFAULT NULL,
  `csa` int DEFAULT NULL,
  `int` int DEFAULT NULL,
  `intcore` int DEFAULT NULL,
  `protein_ptm` int DEFAULT NULL,
  `aa` char(1) DEFAULT NULL,
  `dis` char(1) DEFAULT NULL,
  `variable` int DEFAULT NULL,
  `gnomad` int DEFAULT NULL,
  `maf` double DEFAULT NULL,
  `maf_hc` double DEFAULT NULL,
  `ac` int DEFAULT NULL,
  `ac_hc` int DEFAULT NULL,
  `fast` int DEFAULT NULL,
  `aln_matching` int DEFAULT NULL,
  `aln_species` int DEFAULT NULL,
  `aln_column` varchar(250) DEFAULT NULL,
  `conserved_in_all` varchar(250) DEFAULT NULL,
  `present_in` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Acc` (`acc`),
  KEY `Aln_column` (`aln_column`),
  KEY `Aln_matching` (`aln_matching`),
  KEY `Aln_species` (`aln_species`),
  KEY `Conserved_in_all` (`conserved_in_all`),
  KEY `Evorate` (`evorate`),
  KEY `Hgmd` (`hgmd`),
  KEY `Present_in` (`present_in`),
  KEY `Ptm` (`ptm`),
  KEY `Site` (`site`),
  KEY `Species` (`species`),
  KEY `Variable` (`variable`),
  KEY `Fast` (`fast`),
  KEY `Aa` (`aa`),
  KEY `Protein_ptm` (`protein_ptm`),
  KEY `Ptms` (`ptms`),
  KEY `Dis` (`dis`),
  KEY `Gnomad` (`gnomad`),
  KEY `Ac` (`ac`),
  KEY `Ac_hc` (`ac_hc`),
  KEY `Int` (`int`),
  KEY `Intcore` (`intcore`),
  KEY `Csa` (`csa`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='output-human-rate4site_einsi_tree_1para-unique.txt';

CREATE TABLE `unievo_rate4site_para` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(250) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(250) DEFAULT NULL,
  `site` int DEFAULT NULL,
  `evorate` varchar(250) DEFAULT NULL,
  `rate` double DEFAULT NULL,
  `constraint` double DEFAULT NULL,
  `ptm` int DEFAULT NULL,
  `ptms` varchar(250) DEFAULT NULL,
  `hgmd` int DEFAULT NULL,
  `protein_ptm` int DEFAULT NULL,
  `aa` char(1) DEFAULT NULL,
  `dis` char(1) DEFAULT NULL,
  `csa` int DEFAULT NULL,
  `int` int DEFAULT NULL,
  `m2sg` int DEFAULT NULL,
  `clinvar` int DEFAULT NULL,
  `variable` int DEFAULT NULL,
  `exac` int DEFAULT NULL,
  `maf` double DEFAULT NULL,
  `allelecount` int DEFAULT NULL,
  `fast` int DEFAULT NULL,
  `aln_matching` int DEFAULT NULL,
  `aln_species` int DEFAULT NULL,
  `aln_column` text,
  `conserved_in_all` varchar(250) DEFAULT NULL,
  `present_in` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Site` (`site`),
  KEY `Evorate` (`evorate`),
  KEY `Ptm` (`ptm`),
  KEY `Hgmd` (`hgmd`),
  KEY `Variable` (`variable`),
  KEY `Fast` (`fast`),
  KEY `Clinvar` (`clinvar`),
  KEY `M2sg` (`m2sg`),
  KEY `Int` (`int`),
  KEY `Csa` (`csa`),
  KEY `Aa` (`aa`),
  KEY `Protein_ptm` (`protein_ptm`),
  KEY `Ptms` (`ptms`),
  KEY `Dis` (`dis`),
  KEY `Exac` (`exac`),
  KEY `Allelecount` (`allelecount`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='output-human-rate4site_para-unique.txt';

CREATE TABLE `unifeat` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `start` int DEFAULT NULL,
  `stop` int DEFAULT NULL,
  `type` varchar(50) DEFAULT NULL,
  `description` varchar(2500) DEFAULT NULL,
  `scale` varchar(250) DEFAULT NULL,
  `evidence` varchar(2500) DEFAULT NULL,
  `sources` varchar(2500) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`) USING BTREE,
  KEY `Scale` (`scale`) USING BTREE,
  KEY `Start` (`start`) USING BTREE,
  KEY `Stop` (`stop`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt misc features';

CREATE TABLE `unifunc` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `func` varchar(20000) DEFAULT NULL,
  `evidence` varchar(2500) DEFAULT NULL,
  `sources` varchar(2500) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Species` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt protein function descriptions';

CREATE TABLE `unigo` (
  `id` int NOT NULL AUTO_INCREMENT,
  `name` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `species` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `goid` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `gotype` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `goterm` varchar(250) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `evidencetype` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `evidence` varchar(250) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `project` varchar(250) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `newindex1` (`name`),
  KEY `newindex2` (`goid`),
  KEY `newindex3` (`gotype`),
  KEY `newindex4` (`goterm`),
  KEY `newindex5` (`evidence`),
  KEY `newindex6` (`evidencetype`),
  KEY `newindex7` (`project`),
  KEY `Acc` (`acc`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8mb3 COMMENT='GO terms for SwissProt proteins';

CREATE TABLE `unihom` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name1` varchar(25) DEFAULT NULL,
  `acc1` varchar(25) DEFAULT NULL,
  `ensp1` char(20) NOT NULL,
  `species1` varchar(25) DEFAULT NULL,
  `site1` int DEFAULT NULL,
  `aa1` char(1) DEFAULT NULL,
  `name2` varchar(25) DEFAULT NULL,
  `acc2` varchar(25) DEFAULT NULL,
  `ensp2` char(20) NOT NULL,
  `species2` varchar(25) DEFAULT NULL,
  `site2` int DEFAULT NULL,
  `aa2` char(1) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `alnsite` int DEFAULT NULL,
  `same` int DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name1` (`name1`),
  KEY `Name2` (`name2`),
  KEY `Species1` (`species1`),
  KEY `Species2` (`species2`),
  KEY `Site1` (`site1`),
  KEY `Site2` (`site2`),
  KEY `Ensp1` (`ensp1`),
  KEY `Ensp2` (`ensp2`),
  KEY `Aa1` (`aa1`),
  KEY `Aa2` (`aa2`),
  KEY `Same` (`same`),
  KEY `Acc1` (`acc1`),
  KEY `Acc2` (`acc2`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt residue-level homology (from Compara)';

CREATE TABLE `unihost` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(45) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `hosttax` varchar(25) DEFAULT NULL,
  `host` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Hosttax` (`hosttax`),
  KEY `Host` (`host`),
  KEY `Acc` (`acc`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt viral protein host species';

CREATE TABLE `uniid` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `type` varchar(250) DEFAULT NULL,
  `value` varchar(250) DEFAULT NULL,
  `valueshort` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Type` (`type`),
  KEY `Value` (`value`),
  KEY `Acc` (`acc`),
  KEY `Species` (`species`),
  KEY `Canon` (`canon`),
  KEY `Valueshort` (`valueshort`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt ID mapping (cross-references)';

CREATE TABLE `uniid_counts` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `type` varchar(25) DEFAULT NULL,
  `names` int DEFAULT NULL,
  `values` int DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Type` (`type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt ID mapping counts (for selecting best type)';

CREATE TABLE `uniinterpro` (
  `id` int NOT NULL AUTO_INCREMENT,
  `name` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `acc` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `species` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `domain` varchar(250) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `interproid` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Domain` (`domain`),
  KEY `Interproid` (`interproid`),
  KEY `Acc` (`acc`),
  KEY `Species` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8mb3 COMMENT='InterPro protein domains for SwissProt proteins';

CREATE TABLE `uniiso` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `name` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `event` varchar(250) DEFAULT NULL,
  `sequencetype` varchar(25) DEFAULT NULL,
  `alternative` int DEFAULT NULL,
  `confirmed` int DEFAULT NULL,
  `isoform` varchar(250) DEFAULT NULL,
  `isoform_aliases` varchar(250) DEFAULT NULL,
  `seq` text,
  PRIMARY KEY (`id`),
  KEY `Index_1` (`name`),
  KEY `Index_2` (`acc`),
  KEY `Index_3` (`species`),
  KEY `Index_4` (`sequencetype`),
  KEY `Index_6` (`event`),
  KEY `Canon` (`canon`),
  KEY `Alternative` (`alternative`),
  KEY `Confirmed` (`confirmed`),
  KEY `Isoform` (`isoform`),
  KEY `Isoform_aliases` (`isoform_aliases`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt isoform annotation (sequences from varsplic FASTA)';

CREATE TABLE `unikey` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `keyid` varchar(25) DEFAULT NULL,
  `keyword` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Acc` (`acc`),
  KEY `Species` (`species`),
  KEY `Keyid` (`keyid`),
  KEY `Keyword` (`keyword`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt keywords (these categorize and describe proteins)';

CREATE TABLE `uniloc` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `loc` varchar(255) DEFAULT NULL,
  `evidence` varchar(255) DEFAULT NULL,
  `sources` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Loc` (`loc`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Species` (`species`),
  KEY `Sources` (`sources`),
  KEY `Evidence` (`evidence`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt subcellular localisations';

CREATE TABLE `unimod` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) DEFAULT NULL,
  `acc` varchar(20) DEFAULT NULL,
  `canon` varchar(20) DEFAULT NULL,
  `species` varchar(10) DEFAULT NULL,
  `site` int DEFAULT NULL,
  `aa` char(1) DEFAULT NULL,
  `ptm` varchar(20) DEFAULT NULL,
  `type` varchar(50) DEFAULT NULL,
  `description` varchar(250) DEFAULT NULL,
  `source` varchar(25) DEFAULT NULL,
  `subset` varchar(25) DEFAULT NULL,
  `scale` varchar(10) DEFAULT NULL,
  `evid` varchar(250) DEFAULT NULL,
  `studies` smallint DEFAULT NULL,
  `studies_small` smallint DEFAULT NULL,
  `studies_large` smallint DEFAULT NULL,
  `studies_cst` smallint DEFAULT NULL,
  `pmids` varchar(10000) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`) USING BTREE,
  KEY `Scale` (`scale`) USING BTREE,
  KEY `Site` (`site`) USING BTREE,
  KEY `Ptm` (`ptm`) USING BTREE,
  KEY `AA` (`aa`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Subset` (`subset`),
  KEY `Evidence` (`evid`),
  KEY `Source` (`source`),
  KEY `Pmids` (`pmids`(100)),
  KEY `Studies` (`studies`),
  KEY `Studies_small` (`studies_small`),
  KEY `Studies_large` (`studies_large`),
  KEY `Studies_cst` (`studies_cst`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt PTMs';

CREATE TABLE `unimod_control` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` char(11) DEFAULT NULL,
  `acc` char(12) DEFAULT NULL,
  `species` char(5) DEFAULT NULL,
  `site` mediumint DEFAULT NULL,
  `aa` char(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Site` (`site`),
  KEY `AA` (`aa`),
  KEY `Acc` (`acc`),
  KEY `Species` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt PTM control sites';

CREATE TABLE `unimod_ochoa_full` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `uniprot` varchar(250) DEFAULT NULL,
  `position` int DEFAULT NULL,
  `residue` varchar(250) DEFAULT NULL,
  `mq_siteid` varchar(250) DEFAULT NULL,
  `best_pep` double DEFAULT NULL,
  `best_localization_prob` double DEFAULT NULL,
  `biological_samples` int DEFAULT NULL,
  `spectralcounts` int DEFAULT NULL,
  `is_disopred` varchar(250) DEFAULT NULL,
  `disopred_score` double DEFAULT NULL,
  `disprot_disorder` varchar(250) DEFAULT NULL,
  `coiled_coil_uniprot` int DEFAULT NULL,
  `comp_bias_uniprot` varchar(250) DEFAULT NULL,
  `domain_uniprot` varchar(250) DEFAULT NULL,
  `motif_uniprot` varchar(250) DEFAULT NULL,
  `region_uniprot` varchar(250) DEFAULT NULL,
  `repeat_uniprot` varchar(250) DEFAULT NULL,
  `zn_finger_uniprot` varchar(250) DEFAULT NULL,
  `elm_clv` varchar(250) DEFAULT NULL,
  `elm_deg` varchar(250) DEFAULT NULL,
  `elm_doc` varchar(250) DEFAULT NULL,
  `elm_lig` varchar(250) DEFAULT NULL,
  `elm_mod` varchar(250) DEFAULT NULL,
  `elm_trg` varchar(250) DEFAULT NULL,
  `iselmkinasemotif` varchar(250) DEFAULT NULL,
  `ev_mean_prediction_epistatic` varchar(250) DEFAULT NULL,
  `ev_mean_prediction_independent` double DEFAULT NULL,
  `ev_mean_frequency` double DEFAULT NULL,
  `ev_min_prediction_epistatic` double DEFAULT NULL,
  `ev_min_prediction_independent` double DEFAULT NULL,
  `ev_max_frequency` double DEFAULT NULL,
  `ev_ala_prediction_epistatic` double DEFAULT NULL,
  `ev_ala_prediction_independent` double DEFAULT NULL,
  `ev_ala_frequency` double DEFAULT NULL,
  `ev_acid_prediction_epistatic` double DEFAULT NULL,
  `ev_acid_prediction_independent` double DEFAULT NULL,
  `ev_acid_frequency` double DEFAULT NULL,
  `exac_max_ac_an_hom` double DEFAULT NULL,
  `exac_max_ac_an_adj` double DEFAULT NULL,
  `exac_ala_ac_an_adj` double DEFAULT NULL,
  `exac_acid_ac_an_adj` double DEFAULT NULL,
  `exp3d_max_ddg_prob` double DEFAULT NULL,
  `exp3d_mean_ddg_prob` double DEFAULT NULL,
  `exp3d_ala_ddg_prob` double DEFAULT NULL,
  `exp3d_acid_ddg_prob` double DEFAULT NULL,
  `exp3d_max_ddg_effect` varchar(250) DEFAULT NULL,
  `exp3d_mean_ddg_effect` varchar(250) DEFAULT NULL,
  `exp3d_ala_ddg_effect` varchar(250) DEFAULT NULL,
  `exp3d_acid_ddg_effect` varchar(250) DEFAULT NULL,
  `log10_hotspot_pval_min` varchar(250) DEFAULT NULL,
  `ishotspot` varchar(250) DEFAULT NULL,
  `intfc_max_mean_ddg_prob` varchar(250) DEFAULT NULL,
  `intfc_max_max_ddg_prob` double DEFAULT NULL,
  `intfc_max_ala_ddg_prob` double DEFAULT NULL,
  `intfc_max_acid_ddg_prob` double DEFAULT NULL,
  `isinterface` varchar(250) DEFAULT NULL,
  `adj_ptms_w21` varchar(250) DEFAULT NULL,
  `netpho_max_all` double DEFAULT NULL,
  `netpho_max_kin` double DEFAULT NULL,
  `netpho_max_ptp` double DEFAULT NULL,
  `netpho_max_stdomain` double DEFAULT NULL,
  `netpho_max_ydomain` double DEFAULT NULL,
  `paxdb_abundance_log10` double DEFAULT NULL,
  `paxdb_abundance` double DEFAULT NULL,
  `prot_length` double DEFAULT NULL,
  `w0_mya` double DEFAULT NULL,
  `w0_ancestor_name` varchar(250) DEFAULT NULL,
  `w3_mya` varchar(250) DEFAULT NULL,
  `w3_ancestor_name` varchar(250) DEFAULT NULL,
  `ptmdb_coreg_kinases` varchar(250) DEFAULT NULL,
  `ptmdb_coreg_kinases_len` varchar(250) DEFAULT NULL,
  `ptmdb_maxcoreg_kinase` varchar(250) DEFAULT NULL,
  `pubmed_counts` varchar(250) DEFAULT NULL,
  `quant_top1` varchar(250) DEFAULT NULL,
  `quant_top5` varchar(250) DEFAULT NULL,
  `quant_top10` int DEFAULT NULL,
  `pwm_max_mss` double DEFAULT NULL,
  `pwm_nkintop005` double DEFAULT NULL,
  `pwm_nkintop01` double DEFAULT NULL,
  `pwm_nkintop02` double DEFAULT NULL,
  `accpro` varchar(250) DEFAULT NULL,
  `sspro` varchar(250) DEFAULT NULL,
  `sspro8` varchar(250) DEFAULT NULL,
  `sift_min_score` varchar(250) DEFAULT NULL,
  `sift_mean_score` varchar(250) DEFAULT NULL,
  `sift_ala_score` varchar(250) DEFAULT NULL,
  `sift_acid_score` varchar(250) DEFAULT NULL,
  `topology_uniprot` varchar(250) DEFAULT NULL,
  `transitpep_uniprot` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Accpro` (`accpro`),
  KEY `Adj_ptms_w21` (`adj_ptms_w21`),
  KEY `Biological_samples` (`biological_samples`),
  KEY `Coiled_coil_uniprot` (`coiled_coil_uniprot`),
  KEY `Comp_bias_uniprot` (`comp_bias_uniprot`),
  KEY `Disprot_disorder` (`disprot_disorder`),
  KEY `Domain_uniprot` (`domain_uniprot`),
  KEY `Elm_clv` (`elm_clv`),
  KEY `Elm_deg` (`elm_deg`),
  KEY `Elm_doc` (`elm_doc`),
  KEY `Elm_lig` (`elm_lig`),
  KEY `Elm_mod` (`elm_mod`),
  KEY `Elm_trg` (`elm_trg`),
  KEY `Ev_mean_prediction_epistatic` (`ev_mean_prediction_epistatic`),
  KEY `Exp3d_acid_ddg_effect` (`exp3d_acid_ddg_effect`),
  KEY `Exp3d_ala_ddg_effect` (`exp3d_ala_ddg_effect`),
  KEY `Exp3d_max_ddg_effect` (`exp3d_max_ddg_effect`),
  KEY `Exp3d_mean_ddg_effect` (`exp3d_mean_ddg_effect`),
  KEY `Intfc_max_mean_ddg_prob` (`intfc_max_mean_ddg_prob`),
  KEY `Is_disopred` (`is_disopred`),
  KEY `Iselmkinasemotif` (`iselmkinasemotif`),
  KEY `Ishotspot` (`ishotspot`),
  KEY `Isinterface` (`isinterface`),
  KEY `Log10_hotspot_pval_min` (`log10_hotspot_pval_min`),
  KEY `Motif_uniprot` (`motif_uniprot`),
  KEY `Mq_siteid` (`mq_siteid`),
  KEY `Position` (`position`),
  KEY `Ptmdb_coreg_kinases` (`ptmdb_coreg_kinases`),
  KEY `Ptmdb_coreg_kinases_len` (`ptmdb_coreg_kinases_len`),
  KEY `Ptmdb_maxcoreg_kinase` (`ptmdb_maxcoreg_kinase`),
  KEY `Pubmed_counts` (`pubmed_counts`),
  KEY `Quant_top1` (`quant_top1`),
  KEY `Quant_top5` (`quant_top5`),
  KEY `Quant_top10` (`quant_top10`),
  KEY `Region_uniprot` (`region_uniprot`),
  KEY `Repeat_uniprot` (`repeat_uniprot`),
  KEY `Residue` (`residue`),
  KEY `Sift_acid_score` (`sift_acid_score`),
  KEY `Sift_ala_score` (`sift_ala_score`),
  KEY `Sift_mean_score` (`sift_mean_score`),
  KEY `Sift_min_score` (`sift_min_score`),
  KEY `Spectralcounts` (`spectralcounts`),
  KEY `Sspro` (`sspro`),
  KEY `Sspro8` (`sspro8`),
  KEY `Topology_uniprot` (`topology_uniprot`),
  KEY `Transitpep_uniprot` (`transitpep_uniprot`),
  KEY `Uniprot` (`uniprot`),
  KEY `W0_ancestor_name` (`w0_ancestor_name`),
  KEY `W3_ancestor_name` (`w3_ancestor_name`),
  KEY `W3_mya` (`w3_mya`),
  KEY `Zn_finger_uniprot` (`zn_finger_uniprot`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='input/unimod_ochoa/unimod_ochoa_table_s2_annotated_phosphoproteome_features.csv';

CREATE TABLE `unimod_ochoa_full_clinvar` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `alleleid` int DEFAULT NULL,
  `type` varchar(250) DEFAULT NULL,
  `name` varchar(250) DEFAULT NULL,
  `geneid` double DEFAULT NULL,
  `genesymbol` varchar(250) DEFAULT NULL,
  `hgnc_id` varchar(250) DEFAULT NULL,
  `clinicalsignificance` varchar(250) DEFAULT NULL,
  `clinsigsimple` varchar(250) DEFAULT NULL,
  `lastevaluated` varchar(250) DEFAULT NULL,
  `rs_dbsnp` varchar(250) DEFAULT NULL,
  `nsv_esv_dbvar` varchar(250) DEFAULT NULL,
  `rcvaccession` varchar(250) DEFAULT NULL,
  `phenotypeids` varchar(250) DEFAULT NULL,
  `phenotypelist` varchar(250) DEFAULT NULL,
  `origin` varchar(250) DEFAULT NULL,
  `originsimple` varchar(250) DEFAULT NULL,
  `assembly` varchar(250) DEFAULT NULL,
  `chromosomeaccession` varchar(250) DEFAULT NULL,
  `chromosome` varchar(250) DEFAULT NULL,
  `start` varchar(250) DEFAULT NULL,
  `stop` varchar(250) DEFAULT NULL,
  `referenceallele` varchar(250) DEFAULT NULL,
  `alternateallele` varchar(250) DEFAULT NULL,
  `refaa` varchar(250) DEFAULT NULL,
  `altaa` varchar(250) DEFAULT NULL,
  `uniprot_pos` varchar(250) DEFAULT NULL,
  `uniprot_id` varchar(250) DEFAULT NULL,
  `cytogenetic` varchar(250) DEFAULT NULL,
  `reviewstatus` text,
  `numbersubmitters` text,
  `guidelines` text,
  `testedingtr` text,
  `otherids` text,
  `submittercategories` text,
  `functional_score` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Alleleid` (`alleleid`),
  KEY `Altaa` (`altaa`),
  KEY `Alternateallele` (`alternateallele`),
  KEY `Assembly` (`assembly`),
  KEY `Chromosome` (`chromosome`),
  KEY `Chromosomeaccession` (`chromosomeaccession`),
  KEY `Clinicalsignificance` (`clinicalsignificance`),
  KEY `Clinsigsimple` (`clinsigsimple`),
  KEY `Cytogenetic` (`cytogenetic`),
  KEY `Functional_score` (`functional_score`),
  KEY `Genesymbol` (`genesymbol`),
  KEY `Hgnc_id` (`hgnc_id`),
  KEY `Lastevaluated` (`lastevaluated`),
  KEY `Name` (`name`),
  KEY `Nsv_esv_dbvar` (`nsv_esv_dbvar`),
  KEY `Origin` (`origin`),
  KEY `Originsimple` (`originsimple`),
  KEY `Phenotypeids` (`phenotypeids`),
  KEY `Phenotypelist` (`phenotypelist`),
  KEY `Rcvaccession` (`rcvaccession`),
  KEY `Refaa` (`refaa`),
  KEY `Referenceallele` (`referenceallele`),
  KEY `Rs_dbsnp` (`rs_dbsnp`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Type` (`type`),
  KEY `Uniprot_id` (`uniprot_id`),
  KEY `Uniprot_pos` (`uniprot_pos`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='input/unimod_ochoa/unimod_ochoa_table_s5_clinvar_variants.csv';

CREATE TABLE `unimod_ochoa_full_score` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `uniprot` varchar(250) DEFAULT NULL,
  `position` int DEFAULT NULL,
  `functional_score` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Position` (`position`),
  KEY `Uniprot` (`uniprot`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='input/unimod_ochoa/unimod_ochoa_table_s3_phosphosite_functional_scores.csv';

CREATE TABLE `unimod_ochoa_full_yeast_thermal` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `gene_name` varchar(250) DEFAULT NULL,
  `protein_id` varchar(250) DEFAULT NULL,
  `description` text,
  `found_in_ms_runs` varchar(250) DEFAULT NULL,
  `path` varchar(250) DEFAULT NULL,
  `score` varchar(250) DEFAULT NULL,
  `score_pby5_vs_wt` varchar(250) DEFAULT NULL,
  `score_tdh3_149a_doxo_vs_tdh3_149a` varchar(250) DEFAULT NULL,
  `score_tdh3_149a_doxo_vs_tdh3_149a_vs_tdh3_ko_doxo_vs_wt_doxo` varchar(250) DEFAULT NULL,
  `score_tdh3_149a_doxo_vs_tdh3_control_doxo` double DEFAULT NULL,
  `score_tdh3_149a_vs_tdh3_control` double DEFAULT NULL,
  `score_tdh3_149a_vs_tdh3_control_vs_tdh3_ko_vs_wt` double DEFAULT NULL,
  `score_tdh3_151a_doxo_vs_tdh3_151a` double DEFAULT NULL,
  `score_tdh3_151a_doxo_vs_tdh3_151a_vs_tdh3_ko_doxo_vs_wt_doxo` double DEFAULT NULL,
  `score_tdh3_151a_doxo_vs_tdh3_control_doxo` double DEFAULT NULL,
  `score_tdh3_151a_vs_tdh3_control` double DEFAULT NULL,
  `score_tdh3_151a_vs_tdh3_control_vs_tdh3_ko_vs_wt` double DEFAULT NULL,
  `score_tdh3_control_doxo_vs_tdh3_control` double DEFAULT NULL,
  `score_tdh3_control_vs_pby5` double DEFAULT NULL,
  `score_tdh3_ko_doxo_vs_tdh3_control_doxo` double DEFAULT NULL,
  `score_tdh3_ko_doxo_vs_tdh3_ko` double DEFAULT NULL,
  `score_tdh3_ko_doxo_vs_wt` double DEFAULT NULL,
  `score_tdh3_ko_doxo_vs_wt_doxo` double DEFAULT NULL,
  `score_tdh3_ko_vs_tdh3_control` double DEFAULT NULL,
  `score_tdh3_ko_vs_wt` double DEFAULT NULL,
  `fdr_pby5_vs_wt` double DEFAULT NULL,
  `fdr_tdh3_149a_doxo_vs_tdh3_149a` double DEFAULT NULL,
  `fdr_tdh3_149a_doxo_vs_tdh3_149a_vs_tdh3_ko_doxo_vs_wt_doxo` double DEFAULT NULL,
  `fdr_tdh3_149a_doxo_vs_tdh3_control_doxo` double DEFAULT NULL,
  `fdr_tdh3_149a_vs_tdh3_control` double DEFAULT NULL,
  `fdr_tdh3_149a_vs_tdh3_control_vs_tdh3_ko_vs_wt` double DEFAULT NULL,
  `fdr_tdh3_151a_doxo_vs_tdh3_151a` double DEFAULT NULL,
  `fdr_tdh3_151a_doxo_vs_tdh3_151a_vs_tdh3_ko_doxo_vs_wt_doxo` double DEFAULT NULL,
  `fdr_tdh3_151a_doxo_vs_tdh3_control_doxo` double DEFAULT NULL,
  `fdr_tdh3_151a_vs_tdh3_control` double DEFAULT NULL,
  `fdr_tdh3_151a_vs_tdh3_control_vs_tdh3_ko_vs_wt` double DEFAULT NULL,
  `fdr_tdh3_control_doxo_vs_tdh3_control` double DEFAULT NULL,
  `fdr_tdh3_control_vs_pby5` double DEFAULT NULL,
  `fdr_tdh3_ko_doxo_vs_tdh3_control_doxo` double DEFAULT NULL,
  `fdr_tdh3_ko_doxo_vs_tdh3_ko` double DEFAULT NULL,
  `fdr_tdh3_ko_doxo_vs_wt` double DEFAULT NULL,
  `fdr_tdh3_ko_doxo_vs_wt_doxo` double DEFAULT NULL,
  `fdr_tdh3_ko_vs_tdh3_control` double DEFAULT NULL,
  `fdr_tdh3_ko_vs_wt` double DEFAULT NULL,
  `hit_annotation_pby5_vs_wt` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_149a_doxo_vs_tdh3_149a` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_149a_doxo_vs_tdh3_149a_vs_tdh3_ko_doxo_vs_wt` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_149a_doxo_vs_tdh3_control_doxo` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_149a_vs_tdh3_control` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_149a_vs_tdh3_control_vs_tdh3_ko_vs_wt` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_151a_doxo_vs_tdh3_151a` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_151a_doxo_vs_tdh3_151a_vs_tdh3_ko_doxo_vs_wt` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_151a_doxo_vs_tdh3_control_doxo` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_151a_vs_tdh3_control` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_151a_vs_tdh3_control_vs_tdh3_ko_vs_wt` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_control_doxo_vs_tdh3_control` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_control_vs_pby5` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_ko_doxo_vs_tdh3_control_doxo` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_ko_doxo_vs_tdh3_ko` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_ko_doxo_vs_wt` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_ko_doxo_vs_wt_doxo` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_ko_vs_tdh3_control` varchar(250) DEFAULT NULL,
  `hit_annotation_tdh3_ko_vs_wt` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Found_in_ms_runs` (`found_in_ms_runs`),
  KEY `Gene_name` (`gene_name`),
  KEY `Hit_annotation_pby5_vs_wt` (`hit_annotation_pby5_vs_wt`),
  KEY `Hit_annotation_tdh3_149a_doxo_vs_tdh3_149a` (`hit_annotation_tdh3_149a_doxo_vs_tdh3_149a`),
  KEY `Hit_annotation_tdh3_149a_doxo_vs_tdh3_149a_vs_tdh3_ko_doxo_vs_wt` (`hit_annotation_tdh3_149a_doxo_vs_tdh3_149a_vs_tdh3_ko_doxo_vs_wt`),
  KEY `Hit_annotation_tdh3_149a_doxo_vs_tdh3_control_doxo` (`hit_annotation_tdh3_149a_doxo_vs_tdh3_control_doxo`),
  KEY `Hit_annotation_tdh3_149a_vs_tdh3_control` (`hit_annotation_tdh3_149a_vs_tdh3_control`),
  KEY `Hit_annotation_tdh3_149a_vs_tdh3_control_vs_tdh3_ko_vs_wt` (`hit_annotation_tdh3_149a_vs_tdh3_control_vs_tdh3_ko_vs_wt`),
  KEY `Hit_annotation_tdh3_151a_doxo_vs_tdh3_151a` (`hit_annotation_tdh3_151a_doxo_vs_tdh3_151a`),
  KEY `Hit_annotation_tdh3_151a_doxo_vs_tdh3_151a_vs_tdh3_ko_doxo_vs_wt` (`hit_annotation_tdh3_151a_doxo_vs_tdh3_151a_vs_tdh3_ko_doxo_vs_wt`),
  KEY `Hit_annotation_tdh3_151a_doxo_vs_tdh3_control_doxo` (`hit_annotation_tdh3_151a_doxo_vs_tdh3_control_doxo`),
  KEY `Hit_annotation_tdh3_151a_vs_tdh3_control` (`hit_annotation_tdh3_151a_vs_tdh3_control`),
  KEY `Hit_annotation_tdh3_151a_vs_tdh3_control_vs_tdh3_ko_vs_wt` (`hit_annotation_tdh3_151a_vs_tdh3_control_vs_tdh3_ko_vs_wt`),
  KEY `Hit_annotation_tdh3_control_doxo_vs_tdh3_control` (`hit_annotation_tdh3_control_doxo_vs_tdh3_control`),
  KEY `Hit_annotation_tdh3_control_vs_pby5` (`hit_annotation_tdh3_control_vs_pby5`),
  KEY `Hit_annotation_tdh3_ko_doxo_vs_tdh3_control_doxo` (`hit_annotation_tdh3_ko_doxo_vs_tdh3_control_doxo`),
  KEY `Hit_annotation_tdh3_ko_doxo_vs_tdh3_ko` (`hit_annotation_tdh3_ko_doxo_vs_tdh3_ko`),
  KEY `Hit_annotation_tdh3_ko_doxo_vs_wt` (`hit_annotation_tdh3_ko_doxo_vs_wt`),
  KEY `Hit_annotation_tdh3_ko_doxo_vs_wt_doxo` (`hit_annotation_tdh3_ko_doxo_vs_wt_doxo`),
  KEY `Hit_annotation_tdh3_ko_vs_tdh3_control` (`hit_annotation_tdh3_ko_vs_tdh3_control`),
  KEY `Hit_annotation_tdh3_ko_vs_wt` (`hit_annotation_tdh3_ko_vs_wt`),
  KEY `Path` (`path`),
  KEY `Protein_id` (`protein_id`),
  KEY `Score` (`score`),
  KEY `Score_pby5_vs_wt` (`score_pby5_vs_wt`),
  KEY `Score_tdh3_149a_doxo_vs_tdh3_149a` (`score_tdh3_149a_doxo_vs_tdh3_149a`),
  KEY `Score_tdh3_149a_doxo_vs_tdh3_149a_vs_tdh3_ko_doxo_vs_wt_doxo` (`score_tdh3_149a_doxo_vs_tdh3_149a_vs_tdh3_ko_doxo_vs_wt_doxo`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='input/unimod_ochoa/unimod_ochoa_table_s6_thermal_proteome_profiling.csv';

CREATE TABLE `uniorth` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name1` varchar(25) DEFAULT NULL,
  `ensp1` char(20) NOT NULL,
  `species1` varchar(25) DEFAULT NULL,
  `site1` int DEFAULT NULL,
  `aa1` char(1) DEFAULT NULL,
  `name2` varchar(25) DEFAULT NULL,
  `ensp2` char(20) NOT NULL,
  `species2` varchar(25) DEFAULT NULL,
  `site2` int DEFAULT NULL,
  `aa2` char(1) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `alnsite` int DEFAULT NULL,
  `same` int DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name1` (`name1`),
  KEY `Name2` (`name2`),
  KEY `Species1` (`species1`),
  KEY `Species2` (`species2`),
  KEY `Site1` (`site1`),
  KEY `Site2` (`site2`),
  KEY `Ensp1` (`ensp1`),
  KEY `Ensp2` (`ensp2`),
  KEY `Aa1` (`aa1`),
  KEY `Aa2` (`aa2`),
  KEY `Same` (`same`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt residue-level orthology (from Compara)';

CREATE TABLE `unipdb` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `pdb` char(4) DEFAULT NULL,
  `chains` varchar(250) DEFAULT NULL,
  `start` int DEFAULT NULL,
  `stop` int DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Pdb` (`pdb`),
  KEY `Acc` (`acc`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt-PDB mapping (official, from UniProt XML)';

CREATE TABLE `unipdb_sifts` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `site` int DEFAULT NULL,
  `pdb` char(4) DEFAULT NULL,
  `chain` char(1) DEFAULT NULL,
  `pdbsite` varchar(10) DEFAULT NULL,
  `method` varchar(25) DEFAULT NULL,
  `resolution` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Pdb` (`pdb`),
  KEY `Site` (`site`),
  KEY `Pdbsite` (`pdbsite`),
  KEY `Species` (`species`),
  KEY `Chain` (`chain`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt-PDB mapping (according to EBI SIFTS project)';

CREATE TABLE `unipfam` (
  `id` int NOT NULL AUTO_INCREMENT,
  `name` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `acc` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `canon` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `species` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `domain` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `pfamid` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `occ` int DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `newindex1` (`name`),
  KEY `newindex2` (`domain`),
  KEY `newindex3` (`pfamid`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8mb3 COMMENT='PFAM for SwissProt proteins';

CREATE TABLE `unipfam_pfam` (
  `id` int NOT NULL AUTO_INCREMENT,
  `name` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `pfamid` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `domain` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `description` varchar(250) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `clanid` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `clan` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `start` int DEFAULT NULL,
  `stop` int DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Domain` (`domain`),
  KEY `Pfamid` (`pfamid`),
  KEY `Start` (`start`),
  KEY `Stop` (`stop`),
  KEY `Clanid` (`clanid`),
  KEY `Clan` (`clan`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8mb3 COMMENT='PFAM for SwissProt from Pfam (outdated, but has positions)';

CREATE TABLE `uniprot` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `species` varchar(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `fullname` varchar(200) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `symbol` varchar(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `orderedlocus` varchar(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `orf` varchar(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `evidence` varchar(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `completeproteome` tinyint(1) DEFAULT NULL,
  `referenceproteome` tinyint(1) DEFAULT NULL,
  `trembl` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Orderedlocus` (`orderedlocus`),
  KEY `Orf` (`orf`),
  KEY `Trembl` (`trembl`),
  KEY `Completeproteome` (`completeproteome`),
  KEY `Referenceproteome` (`referenceproteome`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Fullname` (`fullname`),
  KEY `Acc` (`acc`),
  KEY `Symbol` (`symbol`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt database';

CREATE TABLE `uniproteome` (
  `id` int NOT NULL AUTO_INCREMENT,
  `name` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `species` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `proteome` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `newindex1` (`name`),
  KEY `newindex2` (`species`),
  KEY `newindex3` (`proteome`),
  KEY `Acc` (`acc`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8mb3 COMMENT='Reference proteome IDs for SwissProt proteins';

CREATE TABLE `unirefprot` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `species` varchar(25) DEFAULT NULL,
  `tax` int DEFAULT NULL,
  `proteome` varchar(25) DEFAULT NULL,
  `rel` varchar(10) DEFAULT NULL,
  `type` varchar(25) DEFAULT NULL,
  `sptr` char(2) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `name` varchar(25) DEFAULT NULL,
  `symbol` varchar(50) DEFAULT NULL,
  `fullname` varchar(250) DEFAULT NULL,
  `evidence` int unsigned DEFAULT NULL,
  `version` int unsigned DEFAULT NULL,
  `length` int unsigned DEFAULT NULL,
  `seq` text,
  PRIMARY KEY (`id`),
  KEY `Species` (`species`),
  KEY `Tax` (`tax`),
  KEY `Proteome` (`proteome`),
  KEY `Rel` (`rel`),
  KEY `Type` (`type`),
  KEY `Sptr` (`sptr`),
  KEY `Acc` (`acc`),
  KEY `Name` (`name`),
  KEY `Symbol` (`symbol`),
  KEY `Evidence` (`evidence`),
  KEY `Version` (`version`),
  KEY `Canon` (`canon`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt Reference Proteomes';

CREATE TABLE `uniseq` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(25) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `type` varchar(25) DEFAULT NULL,
  `fragments` varchar(8) DEFAULT NULL,
  `precursor` varchar(4) DEFAULT NULL,
  `trembl` tinyint(1) DEFAULT NULL,
  `seq` text,
  PRIMARY KEY (`id`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Seq` (`seq`(200)),
  KEY `Fragments` (`fragments`),
  KEY `Precursor` (`precursor`),
  KEY `Trembl` (`trembl`),
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`),
  KEY `Type` (`type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt sequences + disorder';

CREATE TABLE `unisnp` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `canon` varchar(25) DEFAULT NULL,
  `species` varchar(25) DEFAULT NULL,
  `start` int unsigned DEFAULT NULL,
  `stop` int unsigned DEFAULT NULL,
  `original` varchar(1000) DEFAULT NULL,
  `variant` varchar(1000) DEFAULT NULL,
  `note` text,
  `evidence` varchar(1000) DEFAULT NULL,
  `sources` varchar(1000) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Index_1` (`name`),
  KEY `Index_2` (`species`) USING BTREE,
  KEY `Index_3` (`start`) USING BTREE,
  KEY `Index_4` (`stop`) USING BTREE,
  KEY `Acc` (`acc`),
  KEY `Canon` (`canon`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt Natural Variants';

CREATE TABLE `unisupfam` (
  `id` int NOT NULL AUTO_INCREMENT,
  `name` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `acc` varchar(25) DEFAULT NULL,
  `species` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `domain` varchar(255) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `supfamid` varchar(25) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  `occ` int DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `newindex1` (`name`),
  KEY `newindex2` (`domain`),
  KEY `newindex3` (`supfamid`),
  KEY `Acc` (`acc`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8mb3 COMMENT='SUPERFAMILY for SwissProt proteins';

CREATE TABLE `unitax` (
  `id` int NOT NULL AUTO_INCREMENT,
  `tax` varchar(255) DEFAULT NULL,
  `species` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Tax` (`tax`),
  KEY `Species` (`species`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='UniProt NCBI Taxon ID to UniProt species table';
