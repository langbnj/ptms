#!/home/blang1/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialise

$minsites = 1000;   # Minimum number of sites a PTM type needs to have in table 'unimod' to be analysed
# $minsites = 100;  # Minimum number of sites a PTM type needs to have in table 'unimod' to be analysed

# start

# @uniseq_types = ('AlphaFold');
# @uniseq_types = ('AlphaFold_pLDDT70', 'AlphaFold_pLDDT90');
# @uniseq_types = ('AlphaFold_pLDDT70_PAE2', 'AlphaFold_pLDDT90_PAE1');
@uniseq_types = ('AlphaFold_pLDDT70_PAE2');

# @sources = ('all', 'UniProt', 'Ochoa', 'PhosphoSitePlus', 'dbPTM');
@sources = ('UniProt', 'Ochoa', 'PhosphoSitePlus', 'dbPTM');
# @sources = ('all');
# @modes = ('einsi_tree_1para', 'einsi_tree', 'ginsi_tree', 'linsi_tree', 'einsi', 'ginsi', 'linsi');
# @modes = ('einsi_tree', 'ginsi_tree', 'linsi_tree', 'einsi', 'ginsi', 'linsi');
@modes = ('einsi_tree_1para');
# @types = ('rate4site', 'lichtarge', 'capra', 'norm_rate4site', 'norm_lichtarge', 'norm_capra');
# @types = ('norm_rate4site', 'norm_lichtarge', 'norm_capra');
# @types = ('rate4site', 'lichtarge', 'capra');
@types = ('rate4site');
# @windows = (0, 1, 3, 5);
@windows = (0);

@drawparams = (
    # Use existing statistics (p-values etc.)
    # '',
	'-coresurf',
	'-dismerge',
	'-surfmerge',
	# '-nosurf',
	# '-filtered',

    # # Recalculate statistics (p-values etc.)
    # # '-recalc',
	# '-recalc -coresurf',
	# '-recalc -dismerge',
	# '-recalc -surfmerge',
	# # '-recalc -nosurf',
	# # '-recalc -filtered',
);

# Number of bootstrap resamples
$boots = 10000;
# $boots = 1000;
# $boots = 100;
# $boots = 1;

# @types = ('rate4site');
# @modes = ('linsi_tree');
# @windows = (0);

starttime();

waitforjobs();

foreach $uniseq_type (@uniseq_types)
{
    foreach $source (@sources)
    {
        foreach $mode (@modes)
        {
            foreach $type (@types)
            {
                if (($type eq 'capra') or ($type eq 'norm_capra'))
                {
                    foreach $window (@windows)
                    {
                        # run("Analysis", "~/scripts/qsub.sh analysis.pl human $source $type$window\_$mode AlphaFold $minsites");
                        run("Analysis", "~/scripts/qsub.sh analysis_coresurf.pl human $source $type$window\_$mode $uniseq_type $minsites");
                    }
                }
                else
                {
                    # run("Analysis", "~/scripts/qsub.sh analysis.pl human $source $type\_$mode AlphaFold $minsites");
                    run("Analysis", "~/scripts/qsub.sh analysis_coresurf.pl human $source $type\_$mode $uniseq_type $minsites");
                }
            }
            # waitforjobs();
        }
        # waitforjobs();
    }
}

# run("Analysis", "~/scripts/qsub.sh analysis.pl human capra0_einsi_tree AlphaFold $minsites");
# run("Analysis", "~/scripts/qsub.sh analysis.pl human lichtarge_einsi_tree AlphaFold $minsites");
# run("Analysis", "~/scripts/qsub.sh analysis.pl human rate4site_einsi_tree AlphaFold $minsites");
# 
# run("Analysis", "~/scripts/qsub.sh analysis.pl human norm_capra0_einsi_tree AlphaFold $minsites");
# run("Analysis", "~/scripts/qsub.sh analysis.pl human norm_lichtarge_einsi_tree AlphaFold $minsites");
# run("Analysis", "~/scripts/qsub.sh analysis.pl human norm_rate4site_einsi_tree AlphaFold $minsites");
# 
# run("Analysis", "~/scripts/qsub.sh analysis_coresurf.pl human capra0_einsi_tree AlphaFold $minsites");
# run("Analysis", "~/scripts/qsub.sh analysis_coresurf.pl human lichtarge_einsi_tree AlphaFold $minsites");
# run("Analysis", "~/scripts/qsub.sh analysis_coresurf.pl human rate4site_einsi_tree AlphaFold $minsites");
# 
# run("Analysis", "~/scripts/qsub.sh analysis_coresurf.pl human norm_capra0_einsi_tree AlphaFold $minsites");
# run("Analysis", "~/scripts/qsub.sh analysis_coresurf.pl human norm_lichtarge_einsi_tree AlphaFold $minsites");
# run("Analysis", "~/scripts/qsub.sh analysis_coresurf.pl human norm_rate4site_einsi_tree AlphaFold $minsites");

waitforjobs();
stoptime();
starttime();

# # Draw
# foreach $drawparams (@drawparams)
# {
#     foreach $source (@sources)
#     {
#         foreach $type (@types)
#         {
#             foreach $mode (@modes)
#             {
#                 if (($type eq 'capra') or ($type eq 'norm_capra'))
#                 {
#                     foreach $window (@windows)
#                     {
#                         run("Analysis", "~/scripts/qsub.sh draw.pl human $source $type$window\_$mode AlphaFold $boots $drawparams");
#                     }
#                 }
#                 else
#                 {
#                     run("Analysis", "~/scripts/qsub.sh draw.pl human $source $type\_$mode AlphaFold $boots $drawparams");
#                 }
#             }
#         }
#     }
# }
# 
# # run("Draw", "~/scripts/qsub.sh draw.pl human capra0_einsi_tree AlphaFold");
# # run("Draw", "~/scripts/qsub.sh draw.pl human lichtarge_einsi_tree AlphaFold");
# # run("Draw", "~/scripts/qsub.sh draw.pl human rate4site_einsi_tree AlphaFold");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_capra0_einsi_tree AlphaFold");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_lichtarge_einsi_tree AlphaFold");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_rate4site_einsi_tree AlphaFold");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human capra0_einsi_tree AlphaFold -dismerge");
# # run("Draw", "~/scripts/qsub.sh draw.pl human lichtarge_einsi_tree AlphaFold -dismerge");
# # run("Draw", "~/scripts/qsub.sh draw.pl human rate4site_einsi_tree AlphaFold -dismerge");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_capra0_einsi_tree AlphaFold -dismerge");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_lichtarge_einsi_tree AlphaFold -dismerge");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_rate4site_einsi_tree AlphaFold -dismerge");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human capra0_einsi_tree AlphaFold -surfmerge");
# # run("Draw", "~/scripts/qsub.sh draw.pl human lichtarge_einsi_tree AlphaFold -surfmerge");
# # run("Draw", "~/scripts/qsub.sh draw.pl human rate4site_einsi_tree AlphaFold -surfmerge");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_capra0_einsi_tree AlphaFold -surfmerge");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_lichtarge_einsi_tree AlphaFold -surfmerge");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_rate4site_einsi_tree AlphaFold -surfmerge");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human capra0_einsi_tree AlphaFold -nosurf");
# # run("Draw", "~/scripts/qsub.sh draw.pl human lichtarge_einsi_tree AlphaFold -nosurf");
# # run("Draw", "~/scripts/qsub.sh draw.pl human rate4site_einsi_tree AlphaFold -nosurf");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_capra0_einsi_tree AlphaFold -nosurf");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_lichtarge_einsi_tree AlphaFold -nosurf");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_rate4site_einsi_tree AlphaFold -nosurf");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human capra0_einsi_tree AlphaFold -filtered");
# # run("Draw", "~/scripts/qsub.sh draw.pl human lichtarge_einsi_tree AlphaFold -filtered");
# # run("Draw", "~/scripts/qsub.sh draw.pl human rate4site_einsi_tree AlphaFold -filtered");
# # 
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_capra0_einsi_tree AlphaFold -filtered");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_lichtarge_einsi_tree AlphaFold -filtered");
# # run("Draw", "~/scripts/qsub.sh draw.pl human norm_rate4site_einsi_tree AlphaFold -filtered");
# 
# waitforjobs();
# stoptime();
# 
# done();
