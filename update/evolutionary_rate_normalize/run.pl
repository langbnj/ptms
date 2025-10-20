#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run

@types = ('rate4site', 'lichtarge', 'capra');
# @modes = ('linsi_tree', 'ginsi_tree', 'einsi_tree', 'linsi', 'ginsi', 'einsi', 'einsi_tree_1para', );
@modes = ('einsi_tree_1para', );
# @windows = (0);
@windows = (0, 1, 3, 5);

starttime();

cd('tmp');
foreach $type (@types)
{
    foreach $mode (@modes)
    {
        if ($type eq 'capra')
        {
            foreach $window (@windows)
            {
                run("Normalize", "~/scripts/qsub.sh ../main.pl $type$window\_$mode");
            }
        }
        else
        {
            run("Normalize", "~/scripts/qsub.sh ../main.pl $type\_$mode");
        }
    }
}
cd('..');

waitforjobs();

stoptime();

# # linsi_tree:
# run("Normalize", "main.pl capra0_linsi_tree");
# run("Normalize", "main.pl lichtarge_linsi_tree");
# run("Normalize", "main.pl rate4site_linsi_tree");
# 
# # ginsi_tree:
# run("Normalize", "main.pl capra0_ginsi_tree");
# run("Normalize", "main.pl lichtarge_ginsi_tree");
# run("Normalize", "main.pl rate4site_ginsi_tree");
# 
# # einsi_tree:
# run("Normalize", "main.pl capra0_einsi_tree");
# run("Normalize", "main.pl lichtarge_einsi_tree");
# run("Normalize", "main.pl rate4site_einsi_tree");

# # _para:
# run("Normalize", "main.pl capra0_para");
# run("Normalize", "main.pl lichtarge_para");
# run("Normalize", "main.pl rate4site_para");

done();
