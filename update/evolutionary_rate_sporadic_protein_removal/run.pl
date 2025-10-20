#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run
# run("Split tables", "split_tables.pl");   # No longer necessary (the single giant table, 'evorate', doesn't exist anymore - using individual tables now)
# run("Invert Capra et al. scores", "invert_capra.pl");     # No longer necessary - update/evorate_capra does this now
# run("Remove Capra et al. gaps", "remove_capra_gaps.pl");     # No longer necessary - update/evorate_capra checks for this now

# run("Remove sporadic proteins", "remove_sporadic_proteins.pl 'capra|lichtarge|rate4site' -debug");
run("Remove sporadic proteins", "remove_sporadic_proteins.pl 'capra|lichtarge|rate4site'");

# run("Remove sporadic proteins", "remove_sporadic_proteins.pl capra");
# run("Remove sporadic proteins", "remove_sporadic_proteins.pl lichtarge");
# run("Remove sporadic proteins", "remove_sporadic_proteins.pl rate4site");

# run("Remove sporadic proteins", "remove_sporadic_proteins.pl capra -para");
# run("Remove sporadic proteins", "remove_sporadic_proteins.pl lichtarge -para");
# run("Remove sporadic proteins", "remove_sporadic_proteins.pl rate4site -para");

done();
