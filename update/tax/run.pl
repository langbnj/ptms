#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');



# download
# run("Download", "download.pl");

# run
run("Main", "main.pl");
run("Fill rank, parents and species type fields", "add_parents.pl");

done();
