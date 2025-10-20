#!/usr/bin/env perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');



# Download
run("Download", "wget -O input/taxdump.tar.gz 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'");
chdir("input");
run("Unpack", "tar -xzvpf taxdump.tar.gz");
chdir("..");
run("Clean up", "rm -f input/taxdump.tar.gz");





# show directory
run("ls", "ls input");

done();
