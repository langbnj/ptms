#!/home/blang1/bin/perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $superreservedmysqldb;

our $usage = "$0 [input file] [table name] [-overwrite] [-nonulls] [-nonas] [-nominus] [-noquotes] [-nowhite] [-quiet] [-addheader] [-keepcomments] [-notype] [-double] [-utf8]\n\nImport data from a CSV or tab-separated file into a MySQL table, which gets created automatically.\n -overwrite: Remove MySQL table with this name first, if it already exists.\n -nonulls: Do NOT convert blank fields to NULL.\n -nonas: Do NOT convert NA or NULL fields to NULL.\n -nominus: Do NOT convert '-' fields to NULL.\n -noquotes: Remove double quotes around columns.\n -nowhite: Remove whitespace at the beginning and end of fields.\n -quiet: Don't show progress messages.\n -addheader: Artificially add a header of 'field1', 'field2', etc.\n -keepcomments: Keep comment lines (starting with '#' or '%'), e.g. for when the header line starts with '#' or '%'\n -notype: Leave all fields as text and don't add any indexes (much faster).\n -utf8: Use 'utf8mb4' encoding instead of the default 'latin1' (which is compatible with all existing tables).\n\nNote: The first line of the input file should be a header line with field names.\n\nNote: All columns other than float/double/text will be used as indexes.\n\nExample: $0 file.tsv table_name";
($infile, $table) = args(2);

our $superloudmysql = 1 if (!switch('quiet'));

if ($table =~ /(.+)\.(.+)/)
{
	Change($1);
	$table = $2;
}



# $infile = "input.txt";
# $outfile = "output.txt";

starttime();

# Compressed (.gz) yes/no?
$gz = 0;
# This is still broken (LOAD DATA INFILE won't work with a compressed file)
# if ($infile =~ /\.gz$/)
# {
# 	$gz = 1;
# }

# Open file
if ($gz == 1)
{
	# Compressed (.gz): Use zcat
	open(IN, "zcat $infile|") or die("\nError: Couldn't zcat '$infile'\n\n");
}
else
{
	# Uncompressed
	open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
}
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

die("Error: Input file is empty") if (!-s $infile);

# # Choose linebreak
# $/ = "\r\n";
# @a = <IN>;
# if (scalar(@a) == 1)
# {
# 	$/ = "\r";
# 	seek0();
# 	@a = <IN>;
# 	if (scalar(@a) == 1)
# 	{
# 		$/ = "\n";
# 		seek0();
# 		@a = <IN>;
# 		if (scalar(@a) == 1)
# 		{
# 			die("Error: Couldn't determine linebreak type in file '$infile'");
# 		}
# 	}
# }
# seek0();

# Get first 10 MB of file
# dd if=largefile count=10 bs=1M > largefile.10megsonly
# read(IN, $_, 10000000);
read(IN, $chunk, 10 * 1024 * 1024);

# Choose linebreak
# Try \r\n
$c = 0;
while ($chunk =~ m/\r\n/g) { $c++; }
if ($c > 1) { $/ = "\r\n"; }
else
{
	# Try \r
	$c = 0;
	while ($chunk =~ m/\r/g) { $c++; }
	if ($c > 1) { $/ = "\r"; }
	else
	{
		# Try \n
		$c = 0;
		while ($chunk =~ m/\n/g) { $c++; }
		if ($c > 1) { $/ = "\n"; }
		else
		{
			die("Error: Couldn't determine linebreak type in file '$infile'");
		}
	}
}
seek0();

$linebreak = '';
if ($/ eq "\r\n")
{
	$linebreak = '\r\n';
}
elsif ($/ eq "\r")
{
	$linebreak = '\r';
}
elsif ($/ eq "\n")
{
	$linebreak = '\n';
}
else
{
	die("Error: Unknown linebreak type in file '$infile': \n['$/']");
}
state("Linebreak type: $linebreak") if (!switch('quiet'));

# start

$query = Query("SHOW TABLES LIKE '$table'");
if (Numrows($query) != 0)
{
	if (switch('overwrite'))
	{
		Query("DROP TABLE IF EXISTS `$table`");
		state("Dropped table '$table'") if (!switch('quiet'));
	}
	else
	{
		die("Error: Table '$table' already exists (use -overwrite to remove it first)");
	}
}

if (switch('addheader'))
{
	$headerlines = 0;
}
else
{
	$headerlines = 1;
}
if (!switch('keepcomments'))
{
	while ($_ = <IN>)
	{
		chomp;
		if (/^[#%]/)
		{
			state("Skipping comment line $headerlines: $_", 1) if (!switch('quiet'));
			$headerlines++;
			next;
		}
		else
		{
			last;
		}
	}
	state("Skipping $headerlines header lines in SQL import later", 1) if (!switch('quiet'));
}
else
{
	$_ = <IN>; chomp;
}
@fields = split(/\t/, $_);
if (scalar(@fields) == 1)
{
    @fields = splitcsv($_);
	# @fields = split(/,/);		# This also splits commas within quoted fields (not good)
	if (scalar(@fields) == 1)
	{
		die("Error: Format not CSV or tab-separated (can only find one column in table header)");
	}
	state("Assuming format: CSV") if (!switch('quiet'));
	$format = 'CSV';
}
else
{
	state("Assuming format: Tab-separated") if (!switch('quiet'));
	$format = 'TSV';
}
die if (($format ne 'CSV') and ($format ne 'TSV'));
@fields = removetrailingfields(@fields);

die("Error: First line in '$infile' appears empty (tried to get column headers from there)") if (scalar(@fields) == 0);

$doublequoted = 1;
if (switch('addheader'))
{
	$i = 1;
	foreach (@fields)
	{
		# Check whether ALL fields are surrounded by double quotes
		if (($_ !~ /^"/) or ($_ !~ /"$/))
		{
			$doublequoted = 0;
		}

		$_ = "field$i";
		$i++;
	}
	seek0();
}
else
{
	foreach (@fields)
	{
		# Check whether ALL fields are surrounded by double quotes
		if (($_ !~ /^"/) or ($_ !~ /"$/))
		{
		    state("Field NOT double-quoted: '$_'") if (!switch('quiet'));
			$doublequoted = 0;
		}

		# Make field names lower case and replace anything non-alphanumeric with underscores, then collapse successive underscores to one and remove leading and trailing underscores
		$_ = lc($_);
		s/[^\w]/_/g;
		s/_+/_/g;
		s/_+$//g;
		s/^_+//g;
		
		# Limit them to 64 characters in length (the MySQL maximum)
		if (length($_) > 64)
		{
			$trim = substr($_, 0, 64);
			state("Trimming column name '$_' to '$trim' (64 character maximum)") if (!switch('quiet'));
			$_ = $trim;
		}

		if ($_ eq '')
		{
			state("Fields: ".join(", ", @fields)) if (!switch('quiet'));
			die("Error: Field name can't be empty in file '$infile' (try -addheader if the first line isn't a header)");
		}
		if ($_ eq 'id')
		{
			state("Fields: ".join(", ", @fields)) if (!switch('quiet'));
			die("Error: Field name can't be 'id' in file '$infile' (it's reserved)");
		}
	}

	state("Fields: ".join(", ", @fields)) if (!switch('quiet'));
}
if ($doublequoted == 1)
{
	state("Stripping double quotes (\"\")") if (!switch('quiet'));
}

state("Creating table '$table'") if (!switch('quiet'));

$q = "CREATE TABLE `$table` (
";

$i = 0;
foreach $field (@fields)
{
	# Rename some problematic fields (for MySQL)
	if ($field eq 'set')
	{
		$field = 'myset';
	    state("Renamed field: 'set' to '$field'") if (!switch('quiet'));
	}
	# if ($field eq 'mode')
	# {
	# 	$field = 'mymode';
	#     state("Renamed field: 'mode' to '$field'") if (!switch('quiet'));
	# }
	
	# Always text to begin with
	$q .= "  `$field` text,\n";

	# Next field number
	$i++;
}

$q =~ s/,\n$//s;

state("", 1) if (!switch('quiet'));

if (switch('utf8'))
{
	# UTF8 encoding
	$q .= "
) ENGINE=MyISAM DEFAULT CHARSET=utf8mb4 COMMENT='$infile'";
}
else
{
	# Default encoding
	$q .= "
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='$infile'";
}


# state($q);
Query($q);

# state("Reading '$infile' and inserting into table '$table'...");
# starttime();
# $query = Query("LOAD DATA INFILE '".rel2abs($infile)."' INTO TABLE `$table` FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '$linebreak' IGNORE 1 LINES");
# done();
# state("Inserted ".Numrows($query)." rows");
# stoptime();

nl() if (switch('quiet'));
state("Reading '$infile' and inserting into table '$table'...", switch('quiet'));

# Get path from MySQL server side
$serverpath = `pwd`;
$serverpath =~ s/[\r\n]+$//;
# state("['$serverpath']");
# $serverpath =~ s/^\/home\/lang\/home\///;
# $serverpath =~ s/^\/Users\/lang\///;
$serverpath .= "/$infile";
# state("['$serverpath']");
# exit;

seek0();
# skip headers
<IN>;
$_ = <IN>;
chomp;

# Construct import query

$q = qq(LOAD DATA LOCAL INFILE '$serverpath' INTO TABLE `$table`);

if ($format eq 'CSV')
{
	# $q .= qq( FIELDS TERMINATED BY ',');
	$q .= qq( FIELDS TERMINATED BY ',' ENCLOSED BY '"' ESCAPED BY '"');
}
elsif ($format eq 'TSV')
{
	$q .= qq( FIELDS TERMINATED BY '\\t');
}

if ($doublequoted == 1)
{
	$q .= qq( ENCLOSED BY '"');
}

$q .= qq( LINES TERMINATED BY '$linebreak');

if (!switch('addheader'))
{
	# $q .= q( IGNORE 1 LINES);
	$q .= qq( IGNORE $headerlines LINES);
}

$query = Query($q);

state("Loaded ".commify(Numrows($query))." rows", switch('quiet'));

# Wait for table to appear
state("Waiting for table to appear...") if (!switch('quiet'));
$query = Query("SELECT `TABLE_NAME` FROM `information_schema`.`tables` WHERE `TABLE_SCHEMA`='$superreservedmysqldb' AND `TABLE_NAME`='$table'");
while (Numrows($query) == 0)
{
	sleep(0.5);
}
$query = Query("SELECT * FROM `$table` LIMIT 1");
if (Numrows($query) == 0)
{
	die("Error: Table '$table' is empty");
}
state("OK") if (!switch('quiet'));

if (!switch('nonulls'))
{
	startme("Converting blank fields to NULL", 0, scalar(@fields)) if (!switch('quiet'));
	$affected = 0;
	foreach $field (@fields)
	{
		if (!switch('debug'))
		{
			$query = Query("UPDATE `$table` SET `$field`=NULL WHERE `$field`=''");
		}
		else
		{
			# debug
			$query = Query("SELECT * FROM `$table` WHERE `$field`=''");
		}

		$affected += Numrows($query);

		if (Numrows($query) > 0)
		{
			addme("fields affected", $field);
		}

		stepme(1) if (!switch('quiet'));
	}
	stopme() if (!switch('quiet'));

	state(commify($affected)." rows affected") if (!switch('quiet'));
	showme("fields affected") if (!switch('quiet'));
	clearme("fields affected");
}

Optimize($table, switch('quiet'));

if (!switch('nonas'))
{
	# startme("Converting NA and N/A fields (from R) to NULL", 0, scalar(@fields)) if (!switch('quiet'));
	startme("Converting NA, N/A and 'NULL' fields (from R) to NULL", 0, scalar(@fields)) if (!switch('quiet'));
	$affected = 0;
	foreach $field (@fields)
	{
		if (!switch('debug'))
		{
			# $query = Query("UPDATE `$table` SET `$field`=NULL WHERE `$field`='NA' OR `$field`='N/A'");
			$query = Query("UPDATE `$table` SET `$field`=NULL WHERE `$field`='NA' OR `$field`='N/A' OR `$field`='NULL'");
		}
		else
		{
			# debug
			# $query = Query("SELECT * FROM `$table` WHERE `$field`='NA' OR `$field`='N/A'");
			$query = Query("SELECT * FROM `$table` WHERE `$field`='NA' OR `$field`='N/A' OR `$field`='NULL'");
		}

		$affected += Numrows($query);

		if (Numrows($query) > 0)
		{
			addme("fields affected", $field);
		}

		stepme(1) if (!switch('quiet'));
	}
	stopme() if (!switch('quiet'));

	state(commify($affected)." rows affected") if (!switch('quiet'));
	showme("fields affected") if (!switch('quiet'));
	clearme("fields affected");
}

Optimize($table, switch('quiet'));

if (!switch('nominus'))
{
	startme("Converting '-' fields to NULL", 0, scalar(@fields)) if (!switch('quiet'));
	$affected = 0;
	foreach $field (@fields)
	{
		if (!switch('debug'))
		{
			$query = Query("UPDATE `$table` SET `$field`=NULL WHERE `$field`='-'");
		}
		else
		{
			# debug
			$query = Query("SELECT * FROM `$table` WHERE `$field`='-'");
		}

		$affected += Numrows($query);

		if (Numrows($query) > 0)
		{
			addme("fields affected", $field);
		}

		stepme(1) if (!switch('quiet'));
	}
	stopme() if (!switch('quiet'));

	state(commify($affected)." rows affected") if (!switch('quiet'));
	showme("fields affected") if (!switch('quiet'));
	clearme("fields affected");
}

Optimize($table, switch('quiet'));

if (switch('noquotes'))
{
	startme("Removing double quotes around values", 0, scalar(@fields)) if (!switch('quiet'));
	$affected = 0;
	foreach $field (@fields)
	{
		if (!switch('debug'))
		{
			# $query = Query("UPDATE `$table` SET `$field`=REGEXP_REPLACE(`$field`, '^\"|\"\$', '') WHERE `$field` REGEXP '^\"' AND `$field` REGEXP '\"\$'");
			# Multiple
			# $query = Query("UPDATE `$table` SET `$field`=REGEXP_REPLACE(`$field`, '^\"+|\"+\$', '') WHERE `$field` REGEXP '^\"+' OR `$field` REGEXP '\"+\$'");
			while (1) {
				# Only replace paired double quotes at beginning and end of field
				$query = Query("UPDATE `$table` SET `$field`=REGEXP_REPLACE(`$field`, '^\"|\"\$', '') WHERE `$field` REGEXP '^\"' AND `$field` REGEXP '\"\$'");
				last if (Numrows($query) == 0);
			}
		}
		else
		{
			# debug
			# $query = Query("SELECT * FROM `$table` WHERE `$field` REGEXP '^\"' AND `$field REGEXP '\"\$'");
			$query = Query("SELECT * FROM `$table` WHERE `$field` REGEXP '^\"+' OR `$field REGEXP '\"+\$'");
		}

		$affected += Numrows($query);

		if (Numrows($query) > 0)
		{
			addme("fields affected", $field);
		}

		stepme(1) if (!switch('quiet'));
	}
	stopme() if (!switch('quiet'));

	state(commify($affected)." rows affected") if (!switch('quiet'));
	showme("fields affected") if (!switch('quiet'));
	clearme("fields affected");
}

Optimize($table, switch('quiet'));

if (switch('nowhite'))
{
	startme("Removing whitespace (including linebreaks) around values", 0, scalar(@fields)) if (!switch('quiet'));
	$affected = 0;
	foreach $field (@fields)
	{
		if (!switch('debug'))
		{
			$query = Query("UPDATE `$table` SET `$field`=REGEXP_REPLACE(`$field`, '^\\\\s+|\\\\s+\$', '') WHERE `$field` REGEXP '^\\\\s+' OR `$field` REGEXP '\\\\s+\$'");
		}
		else
		{
			# debug
			$query = Query("SELECT * FROM `$table` WHERE `$field` REGEXP '^\\\\s+' OR `$field REGEXP '\\\\s+\$'");
		}

		$affected += Numrows($query);

		if (Numrows($query) > 0)
		{
			addme("fields affected", $field);
		}

		stepme(1) if (!switch('quiet'));
	}
	stopme() if (!switch('quiet'));

	state(commify($affected)." rows affected") if (!switch('quiet'));
	showme("fields affected") if (!switch('quiet'));
	clearme("fields affected");
}

Optimize($table, switch('quiet'));

if (!switch('nonulls'))
{
	startme("Converting blank fields to NULL", 0, scalar(@fields)) if (!switch('quiet'));
	$affected = 0;
	foreach $field (@fields)
	{
		if (!switch('debug'))
		{
			$query = Query("UPDATE `$table` SET `$field`=NULL WHERE `$field`=''");
		}
		else
		{
			# debug
			$query = Query("SELECT * FROM `$table` WHERE `$field`=''");
		}

		$affected += Numrows($query);

		if (Numrows($query) > 0)
		{
			addme("fields affected", $field);
		}

		stepme(1) if (!switch('quiet'));
	}
	stopme() if (!switch('quiet'));

	state(commify($affected)." rows affected") if (!switch('quiet'));
	showme("fields affected") if (!switch('quiet'));
	clearme("fields affected");
}

Optimize($table, switch('quiet'));

# $q = "CREATE TABLE `$table` (
#   `id` int(10) unsigned NOT NULL auto_increment,
# ";

$q = "ALTER TABLE `$table`
ADD COLUMN `id` INT(10) UNSIGNED NOT NULL AUTO_INCREMENT FIRST,\n";


if (!switch('notype'))
{
	nl() if (switch('quiet'));
	state("Determining best field types and selecting indexes...", switch('quiet'));
	%fieldtype = ();
	# @indexes = ();
	# @badindexes = ();
	@allindexes = ();
	$i = 0;
	foreach $field (@fields)
	{
		# @values = ();

		# $mainquery = Query("SELECT * FROM `$table`");

		# startme("$field [1 of ".scalar(@fields)."]", 1, Numrows($mainquery)) if (!switch('quiet'));
		state(">> $field", 1) if (!switch('quiet'));

		$maxlen = 0;
		$numeric = 1;
		$double = 1;
		# $uniq = 0;
		# $baduniq = 0;

		# if ((defined($a[$i])) and ($a[$i] =~ /[^\d]/))
		$query = Query("SELECT 1 FROM `$table` WHERE `$field` REGEXP '[^0-9]' LIMIT 1");
		if (Numrows($query) != 0)
		{
			$numeric = 0;
		}

		# if ((defined($a[$i])) and (($a[$i] =~ /[^\d\.\-]/) or ($a[$i] =~ /^.+-/)))
		$query = Query("SELECT 1 FROM `$table` WHERE `$field` REGEXP '^[^0-9]+\$' OR `$field` REGEXP '[^\\-0-9eE.+]' OR `$field` REGEXP '-[^0-9]' OR `$field` REGEXP '^[eE]' OR `$field` REGEXP '[eE]\$' OR `$field` REGEXP '-[^eE]+-' OR `$field` REGEXP '[+-]\$' OR `$field` REGEXP '^[0-9][+-][0-9]\$' OR `$field` REGEXP '\\\\..*\\\\.' LIMIT 1");
		if (Numrows($query) != 0)
		{
			$double = 0;
		}

		# if ((defined($a[$i])) and (length($a[$i]) > $maxlen))
		$query = Query("SELECT MAX(LENGTH(`$field`)) FROM `$table`");
		# {
		# 	$maxlen = length($a[$i]);
		# }
		($maxlen) = FetchOne($query);
		$maxlen = 0 if (!defined($maxlen));

		# if (defined($a[$i]))
		# {
		# 	push(@values, $a[$i]);
		# }
		# else
		# {
		# 	push(@values, '');
		# }

		# # See if field content is fairly unique (# unique values >= 1% of #total values), and therefore informative as an index:
		# # if (scalar(unique(@values)) >= (scalar(@values) / 100))
		# # {
		# # 	$uniq = 1;
		# # }
		# $query = Query("SELECT COUNT(`$field`), COUNT(DISTINCT `$field`) FROM `$table`");
		# ($values, $uniquevalues) = FetchOne($query);
		# if ($uniquevalues >= ($values / 100))
		# {
		# 	$uniq = 1;
		# }
		# # See if field content is vaguely unique (# unique values >= 0.01% of #total values), and therefore informative as an index:
		# # elsif (scalar(unique(@values)) >= (scalar(@values) / 10000))
		# # {
		# # 	$baduniq = 1;
		# # }
		# elsif ($uniquevalues >= ($values / 10000))
		# {
		# 	$baduniq = 1;
		# }

		if ($maxlen == 0)
		{
			# Shouldn't happen
			$maxlen = 1;
		}

		# Choose 'best' field type
		if ($numeric == 1)
		{
			# Int
			# $q .= "  `$field` int(10) default NULL,\n";
			$q .= "CHANGE COLUMN `$field` `$field` INT(10) DEFAULT NULL,\n";
			$fieldtype{$field} = 'numeric';

			# # as index
			# if ($uniq == 1)
			# {
			# 	push(@indexes, $field);
			# }
			# 
			# # as bad index
			# if ($baduniq == 1)
			# {
			# 	push(@badindexes, $field);
			# }

			# as non-float/double/text index
			push(@allindexes, $field);
		}
		elsif ($double == 1)
		{
			# Double
			# $q .= "  `$field` double default NULL,\n";
			# if (switch('double'))
			# {
			# 	$q .= "CHANGE COLUMN `$field` `$field` DOUBLE DEFAULT NULL,\n";
			# }
			# else
			# {
				$q .= "CHANGE COLUMN `$field` `$field` DOUBLE DEFAULT NULL,\n";
			# }
			$fieldtype{$field} = 'double';
		}
		elsif ($maxlen <= 250)
		{
			# Varchar
			# $q .= "  `$field` varchar($maxlen) default NULL,\n";
			# $q .= "  `$field` varchar(250) default NULL,\n";
			$q .= "CHANGE COLUMN `$field` `$field` VARCHAR(250) DEFAULT NULL,\n";
			$fieldtype{$field} = 'varchar';

			# # as index
			# if ($uniq == 1)
			# {
			# 	push(@indexes, $field);
			# }
			# 
			# # as bad index
			# if ($baduniq == 1)
			# {
			# 	push(@badindexes, $field);
			# }

			# as non-float/double/text index
			push(@allindexes, $field);
		}
		else
		{
			# Text
			# $q .= "  `$field` text,\n";
			# $q .= "CHANGE COLUMN `$field` `$field` TEXT,\n";		# already text
			$fieldtype{$field} = 'text';
		}

		state("  >> $fieldtype{$field}", 1) if (!switch('quiet'));

		# Next field number
		$i++;
	}

	# $q .= "
	#   PRIMARY KEY  (`id`)
	# ";

	state("", 1) if (!switch('quiet'));

	# If there are no decent 1% unique indexes, get 0.01% unique indexes...or just use all (except float/double) if the switch is set
	# if (switch('allindexes') or switch('allindices'))
	# {
		# state("Using all indexes (-allindexes)", switch('quiet'));
		# @indexes = @fields;

		# state("Using all indexes (except float/double/text) (-allindexes)", switch('quiet'));
		state("Using all indexes (except float/double/text)", switch('quiet'));
		@indexes = @allindexes;
	# }
	# elsif ((@indexes == 0) or switch('moreindexes') or switch('moreindices'))
	# {
	# 	state("No ≥1% unique indexes -- using ≥0.01% unique indexes (bad)", switch('quiet'));
	# 	@indexes = @badindexes;
	# }
	# else
	# {
	# 	state("Using 10% unique indexes (good)", switch('quiet'));
	# }

	foreach $index (unique(@indexes))
	{
		$indextitle = ucfirst($index);

		state(">> $indextitle", 1) if (!switch('quiet'));

		if ($fieldtype{$index} ne 'text')
		{
			# $q .= ",
			#   KEY `$indextitle` (`$index`)";
			$q .= "ADD INDEX `$indextitle` (`$index`),\n";
		}
		else
		{
			# $q .= ",
			#   KEY `$indextitle` USING HASH (`$index`(200))";
			$q .= "ADD INDEX `$indextitle` (`$index`(200)),\n";
		}
	}
}

$q .= q( ADD PRIMARY KEY (`id`));

# $q =~ s/,\n$//s;


state("Changing field types and adding indexes...") if (!switch('quiet'));
Query($q);
done() if (!switch('quiet'));

nl() if (!switch('quiet'));
nl() if (!switch('quiet'));
nl() if (switch('quiet'));

stoptime();

Optimize($table);



# Reset file: seek(IN, 0, 0);
sub seek0
{
	our $gz;
	our $infile;
	# our IN;
	
	# Reset file
	if ($gz == 1)
	{
		# Compressed (.gz): Use zcat
		close(IN);
		open(IN, "zcat $infile|") or die("\nError: Couldn't zcat '$infile'\n\n");
	}
	else
	{
		# Uncompressed
		seek(IN, 0, 0);
	}
}