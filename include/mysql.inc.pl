#!/home/blang1/bin/perl -w
use DBD::mysql;
use DBI;
use DBD::Proxy; # For TSV export

# Print out each MySQL query before running it?
# $superloudmysql = 1;

# Connect to default database
our $superreservedmysqldb = 'blang';
our $superreserverdmysqlhost = "[server name]"; # Insert MySQL server name
Connect($superreservedmysqldb);

sub Connect
{
	my ($database) = @_;

	instakill();
	
	$host = $superreserverdmysqlhost;

	our $superreservedmysqldb = $database;

	# If connection doesn't work straight away, retry a bit later (probably too many connections open right now)
	while (!defined($superreservedlink))
	{
		eval { our $superreservedlink = DBI -> connect("DBI:mysql:database=$database;host=$host;mysql_read_default_file=~/.my.cnf;mysql_multi_statements=1;mysql_local_infile=1", '[mysql account]') }; # Insert [mysql account] name
	}
    if ($@)
    {
        warn("Warning: ".$@);
        $superreservedlink = undef;
    }

	if (!defined($superreservedlink))
	{
		sleep(10);
	}
	# else
	# {
	# 	$superreservedlink -> selectdb($database) or $superreservedlink = undef;
	# }
}

sub Disconnect
{
	our $superreservedlink;

	if (defined($superreservedlink))
	{
		$superreservedlink -> disconnect();
		$superreservedlink = undef;
	}
}

sub Change
{
	my ($database) = @_;

	instakill();

	Disconnect();

	Connect($database);

	our $superreservedmysqldb = $database;
}

sub Query
{
    my ($query) = @_;

    our $superreservedlink;

    instakill();

    print "$query\n" if (defined($superloudmysql));

	$result = $superreservedlink -> prepare($query);
	$result -> execute;

    # my $error = $superreservedlink -> errmsg();
    my $error = $superreservedlink -> {'mysql_error'};
    if ($error ne '')
    {
        # die("MySQL error: $error");
        # die("Error: MySQL query failed!\n\n");
        die("Error: MySQL query failed:\n\n$query\n\n");
    }


    return $result;
}

sub Megaquery
{
    # Run a (huge) query without binary logging (so as not to waste disk space)

	my ($query) = @_;

    Query("SET sql_log_bin=0");

    my $result = Query($query);

    Query("SET sql_log_bin=1");

    return $result;
}

sub Fetch
{
	my ($link) = @_;

	instakill();

	# return $link -> fetchrow;
	return $link -> fetchrow_array;
}

sub FetchOne
{
	my ($link) = @_;

	instakill();

	# if (Numrows($$1)() != 1)
	if ($link -> rows != 1)
	{
		die("Error: Numrows ".$link -> rows." instead of 1\n");
	}

	# return $link -> fetchrow;
	return $link -> fetchrow_array;
}

sub DumpTable
{
	my ($table) = @_;

	state("Dumping table '$table':");

	my $query = Query("SELECT * FROM $table");
	while (my @a = Fetch($query))
	{
		print join("\t", @a)."\n";
	}

	nl();
}
sub Numrows
{
	my ($link) = @_;

	instakill();

	# return Numrows($$1);
	return $link -> rows;
}

sub Clear
{
	my ($table) = @_;

	instakill();

	$superreservedlink -> do("TRUNCATE $table;");

	print "\nCleared table '$table'\n";
}

sub Optimize
{
	my ($table, $silent) = @_;

	instakill();

	if (!defined($silent)) { $silent = 0; }

	if ($silent == 0)
	{
		print "Optimizing table '$table'...\n";
	}

	$superreservedlink -> do("ANALYZE TABLE $table;");
	$superreservedlink -> do("OPTIMIZE TABLE $table;");

	if ($silent == 0)
	{
		print "Done!\n\n";
	}
}

sub Exists
{
	my ($table) = @_;
	
	my $query = Query("SHOW TABLES LIKE '$table'");
	
	return Numrows($query);
}

return 1;
