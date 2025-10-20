#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

startme("Filling table 'uniid_counts' with the number of associated names and values for each ID type in 'uniid'");
starttime();
$mainquery = Query("SELECT DISTINCT type FROM uniid ORDER BY type");
while (($type) = Fetch($mainquery))
{
	addme("types", $type);
	
	$query = Query("SELECT COUNT(DISTINCT name), COUNT(DISTINCT value) FROM uniid WHERE type='$type'");
	($names, $values) = FetchOne($query);
	
	Query("INSERT INTO uniid_counts SET type='$type', names=$names, `values`=$values");
	
	stepme(10);
}
stopme();
stoptime();

showmeall(1);

done();
