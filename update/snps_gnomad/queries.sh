#!/bin/bash

echo "Counts at different allele frequencies:"

echo "Common variants (>5%)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE af > 0.05;" -h
echo ">> 27,514 (AC 1-251,463)"

echo "Low-frequency variants (between 1 and 5%)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE af >= 0.01 AND af<=0.05;" -h
echo ">> 31,066 (AC 1-12,564)"

echo "Very low-frequency variants (between 0.1 and 1%)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE af >= 0.001 AND af<=0.01;" -h
echo ">> 120,486 (AC 1-2,514)"

echo "Rare variants (<1%)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE af < 0.01;" -h
echo ">> 9,412,578 (AC 1-2,514)"

echo "Very rare variants (<0.1%)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE af < 0.001;" -h
echo ">> 9,347,027 (AC 1-251)"

echo "Super-rare variants (<0.01%)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE af < 0.0001;" -h
echo ">> 9,147,048 (AC 1-25)"

echo "Tripletons and above (AC>=3)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE ac >= 3;" -h
echo ">> 2,777,135 (AC 3-avg 744)"

echo "Ultra-rare variants (<0.001%)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE af < 0.00001;" -h
echo ">> 6,885,555 (AC 1 or 2)"

echo "Homozygous doubletons (AC=2)"
time ~/scripts/query.pl "SELECT DISTINCT nhomalt FROM snps_gnomad WHERE ac = 2;" -h
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE ac = 2 AND nhomalt = 1;" -h
echo ">> 15,670"

echo "Heterozygous doubletons (AC=2)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE ac = 2 AND nhomalt = 0;" -h
echo ">> 1,989,592"

echo "Singletons (AC=1)"
time ~/scripts/query.pl "SELECT COUNT(DISTINCT canon, site, original, variant), MIN(ac), AVG(ac), MAX(ac) FROM snps_gnomad WHERE ac = 1;" -h
echo ">> 6,769,045"
