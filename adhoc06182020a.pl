#!/usr/bin/env perl
use warnings;
#This script extracted the full strings for phrapsettings, etc., from a list of SLAG configuration files.
open (LIS, "< $ARGV[0]"); #sampleSLAGconfigfiles.txt
$searchstring = $ARGV[1]; #phrapsettings
while ($file = <LIS>) {
  chomp $file;
  open (INP, "< $file");
  while ($line = <INP>) {
    if ($line =~ m/$searchstring/) {print $line;}
  }
}
