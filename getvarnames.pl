#!/usr/bin/env perl
use warnings;
#This script prepared a list of all variable names that appear in SLAG's configuration files.
open (LIS, "< $ARGV[0]"); #sampleSLAGconfigfiles.txt
open (OUT, "> $ARGV[1]");
while ($file = <LIS>) {
  chomp $file;
  open (INP, "< $file");
  while ($line = <INP>) {
    if ($line =~ m/\$/) {
      chomp $line;
      @vars = split(/\s+/, $line);
      $names{$vars[0]} = $vars[2];
    }
  }
}
for $key (sort(keys(%names))) {print OUT "$key\t$names{$key}\n";}
