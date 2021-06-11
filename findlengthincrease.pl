#!/usr/bin/env perl
use warnings;
#This script finds the ratio of longest contig in any cycle to longest contig in first cycle for SLAG table 3.
open (LST, "< $ARGV[0]"); #TraesSRcap3listofsubsetcontigslists.txt
open (OUT, "> $ARGV[1]");
$maxratio = $ARGV[2]; #1.6
$increment = $ARGV[3]; #0.2
$limit = int($maxratio / $increment) + 1;
for ($i = 0; $i < $limit; $i++) {$brackets[$i] = 1 + $i * $increment;}
for ($i = 0; $i < $limit - 1; $i++) {$totals[$i] = 0;} #n - 1 intervals for n bracket boundaries
$nfiles = 0;
while ($list = <LST>) {
  chomp $list;
  open (LIS, "< $list"); #Each list MUST be in chronological order of cycles.
  $k = 0;
  $longest = 0;
  $firstlongest = 0;
  $nfiles++;
  while ($file = <LIS>) {
    chomp $file;
    open (SEQ, "< $file");
    $currlongest = 0;
    $currlength = 0;
    while ($line = <SEQ>) {
      if ($line =~ m/>/) {
        if ($currlength > $currlongest) {$currlongest = $currlength;}
        $currlength = 0;
      }
      else {
        chomp $line;
        $currlength += length($line);
      }
    }
    if ($currlength > $currlongest) {$currlongest = $currlength;}
    if ($k == 0) {$firstlongest = $currlongest;}
    if ($currlongest > $firstlongest) {$longest = $currlongest;}
    $k++;
  }
  $ratio = $longest / $firstlongest;
  for ($i = 1; $i < scalar(@brackets); $i++) {
    if ($ratio >= $brackets[$i-1] && $ratio < $brackets[$i]) {$totals[$i-1]++;}
  }
}
print OUT "There were $nfiles files.\n";
for ($i = 0; $i < scalar(@brackets); $i++) {
  if ($i == 0) {print OUT $brackets[0];}
  else {print OUT "\t$brackets[$i]";}
}
print OUT "\n";
$sum = 0;
for ($i = 0; $i < scalar(@totals); $i++) {
  $sum += $totals[$i];
  if ($i == 0) {print OUT $totals[0];}
  else {print OUT "\t$totals[$i]";}
}
$v = $nfiles - $sum;
print OUT "\t$v";
