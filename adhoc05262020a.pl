#!/usr/bin/env perl
use warnings;
#This script build SLAG table 3.
open (OUT, "> $ARGV[0]");
$increment = $ARGV[1]; #0.4
$limit = $ARGV[2]; #3.0
@files = ("TraesSRcap3stats0525.txt", "TraesSRphrapstats0525.txt", "TraesSRspadesstats0525.txt", "TraesLR1Acap3stats0525.txt", "TraesLR1phrapstats0525.txt", "TraesLR2cap3stats.txt", "TraesLR2phrapstats0525.txt");
@programs = ("cap3", "phrap", "SPADes", "cap3", "phrap", "cap3", "phrap");
@readlengths = ("2x150", "2x150", "2x150", "7-14kf", "7-14kf", "7-14kd", "7-14kd");
@readdepths = ("5x", "5x", "60x", "5x", "5x", "5x", "5x");
$nboundaries = 0;
print OUT "Program\tRead Length\tRead Depth";
$k = int($limit / $increment);
for ($i = 0; $i < $k; $i++) {
  $limits[$i] = 1 + $i * $increment;
  $nboundaries++;
  print OUT "\t$limits[$i]";
}
print OUT "\n";
for ($i = 0; $i < scalar(@files); $i++) {
  open (INP, "< $files[$i]");
  for ($j = 0; $j < $k; $j++) {$counts[$j] = 0;}
  while ($line = <INP>) {
    if ($line =~ m/=/) {next;}
    chomp $line;
    #TraesCS1A01G433600.1_SRcap3_subsetcontigsfiles.txt      21      438     1077    2.45890410958904        4
    @vars = split(/\t/, $line);
    $ratio = $vars[4];
    for ($j = 1; $j < $k; $j++) {
      if ($ratio >= $limits[$j-1]) {
        if ($ratio < $limits[$j]) {$counts[$j-1]++;}
      }
    }
    if ($ratio >= $limits[$k-1]) {$counts[$k-1]++;}
  }
  print OUT "$programs[$i]\t$readlengths[$i]\t$readdepths[$i]";
  for ($j = 0; $j < $k; $j++) {print OUT "\t$counts[$j]";}
  print OUT "\n";
}
