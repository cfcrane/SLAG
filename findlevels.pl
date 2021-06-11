#!/usr/bin/env perl
use warnings;
#This script counts number of values of contig length over SLAG cycles to find instances of stepwise increase.
open (LIS, "< $ARGV[0]"); #allZeasubsetcontigsfilelists.txt
open (OUT, "> $ARGV[1]");
while ($file = <LIS>) {
  chomp $file;
  push(@files, $file);
}
for ($k = 0; $k < scalar(@files); $k++) {
  open (INP, "< $files[$k]");
  @matchinglongests = ();
  while ($cyclefile = <INP>) { #This file must be in chronological order of cycles.
    chomp $cyclefile;
    open (SEQ, "< $cyclefile");
    $matchinglongest = 0;
    $currlength = 0;
    while ($line = <SEQ>) {
      chomp $line;
      if ($line =~ m/>/) {
        if ($currlength > $matchinglongest) {$matchinglongest = $currlength;}
        $currlength = 0;
      }
      else {$currlength += length($line);}
    }
    if ($currlength > $matchinglongest) {$matchinglongest = $currlength;}
    push(@matchinglongests, $matchinglongest);
    close(SEQ);
  }
  @distinctlengths = ();
  push(@distinctlengths, $matchinglongests[0]);
  for ($i = 1; $i < scalar(@matchinglongests); $i++) {
    if ($matchinglongests[$i] != $matchinglongests[$i-1]) {
      $found = 0;
      for ($j = 0; $j < scalar(@distinctlengths); $j++) {
        if ($matchinglongests[$i] == $distinctlengths[$j]) {$found = 1;}
      }
      if ($found == 0) {push(@distinctlengths, $matchinglongests[$i]);}
    }
  }
  $nlevels = scalar(@distinctlengths);
  print OUT "$files[$k] has $nlevels distinct longest contig lengths.\n";
}
