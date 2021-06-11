#!/usr/bin/env perl
use warnings;
#This script looks for peaks and holes in contig length over SLAG cycles.
open (LIS, "< $ARGV[0]"); #allZeasubsetcontigsfilelists.txt
open (OUT, "> $ARGV[1]");
$elevation = $ARGV[2]; #1.2
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
  $bottoms = 0;
  $tops = 0;
  for ($i = 2; $i < scalar(@matchinglongests); $i++) {
    if ($matchinglongests[$i] > $elevation * $matchinglongests[$i-1]) {
      if ($matchinglongests[$i-2] > $elevation * $matchinglongests[$i-1]) {$bottoms++;}
    }
    if ($matchinglongests[$i-1] > $elevation * $matchinglongests[$i-2]) {
      if ($matchinglongests[$i-1] > $elevation * $matchinglongests[$i]) {$tops++;}
    }
  }
  print OUT "$files[$k] bottoms = $bottoms tops = $tops\n";
}
