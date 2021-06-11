#!/usr/bin/env perl
use warnings;
#This script looks for peaks and holes in contig length over SLAG cycles.
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
    open (SEQ, "< $cyclefile") || die "Where is $cyclefile ?\n";
    $matchinglongest = 0;
    $currlength = 0;
    $nseqs = 0;
    while ($line = <SEQ>) {
      chomp $line;
      if ($line =~ m/>/) {
        if ($currlength > $matchinglongest) {$matchinglongest = $currlength;}
        $currlength = 0;
        $nseqs++;
      }
      else {$currlength += length($line);}
    }
    if ($currlength > $matchinglongest) {$matchinglongest = $currlength;}
    push(@matchinglongests, $matchinglongest);
    close(SEQ);
  }
  $increases = 0;
  $decreases = 0;
  $flats = 0;
  for ($i = 1; $i < scalar(@matchinglongests); $i++) {
    if ($matchinglongests[$i] > $matchinglongests[$i-1]) {$increases++;}
    elsif ($matchinglongests[$i] == $matchinglongests[$i-1]) {$flats++;}
    else {$decreases++;}
  }
  print OUT "$files[$k] decreases = $decreases flats = $flats increases = $increases\n";
  print OUT "$matchinglongests[0]";
  for ($i = 1; $i < scalar(@matchinglongests); $i++) {print OUT "\t$matchinglongests[$i]";}
  print OUT "\n";
}
