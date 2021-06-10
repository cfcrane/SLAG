#!/usr/bin/env perl
use warnings;
#This script gathers statistics from .subset.contigs files or aTRAM2 xx files.
open (LIS, "< $ARGV[0]"); #filelist03082021a.txt
open (OUT, "> $ARGV[1]");
$cyclescompleted = 0;
while ($file = <LIS>) {
  chomp $file;
  push(@files, $file);
  $cyclescompleted++;
  open (INP, "< $file");
  while ($line = <INP>) {
    chomp $line;
    if ($line =~ m/>/) {
      if (exists($ncontigs{$file})) {$ncontigs{$file}++;}
      else {$ncontigs{$file} = 1;}
      ($name = $line) =~ s/>//;
      $contiglengths{$file}{$name} = 0;
    }
    else {
      $contiglengths{$file}{$name} += length($line);
    }
  }
}
$longestcontig = 0;
$longestmedian = 0;
$shortestmedian = 2**30;
$mostcontigs = 0;
$fewestcontigs = 2**30;
for $key (keys(%contiglengths)) {
  $totalbases{$key} = 0;
  @sizes = ();
  for $index (keys(%{$contiglengths{$key}})) {
    $totalbases{$key} += $contiglengths{$key}{$index};
    push(@sizes, $contiglengths{$key}{$index});
    if ($contiglengths{$key}{$index} > $longestcontig) {
      $longestcontig = $contiglengths{$key}{$index};
      $longestcontigfile = $key;
      $longestcontigname = $index;
    }
  }
  @sortedsizes = sort {$a <=> $b} @sizes;
  $m = int(scalar(@sizes) / 2);
  if (scalar(@sizes) % 2 == 1) {$mediansize{$key} = $sortedsizes[$m];}
  else {$mediansize{$key} = ($sortedsizes[$m-1] + $sortedsizes[$m]) / 2;}
  $longestinfile{$key} = $sortedsizes[-1];
  if ($mediansize{$key} > $longestmedian) {$longestmedian = $mediansize{$key};}
  if ($mediansize{$key} < $shortestmedian) {$shortestmedian = $mediansize{$key};}
  if ($ncontigs{$key} > $mostcontigs) {$mostcontigs = $ncontigs{$key};}
  if ($ncontigs{$key} < $fewestcontigs) {$fewestcontigs = $ncontigs{$key};}
}
print OUT "Number of cycles completed = $cyclescompleted\n";
print OUT "Longest contig length = $longestcontig for contig $longestcontigname in file $longestcontigfile\n";
print OUT "Shortest median contig length = $shortestmedian longest median length = $longestmedian\n";
print OUT "Fewest contigs = $fewestcontigs most contigs = $mostcontigs\n";
for ($i = 0; $i < scalar(@files); $i++) {
  print OUT "$files[$i]\t$ncontigs{$files[$i]}\t$mediansize{$files[$i]}\t$longestinfile{$files[$i]}\t$totalbases{$files[$i]}\n";
}
