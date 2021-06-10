#!/usr/bin/perl
use warnings;
#This script attempts to identify the largest *subset.contigs file associated with a GenBank id, which is included in the file name.
if (scalar(@ARGV) != 5) {die "This script requires five arguments.\n";}
open (LIS, "< $ARGV[0]"); #alreadyrunaccessions11032014.txt
open (OUT, "> $ARGV[1]"); #full output
open (OUB, "> $ARGV[2]"); #succinct output
$dir = $ARGV[3]; #/scratch/lustreC/c/ccrane, without the terminal slash
$secstring = $ARGV[4]; #1023A
while ($line = <LIS>) {
  chomp $line;
  push(@accessions, $line);
}
for ($i = 0; $i < scalar(@accessions); $i++) {
  @files = ();
  @sizes = ();
  @files = <$dir/*$accessions[$i]*$secstring*subset.contigs>;
  $ncycles = scalar(@files); #THIS ASSUMES THAT THE FILES COME FROM ONE RUN! $secstring must distinguish this run from all others of the same accession.
  $maxsize = 0;
  for ($j = 0; $j < scalar(@files); $j++) {
    @filestats = stat($files[$j]);
    push(@sizes, $filestats[7]);
    if ($filestats[7] > $maxsize) {
      $maxsize = $filestats[7];
      $biggestfile = $files[$j];
    }
    print OUT "files[$j] = $files[$j] sizes[$j] = $sizes[$j]\n";
  }
  print OUT "Biggest one is $biggestfile at $maxsize\n";
  if ($maxsize > 1) {
    open (SEQ, "< $biggestfile");
    $ncontigs = 0;
    @clens = ();
    $totcontiglength = 0;
    $maxcontiglength = 0;
    $mincontiglength = 2000000000;
    $longestcontig = "NULL";
    $sequence = "";
    while ($line = <SEQ>) {
      chomp $line;
      if ($line =~ m/>/ || eof(SEQ)) {
        $contiglength = length($sequence);
        if ($contiglength > 3) {push(@clens, $contiglength);} #Don't push on the initial zero, it messes up calculation of the median.
        if ($contiglength > 0) {
          $ncontigs++;
          $totcontiglength += $contiglength;
          if ($contiglength > $maxcontiglength) {
            $maxcontiglength = $contiglength;
            $longestcontig = $defline;
          }
          if ($contiglength < $mincontiglength) {$mincontiglength = $contiglength;}
        }
        $defline = $line;
        $sequence = "";
      }
      else {$sequence .= $line;}
    }
    $meancontiglength = $totcontiglength / $ncontigs;
    @sortedclens = sort {$a <=> $b} @clens;
    for ($z = 0; $z < scalar(@sortedclens); $z++) {print OUT "DIAG for median: z = $z clens[$z] = $clens[$z] sortedclens[$z] = $sortedclens[$z]\n";}
    $k = int($ncontigs / 2);
    if ($ncontigs % 2 == 1) {$median = $sortedclens[$k];}
    else {$median = ($sortedclens[$k-1] + $sortedclens[$k]) / 2;}
    @revsortedclens = reverse(@sortedclens);
    $totclenslength = 0;
    for ($z = 0; $z < scalar(@clens); $z++) {$totclenslength += $clens[$z];}
    $halftotclenslength = $totclenslength / 2;
    $cumulength = 0;
    $n50 = "NULL";
    for ($z = 0; $z < scalar(@revsortedclens); $z++) {
      $cumulength += $revsortedclens[$z];
      if ($cumulength >= $halftotclenslength) {
        $n50 = $revsortedclens[$z];
        last;
      }
    }
    if ($ncontigs > 1) {
      $sdcontiglength = 0;
      for ($j = 0; $j < $ncontigs; $j++) {$sdcontiglength += ($clens[$j] - $meancontiglength)**2;}
      $sdcontiglength /= ($ncontigs - 1);
      $sdcontiglength = sqrt($sdcontiglength);
      $skewness = 0; 
      for ($j = 0; $j < $ncontigs; $j++) {$skewness += ($clens[$j] - $meancontiglength)**3;}
      $skewness /= (($ncontigs - 1) * $sdcontiglength**3);
    }
    else {
      $sdcontiglength = "undef";
      $skewness = "undef";
    }
    print OUT "For $biggestfile cycles = $ncycles ncontigs = $ncontigs shortest length = $mincontiglength mean length = $meancontiglength median length = $median n50 = $n50 longest length = $maxcontiglength stddev = $sdcontiglength skewness = $skewness longest id = $longestcontig\n";
    print OUB "$biggestfile\t$ncycles\t$ncontigs\t$mincontiglength\t$meancontiglength\t$median\t$n50\t$maxcontiglength\t$sdcontiglength\t$skewness\t$longestcontig\n";
  }
}
