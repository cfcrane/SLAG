#!/usr/bin/env perl
use warnings;
#This script examines contig lengths and accuracy over a set of cycles.
open (LIA, "< $ARGV[0]"); #ZeaferredoxinsSRcap3subsetcontigsfiles.txt, presumed to be in chronological order
open (LIB, "< $ARGV[1]"); #ZeaferredoxinsSRcap3basecontigsfiles.txt, presumed to be in chronological order
open (OUT, "> $ARGV[2]"); #ZeaferredoxinsSRcap3subsetcontigseries3.txt
$rscript = $ARGV[3]; #ZeaferredoxinsSRcap3graphics3.R
$pdffile = $ARGV[4]; #ZeaferredoxinsSRcap3graph3.pdf
$pdfheight = $ARGV[5]; #3
$pdfwidth = $ARGV[6]; #4
while ($file = <LIA>) { #LIA MUST be in chronological order of cycles.
  chomp $file;
  open (SEQ, "< $file");
  $matchingcontigs = 0;
  $matchinglength = 0;
  $matchinglongest = 0;
  $currlength = 0;
  while ($line = <SEQ>) {
    chomp $line;
    if ($line =~ m/>/) {
      $matchinglength += $currlength;
      if ($currlength > $matchinglongest) {$matchinglongest = $currlength;}
      $currlength = 0;
      $matchingcontigs++;
    }
    else {$currlength += length($line);}
  }
  if ($currlength > $matchinglongest) {$matchinglongest = $currlength;}
  $matchinglength += $currlength;
  $meanmatchinglength = $matchinglength / $matchingcontigs;
  push(@rmatchingcontigs, $matchingcontigs);
  push(@rmatchinglengths, $matchinglength);
  push(@rmeanmatchinglengths, $meanmatchinglength);
  push(@rmatchinglongests, $matchinglongest);
  close(SEQ);
}
while ($file = <LIB>) { #LIB MUST be in chronological order of cycles.
  chomp $file; 
  open (SEQ, "< $file");
  $allcontigs = 0;
  $alllength = 0;
  $currlength = 0;
  while ($line =  <SEQ>) {
    chomp $line;
    if ($line =~ m/>/) {
      $alllength += $currlength;
      $currlength = 0;
      $allcontigs++;
    }
    else {$currlength += length($line);}
  }
  $alllength += $currlength;
  $meanalllength = $alllength / $allcontigs;
  push(@rallcontigs, $allcontigs);
  push(@ralllengths, $alllength);
  push(@rmeanalllengths, $meanalllength);
  close(SEQ);
}
open (RSC, "> $rscript");
print RSC "rallcontigs <- c($rallcontigs[0]";
for ($i = 1; $i < scalar(@rallcontigs); $i++) {print RSC ", $rallcontigs[$i]";}
print RSC ")\nralllengths <- c($ralllengths[0]";
for ($i = 1; $i < scalar(@ralllengths); $i++) {print RSC ", $ralllengths[$i]";}
print RSC ")\nrmeanalllengths <- c($rmeanalllengths[0]";
for ($i = 1; $i < scalar(@rmeanalllengths); $i++) {print RSC ", $rmeanalllengths[$i]";}
print RSC ")\nrmatchingcontigs <- c($rmatchingcontigs[0]";
for ($i = 1; $i < scalar(@rmatchingcontigs); $i++) {print RSC ", $rmatchingcontigs[$i]";}
print RSC ")\nrmatchinglengths <- c($rmatchinglengths[0]";
for ($i = 1; $i < scalar(@rmatchinglengths); $i++) {print RSC ", $rmatchinglengths[$i]";}
print RSC ")\nrmeanmatchinglengths <- c($rmeanmatchinglengths[0]";
for ($i = 1; $i < scalar(@rmeanmatchinglengths); $i++) {print RSC ", $rmeanmatchinglengths[$i]";}
print RSC ")\nrmatchinglongests <- c($rmatchinglongests[0]";
for ($i = 1; $i < scalar(@rmatchinglongests); $i++) {print RSC ", $rmatchinglongests[$i]";}
print RSC ")\nx <-c(0";
for ($i = 1; $i < scalar(@rallcontigs); $i++) {print RSC ", $i";}
if ($pdffile =~ m/pdf$/) {print RSC ")\npdf(\"$pdffile\", height = $pdfheight, width = $pdfwidth)\n";}
elsif ($pdffile =~ m/svg$/) {print RSC ")\nsvg(\"$pdffile\", height = $pdfheight, width = $pdfwidth)\n";}
else {print RSC ")\npng(\"$pdffile\", height = $pdfheight, width = $pdfwidth, units = \"in\")\n";}
#print RSC "plot\(x, means, ylim = c\($lbound, $ubound\), type = \"n\", xlab = \"Rounds of Polishing\", ylab = \"Fraction Identical to Reference\"\)\n";
print RSC "par(mar = c(5, 5, 3, 5))\n";
$ubound = $ralllengths[0];
for ($i = 1; $i < scalar(@ralllengths); $i++) {
  if ($ralllengths[$i] > $ubound) {$ubound = $ralllengths[$i];}
}
$ubound = int(1.05 * $ubound);
print RSC "plot(x, rmatchinglengths, ylim = c(0, $ubound), type = \"l\", xlab = \"Cycle\", ylab = \"Length\", lty = 1, col = \"blue\")\n";
print RSC "lines(x, rmatchinglongests, lty = 4, col = \"blue\")\n";
print RSC "lines(x, rmeanmatchinglengths, lty = 2, col = \"blue\")\n";
print RSC "lines(x, ralllengths, lty = 1, col = \"black\")\n";
print RSC "lines(x, rmeanalllengths, lty = 2, col = \"black\")\n";
print RSC "par(new = TRUE)\n";
$lowbound = $rallcontigs[0]; $highbound = $rallcontigs[0];
for ($i = 1; $i < scalar(@rallcontigs); $i++) {
  if ($rallcontigs[$i] < $lowbound) {$lowbound = $rallcontigs[$i];}
  if ($rallcontigs[$i] > $highbound) {$highbound = $rallcontigs[$i];}
}
for ($i = 0; $i < scalar(@rmatchingcontigs); $i++) {
  if ($rmatchingcontigs[$i] < $lowbound) {$lowbound = $rmatchingcontigs[$i];}
  if ($rmatchingcontigs[$i] > $highbound) {$highbound = $rmatchingcontigs[$i];} #This should never happen.
}
if ($highbound > 5) {
  $lowbound = int(0.7 * $lowbound);
  $highbound = int(1.3 * $highbound);
}
else {$lowbound = 0; $highbound += 2;}
print RSC "plot(x, rallcontigs, ylim = c($lowbound, $highbound), type = \"l\", xlab = \"\", ylab = \"\", xaxt = \"n\", yaxt = \"n\", lty = 3, col = \"black\")\n";
print RSC "lines(x, rmatchingcontigs, lty = 3, col = \"blue\")\n";
print RSC "axis(side = 4)\nmtext(\"Count\", side = 4, line = 3)\ndev.off()\n";
close (RSC);
($shellname = $rscript) =~ s/R$/sh/;
open (SHL, "> $shellname");
print SHL "source /etc/profile\nmodule load bioinfo\nmodule load R\nRscript $rscript\n";
close (SHL);
$shellname = "./".$shellname;
system("chmod 755 $shellname");
system("$shellname");
