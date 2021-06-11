#!/usr/bin/env perl
use warnings;
#This script examines contig lengths and accuracy over a set of cycles.
open (LIS, "< $ARGV[0]"); #ZeaferredoxinsSRcap3subsetcontigsfiles.txt, presumed to be in chronological order
open (OUT, "> $ARGV[1]"); #ZeaferredoxinsSRcap3subsetcontigseries.txt
$blastexe = $ARGV[2]; #/group/bioinfo/apps/apps/blast-2.10.0+/bin/blastn
$blastdb = $ARGV[3]; #/scratch/halstead/c/ccrane/GCF_000005005.2_B73_RefGen_v4_genomic_db
$blastoutstem = $ARGV[4]; #/scratch/halstead/c/ccrane/slagtests/ZeaferredoxinsSRcap3ssctgs
$rscript = $ARGV[5]; #ZeaferredoxinsSRcap3graphics.R
$pdffile = $ARGV[6]; #ZeaferredoxinsSRcap3graph.pdf
$pdfheight = $ARGV[7]; #3
$pdfwidth = $ARGV[8]; #4
$evalue = $ARGV[9]; #1e-10
@tars = split(/\//, $blastdb);
$blastdbkernel = $tars[-1];
$k = 0;
while ($file = <LIS>) { #LIS MUST be in chronological order of cycles.
  chomp $file;
  open (SEQ, "< $file");
  $ncontigs = 0;
  $totallength = 0;
  $greatestlength = 0;
  $currlength = 0;
  while ($line = <SEQ>) {
    chomp $line;
    if ($line =~ m/>/) {
      if ($currlength > $greatestlength) {$greatestlength = $currlength;}
      $totallength += $currlength;
      $currlength = 0;
      $ncontigs++;
    }
    else {$currlength += length($line);}
  }
  if ($currlength > $greatestlength) {$greatestlength = $currlength;}
  $totallength += $currlength;
  $meanlength = $totallength / $ncontigs;
  $blastoutfile = $blastoutstem."vs".$blastdbkernel."_".$k."_".$evalue.".txt";
  `$blastexe -query $file -db $blastdb -out $blastoutfile -evalue $evalue -dust no`;
  open (BLA, "< $blastoutfile");
  $numeratorsum = 0;
  $denominatorsum = 0;
  $nidentities = 0;
  while ($line = <BLA>) {
    if ($line =~ m/uery=/) {$nidentities++; $curridentities = 0;}
    if ($line =~ m/Identities =/) {
      chomp $line;
      $curridentities++;
      # Identities = 1704/1709 (99%), Gaps = 1/1709 (0%)
      if ($curridentities == 1) {
        @vars = split(/\s+/, $line);
        ($num, $denom) = split(/\//, $vars[3]);
        $numeratorsum += $num;
        $denominatorsum += $denom;
      }
    }
  }
  $meanidentities = 100 * $numeratorsum / $denominatorsum;
  print OUT "$ncontigs\t$totallength\t$meanlength\t$greatestlength\t$meanidentities\n";
  push(@rncontigs, $ncontigs);
  push(@rtotallengths, $totallength);
  push(@rmeanlengths, $meanlength);
  push(@rgreatestlengths, $greatestlength);
  push(@rmeanidentities, $meanidentities);
  $k++;
}
open (RSC, "> $rscript");
print RSC "rncontigs <- c($rncontigs[0]";
for ($i = 1; $i < scalar(@rncontigs); $i++) {print RSC ", $rncontigs[$i]";}
print RSC ")\nrtotallengths <- c($rtotallengths[0]";
for ($i = 1; $i < scalar(@rtotallengths); $i++) {print RSC ", $rtotallengths[$i]";}
print RSC ")\nrmeanlengths <- c($rmeanlengths[0]";
for ($i = 1; $i < scalar(@rmeanlengths); $i++) {print RSC ", $rmeanlengths[$i]";}
print RSC ")\nrgreatestlengths <- c($rgreatestlengths[0]";
for ($i = 1; $i < scalar(@rgreatestlengths); $i++) {print RSC ", $rgreatestlengths[$i]";}
print RSC ")\nrmeanidentities <- c($rmeanidentities[0]";
for ($i = 1; $i < scalar(@rmeanidentities); $i++) {print RSC ", $rmeanidentities[$i]";}
print RSC ")\nx <-c(0";
for ($i = 1; $i < scalar(@rncontigs); $i++) {print RSC ", $i";}
print RSC ")\npdf(\"$pdffile\", height = $pdfheight, width = $pdfwidth)\n";
#print RSC "plot\(x, means, ylim = c\($lbound, $ubound\), type = \"n\", xlab = \"Rounds of Polishing\", ylab = \"Fraction Identical to Reference\"\)\n";
print RSC "par(mar = c(5, 5, 3, 5))\n";
$ubound = $rgreatestlengths[0];
for ($i = 1; $i < scalar(@rgreatestlengths); $i++) {
  if ($rgreatestlengths[$i] > $ubound) {$ubound = $rgreatestlengths[$i];}
}
$ubound = int(1.05 * $ubound);
print RSC "plot(x, rmeanlengths, ylim = c(0, $ubound), type = \"l\", xlab = \"Cycle\", ylab = \"Length\", lty = 2)\n";
print RSC "lines(x, rgreatestlengths, lty = 1)\npar(new = TRUE)\n";
$lowbound = $rncontigs[0]; $highbound = $rncontigs[0];
for ($i = 1; $i < scalar(@rncontigs); $i++) {
  if ($rncontigs[$i] < $lowbound) {$lowbound = $rncontigs[$i];}
  if ($rncontigs[$i] > $highbound) {$highbound = $rncontigs[$i];}
}
if ($highbound > 5) {
  $lowbound = int(0.7 * $lowbound);
  $highbound = int(1.3 * $highbound);
}
else {$lowbound = 0; $highbound += 2;}
print RSC "plot(x, rncontigs, ylim = c($lowbound, $highbound), type = \"l\", xlab = \"\", ylab = \"\", xaxt = \"n\", yaxt = \"n\", lty = 3)\n";
print RSC "axis(side = 4)\nmtext(\"Count\", side = 4, line = 3)\ndev.off()\n";
$pdffileb = "sec".$pdffile;
print RSC "pdf(\"$pdffileb\", height = $pdfheight, width = $pdfwidth)\n";
print RSC "par(mar = c(5, 5, 3, 5))\n";
$ubound = $rtotallengths[0];
for ($i = 1; $i < scalar(@rtotallengths); $i++) {
  if ($rtotallengths[$i] > $ubound) {$ubound = $rtotallengths[$i];}
}
$ubound = int(1.05 * $ubound);
print RSC "plot(x, rtotallengths, ylim = c(0, $ubound), type = \"l\", xlab = \"Cycle\", ylab = \"Total contig length\", lty = 1)\n";
$highbound = $rmeanidentities[0];
$lowbound = $rmeanidentities[0];
for ($i = 1; $i < scalar(@rmeanidentities); $i++) {
  if ($rmeanidentities[$i] < $lowbound) {$lowbound = $rmeanidentities[$i];}
  if ($rmeanidentities[$i] > $highbound) {$highbound = $rmeanidentities[$i];}
}
$lowbound = int($lowbound - 2); $highbound = int($highbound + 2);
print RSC "par(new = TRUE)\n";
print RSC "plot(x, rmeanidentities, ylim = c($lowbound, $highbound), type = \"l\", xlab = \"\", ylab = \"\", xaxt = \"n\", yaxt = \"n\", lty = 2)\n";
print RSC "axis(side = 4)\nmtext(\"Percent identity\", side = 4, line = 3)\ndev.off()\n";
close (RSC);
($shellname = $rscript) =~ s/R$/sh/;
open (SHL, "> $shellname");
print SHL "source /etc/profile\nmodule load bioinfo\nmodule load R\nRscript $rscript\n";
close (SHL);
$shellname = "./".$shellname;
system("$shellname");
