#!/usr/bin/env perl
use warnings;
#This script calculated correlation coefficients of aTRAM2 runtime to SLAG runtime.
open (INP, "< $ARGV[0]"); #cycletimestable0531.txt
$outfile = $ARGV[1];
$rscript = $ARGV[2];
$soughtreadset = $ARGV[3]; #Zea, Stanley, hemiStanley, must equal entry in INP.
while ($line = <INP>) {
  chomp $line;
  #Zea     transaminase    2055.24 122.00  16.846
  @vars = split(/\t/, $line);
  if ($vars[0] eq $soughtreadset) {
    push(@x, $vars[2]);
    push(@y, $vars[3]);
  }
}
open (RSC, "> $rscript");
print RSC "x <- c($x[0]";
for ($i = 1; $i < scalar(@x); $i++) {print RSC ", $x[$i]";}
print RSC ")\ny <- c($y[0]";
for ($i = 1; $i < scalar(@y); $i++) {print RSC ", $y[$i]";}
print RSC ")\ncapture.output(cor.test(x, y, alternative=\"two.sided\", method=\"pearson\"), file = \"$outfile\")\n";
close (RSC);
`Rscript $rscript`;
