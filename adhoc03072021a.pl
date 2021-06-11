#!/usr/bin/env perl
use warnings;
#This script reformatted fasta headers so that localassembly1115.pl can use the resulting database.
open (SEQ, "< $ARGV[0]");
open (OUT, "> $ARGV[1]");
while ($line = <SEQ>) {
  if ($line =~ m/>/) {
    chomp $line;
    #>SRR9125476.sra.17_R1 17 length=264
    $line =~ s/>//;
    @vars = split(/ /, $line);
    @tars = split(/\./, $vars[0]);
    $tars[2] =~ s/_R/./;
    print OUT ">$tars[0]", ".$tars[2] $vars[-1]\n";
  }
  else {print OUT $line;}
}
