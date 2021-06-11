#!/usr/bin/env perl
use warnings;
#This script collected run times per cycyle for contig files.
open (INP, "< $ARGV[0]"); #cycletimesfromlogout0517.txt
open (OUT, "> $ARGV[1]");
while ($line = <INP>) {
  chomp $line;
  #Processing slurm-10426128.out with code aTRAMfullStanley leading to atramout/Stanleyhistonedeacetylasesbyatram2.SRR9125476_htsclepurrenrefforatram2_Zeahistonedeacetylases.all_contigs.fasta
  #Run ended at 933785 with 8 cycles completed in 2466 seconds (308.25 sec/cycle)
  if ($line =~ m/Processing/) {
    @vars = split(/\s+/, $line);
    $file = $vars[-1];
  }
  elsif ($line =~ m/Run ended at/) {
    @vars = split(/\s+/, $line);
    $vars[-2] =~ s/\(//;
    $secpercycle{$file} = $vars[-2];
  }
}
for $key (sort(keys(%secpercycle))) {print OUT "$key\t$secpercycle{$key}\n";}
