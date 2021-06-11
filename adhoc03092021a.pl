#!/usr/bin/env perl
use warnings;
#This script collected contig counts and lengths from files in a list.
open (LIS, "< $ARGV[0]"); #filelist03092021a.txt
open (OUT, "> $ARGV[1]");
while ($pile = <LIS>) {
  chomp $pile;
  ($file, $name) = split(/\t/, $pile);
  if ($file =~ m/atram/) {push(@atramfiles, $file); push(@atramnames, $name);}
  else {push(@slagfiles, $file); push (@slagnames, $name);}
}
for ($i = 0; $i < scalar(@atramfiles); $i++) {
  open (INP, "< $atramfiles[$i]");
  $cyclescompleted{$atramnames[$i]} = "NULL";
  while ($line = <INP>) {
    chomp $line;
    if ($line =~ m/\t/) {
      @vars = split(/\t/, $line);
      #/scratch/brown/ycrane/atramout/Zeaisocitratedehydrogenasesbyatram2.ERR328821xcleanedforatram2_Zeaisocitratedehydrogenases.all_contigs.fasta     198     995     6886    216985
      $atramcounts{$atramnames[$i]} = $vars[1];
      $atrammedians{$atramnames[$i]} = $vars[2];
      $atramlengths{$atramnames[$i]} = $vars[3];
      $atrambases{$atramnames[$i]} = $vars[4];
    }
  }
}
for ($i = 0; $i < scalar(@slagfiles); $i++) {
  open (INP, "< $slagfiles[$i]");
  #Number of cycles completed = 21
  #Longest contig length = 6430 for contig NODE_1_length_6430_cov_8.052101 in file Zeaisocitratedehydrogenasesbm_assembly19.subset.contigs
  #Shortest median contig length = 1095.5 longest median length = 1735
  #Fewest contigs = 15 most contigs = 24
  $mintotalbases = 2**30;
  $maxtotalbases = 0;
  while ($line = <INP>) {
    chomp $line;
    @vars = split(/\s+/, $line); 
    if ($line =~ m/cycles/) {$cyclescompleted{$slagnames[$i]} = $vars[5];}
    elsif ($line =~ m/Longest contig/) {$longestcontig{$slagnames[$i]} = $vars[4];}
    elsif ($line =~ m/median contig/) {$shortestmedian{$slagnames[$i]} = $vars[5]; $longestmedian{$slagnames[$i]} = $vars[10];}
    elsif ($line =~ m/most contigs/) {$fewestcontigs{$slagnames[$i]} = $vars[3]; $mostcontigs{$slagnames[$i]} = $vars[7];}
    else {
      if ($vars[-1] < $mintotalbases) {$mintotalbases = $vars[-1];}
      if ($vars[-1] > $maxtotalbases) {$maxtotalbases = $vars[-1];}
    }
  }
  $minbases{$slagnames[$i]} = $mintotalbases;
  $maxbases{$slagnames[$i]} = $maxtotalbases;
}
for ($i = 0; $i < scalar(@atramnames); $i++) {
  print OUT "$atramnames[$i]\tCC";
  print OUT "\t$atrammedians{$atramnames[$i]}\t$atramlengths{$atramnames[$i]}\t$atrambases{$atramnames[$i]}";
  print OUT "\t$cyclescompleted{$atramnames[$i]}\t$shortestmedian{$atramnames[$i]}-$longestmedian{$atramnames[$i]}";
  print OUT "\t$longestcontig{$atramnames[$i]}\t$fewestcontigs{$atramnames[$i]}-$mostcontigs{$atramnames[$i]}";
  print OUT "\t$minbases{$atramnames[$i]}-$maxbases{$atramnames[$i]}\n";
}
