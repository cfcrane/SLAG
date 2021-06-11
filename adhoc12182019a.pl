#!/usr/bin/env perl
use warnings;
#This script was the first to build a section of the SLAG results table.
open (LIS, "< $ARGV[0]"); #filelist12182019a.txt
open (OUT, "> $ARGV[1]");
$particle = $ARGV[2]; #Zt
$assembler = $ARGV[3]; #unicycler
$readlength = $ARGV[4]; #7-14kb, no spaces
while ($file = <LIS>) {
  chomp $file;
  open (INP, "< $file");
  while ($line = <INP>) {
    chomp $line;
    @vars = split(/\s+/, $line);
    if ($line =~ m/longest matching target is/) {
      $vars[0] =~ s/$particle//;
      $currid = $vars[0];
    }
    if ($line =~ m/Number of matching contigs/) {$numbermatching{$currid} = $vars[5]; $numberall{$currid} = $vars[11];}
    elsif ($line =~ m/Mean longest matching contig length/) {$meanlongestmatching{$currid} = $vars[-1];}
    elsif ($line =~ m/Mean longest any contig length/) {$meanlongestany{$currid} = $vars[-1];}
  }
  close(INP);
}
$extmeanlongestmatching = 0;
$N = 0;
for $key (sort(keys(%numbermatching))) {
  $extmeanlongestmatching += $meanlongestmatching{$key}; $N++;
  print OUT "$key\t$assembler\t$readlength\t$numbermatching{$key}\t$meanlongestmatching{$key}\t$numberall{$key}";
  print OUT "\t$meanlongestany{$key}\n";
}
$extmeanlongestmatching /= $N;
$stdevlongestmatching = 0;
for $key (sort(keys(%numbermatching))) {$stdevlongestmatching += ($meanlongestmatching{$key} - $extmeanlongestmatching)**2;}
$stdevlongestmatching = sqrt($stdevlongestmatching / ($N - 1));
$sterrlongestmatching = $stdevlongestmatching / sqrt($N);
print OUT "\nOverall mean = $extmeanlongestmatching stddev = $stdevlongestmatching stderr = $sterrlongestmatching\n";

#Number of matching contigs = 1 number of all contigs = 1 number of identities = 1
#ANQ91929.1Zt longest matching target is 147625
#Mean longest matching contig length for ANQ91929.1_unicycler_increment_intact is 147625
#Standard deviation of longest matching contig length for ANQ91929.1_unicycler_increment_intact is undef
#Standard error of the mean for longest matching contig length for ANQ91929.1_unicycler_increment_intact is undef
#Median longest matching contig length for ANQ91929.1_unicycler_increment_intact is mean
#Range of longest matching contig lengths is 147625 to 147625
#Mean longest any contig length for ANQ91929.1_unicycler_increment_intact is 147625
#Standard deviation of longest any contig length for ANQ91929.1_unicycler_increment_intact is undef
#Standard error of the mean for longest any contig length for ANQ91929.1_unicycler_increment_intact is undef
#Median longest any-contig length for ANQ91929.1_unicycler_increment_intact is mean
#Range of longest any contig length for ANQ91929.1_unicycler_increment_intact is 147625 to 147625
#Mean highest identity for ANQ91929.1_unicycler_increment_intact is 74.286
#Median highest identity for ANQ91929.1_unicycler_increment_intact is mean
#Range of highest identities is 74.286 to 74.286
#There are 21 .subset.contigs files, 21 parent .contigs files, and 4 blast vs retrieved files.
#diag: totncontigs = 75 nkeys = 21 totaltncontigs = 21 nsisters = 21
#Mean number of contigs produced = 1 mean that match the target = 3.57142857142857
#Longest matching contig lengths ordered by target identifier for comparison with other assemblies:
#ANQ91929.1Zt    147625
