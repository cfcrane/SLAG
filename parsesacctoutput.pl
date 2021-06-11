#!/usr/bin/env perl
use warnings;
#This script parses sacct output formatted as --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,MaxPages,TotalCPU,nnodes,ncpus,nodelist
#This script finds three sacct lines for the last job of the searched name to be completed, unless the job was cancelled.
open (LIS, "< $ARGV[0]"); #searchstrings03082021a.txt
open (INP, "< $ARGV[1]"); #attemptatjobstatsfrom01Feb.txt
open (OUT, "> $ARGV[2]"); 
while ($line = <LIS>) {
  chomp $line;
  @duo = split(/\s+/, $line, 2);
  push (@searchstringsa, $duo[0]);
  push (@searchstringsb, $duo[1]);
}
while ($line = <INP>) {
  chomp $line;
  #   ycrane 10403665     rsZmisosp+    brown-a     FAILED 4-04:00:00 2021-03-08T09:25:19 2021-03-08T20:04:54   10:39:35                              1         24      brown-a109
  #             10403665.ba+      batch                FAILED            2021-03-08T09:25:19 2021-03-08T20:04:54   10:39:35  13558984K  57792524K        1         24      brown-a109
  #                       10403665.ex+     extern             COMPLETED            2021-03-08T09:25:19 2021-03-08T20:04:56   10:39:37       928K    140376K        1         24      brown-a109
  #   batch                FAILED            2021-03-08T09:25:19 2021-03-08T20:04:54   10:39:35  13558984K  57792524K        1         24      brown-a109
  $found = 0;
  for ($i = 0; $i < scalar(@searchstringsa); $i++) {
    if ($line =~ m/$searchstringsa[$i]/) {$found = 1; $savedname = $searchstringsa[$i];}
  }
  for ($i = 0; $i < scalar(@searchstringsb); $i++) {
    if ($line =~ m/$searchstringsb[$i]/) {$found = 1; $savedname = $searchstringsb[$i];}
  }
  if ($found == 1) {
    @vars = split(/\s+/, $line);
    $jobid{$savedname} = $vars[2];
  }
}
for $key (keys(%jobid)) {$master{$jobid{$key}} = $key;}
for $key (sort(keys(%jobid))) {print OUT "$key\t$jobid{$key}\n";}
seek (INP, 0, 0);
while ($line = <INP>) {
  chomp $line;
  @vars = split(/\s+/, $line);
  ($parta, $partb) = split(/\./, $vars[1], 2);
  if (exists($master{$parta})) {
    #print OUT "diag: vars[1] = $vars[1] parta = $parta partb = $partb master{$parta} = $master{$parta}\n";
    if ($vars[2] eq "batch") {
      #Get memory usage as vars[7]
      $maxrmem{$master{$parta}} = $vars[7];
      $maxvmem{$master{$parta}} = $vars[8];
      $maxpages{$master{$parta}} = $vars[9];
    }
    if ($vars[2] eq "extern") {
      #Get runtime as vars[4] and convert to seconds.
      $runtime{$master{$parta}} = $vars[6];
    }
  }
}
for $key (sort(keys(%maxrmem))) {
  $maxrmem{$key} =~ s/M/000000/;
  $maxrmem{$key} =~ s/K/000/;
  $maxvmem{$key} =~ s/M/000000/;
  $maxvmem{$key} =~ s/K/000/;
  $days = 0;
  if ($runtime{$key} =~ m/\-/) {
    ($days, $filler) = split(/\-/, $runtime{$key});
    ($hours, $minutes, $seconds) = split(/:/, $filler);
  }
  else {($hours, $minutes, $seconds) = split(/:/, $runtime{$key});}
  $time = $days * 24 * 3600 + $hours * 3600 + $minutes * 60 + $seconds;
  $runtime{$key} = $time;
}
for ($i = 0; $i < scalar(@searchstringsa); $i++) {
  print OUT "$searchstringsa[$i]\t$runtime{$searchstringsa[$i]}\t$maxrmem{$searchstringsa[$i]}\t$maxvmem{$searchstringsa[$i]}\t$maxpages{$searchstringsa[$i]}\t$searchstringsb[$i]\t$runtime{$searchstringsb[$i]}\t$maxrmem{$searchstringsb[$i]}\t$maxvmem{$searchstringsb[$i]}\t$maxpages{$searchstringsb[$i]}\n";
}
