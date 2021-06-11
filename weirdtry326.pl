#!/usr/bin/env perl
use warnings;
#This script attempts to pair fasta and slurm files.
open (LOL, "< $ARGV[0]"); #listofcontigfilelists.txt
open (SLL, "< $ARGV[1]"); #slurmoutfiles03262021.txt
open (ATT, "< $ARGV[2]"); #attemptatjobstatsfrom01Febto26Mar.txt
open (OUT, "> $ARGV[3]"); #
$user = $ARGV[4];
$tol = $ARGV[4]; #3
while ($listline = <LOL>) {
  chomp $listline;
  ($listfile, $code) = split(/\s+/, $listline);
  push (@{$listfiles{$code}}, $listfile);
  open (LIS, "< $listfile");
#  print OUT "list file is $listfile\n";
  $latesttime = 0;
  $linecount = 0;
  while ($file = <LIS>) {
    chomp $file;
    @stats = stat($file);
    if ($stats[9] > $latesttime) {$latesttime = $stats[9]; $savedsize = $stats[7]; $savedfile = $file;}
    $linecount++;
  }
  close(LIS);
  push (@savedfiles, $savedfile);
  $savedtimes{$savedfile} = $latesttime; #in seconds since the Epoch
  $savedsizes{$savedfile} = $savedsize;
#  print OUT "$linecount lines in $listfile\n";
  $codes{$savedfile} = $code;
  push(@{$children{$code}}, $savedfile);
}
for $key (keys(%codes)) {print OUT "$key code = $codes{$key}\n";}
#for ($i = 0; $i < scalar(@savedfiles); $i++) {
#  print OUT "saved times, sizes, files, i = $savedtimes{$savedfiles[$i]}\t$savedsizes{$savedfiles[$i]}\t$savedfiles[$i]\t$i\n";
#}
while ($slurmline = <SLL>) {
  chomp $slurmline;
  @sars = split(/\s+/, $slurmline);
  $slurmfile = $sars[-1];
  @stats = stat($slurmfile);
  $slurmtimes{$slurmfile} = $stats[9];
}
#for $key (sort(keys(%slurmtimes))) {print OUT "slurm time = $slurmtimes{$key}\t$key\n";}
for ($i = 0; $i < scalar(@savedfiles); $i++) {
  for $key (keys(%slurmtimes)) {
    if (abs($slurmtimes{$key} - $savedtimes{$savedfiles[$i]}) < $tol) {
      $mates{$savedfiles[$i]} = $key;
      ($jobnumber = $key) =~ s/slurm-//;
      $jobnumber =~ s/.out$//;
      $jobnumbers{$savedfiles[$i]} = $jobnumber;
      $filesbyjob{$jobnumber} = $savedfiles[$i];
#      print OUT "line 48: $savedtimes{$savedfiles[$i]} $slurmtimes{$key} $savedfiles[$i] $key $jobnumber\n";
    }
  }
}
for ($i = 0; $i < scalar(@savedfiles); $i++) {
#  if ($savedfiles[$i] =~ m/atramtemp/) {
#    print OUT "line 54: savedfiles[$i] = $savedfiles[$i] saved time = $savedtimes{$savedfiles[$i]} mate = $mates{$savedfiles[$i]}\n";
#  }
  if (!exists($mates{$savedfiles[$i]})) {
#    print OUT "line 57, $savedfiles[$i] has no mate.\n";
    $mindist = 2**30;
    for $index (keys(%slurmtimes)) {
      if (abs($savedtimes{$savedfiles[$i]} - $slurmtimes{$index}) < $mindist) {
        $mindist = abs($savedtimes{$savedfiles[$i]} - $slurmtimes{$index});
#        print OUT "DIAG: index = $index savedtimes{$savedfiles[$i]} = $savedtimes{$savedfiles[$i]} slurmtimes{$index} = $slurmtimes{$index} mindist is $mindist\n";
        $savedindex = $index;
      }
    }
    $mates{$savedfiles[$i]} = $savedindex;
    ($jobnumber = $savedindex) =~ s/slurm-//;
    $jobnumber =~ s/.out$//;
    $jobnumbers{$savedfiles[$i]} = $jobnumber;
    $filesbyjob{$jobnumber} = $savedfiles[$i];
 
#    print OUT "$savedfiles[$i] with $savedtimes{$savedfiles[$i]} is closest to $savedindex with $slurmtimes{$savedindex}\n";
  }
}
for $key (sort(keys(%mates))) {
  open (SLM, "< $mates{$key}");
  $lastcycle = -1;
  $ncycles{$key} = 0;
  while ($line = <SLM>) {
    chomp $line;
    if ($line =~ m/Creating/) {
      @vars = split(/\s+/, $line);
      $lastcycle = $vars[-1];
    }
    elsif ($line =~ m/retrieval/) {
      @vars = split(/\s+/, $line); 
      ($lastcycle = $vars[-1]) =~ s/\.//;
      $lastcycle++; #Because slag cycles begin with 0.
    }
    $ncycles{$key} = $lastcycle;
  }
}
for $key (sort {$a <=> $b} (keys(%children))) {
  #The files tested here all produced at least one contig. First-round assembly failures are not counted.
  $meancyclestried = 0;
  $n = 0;
  for $i (0 .. $#{$children{$key}}) {
#    print OUT "diag: i = $i cycles tried = $ncycles{$children{$key}[$i]} in file $children{$key}[$i]\n";
    $meancyclestried += $ncycles{$children{$key}[$i]};
    $n++;
  }
  $totalcyclestried{$key} = $meancyclestried;
  $meancyclestried /= $n;
  print OUT "Group $key mean cycles tried = $meancyclestried over $n files. Total cycles tried = $totalcyclestried{$key}\n";
  $totaleltime{$key} = 0; $totalcotime{$key} = 0; $totalresmem{$key} = 0; $totalvirmem{$key} = 0; $totalpfaults{$key} = 0;
  $nobs{$key} = 0;
}
while ($line = <ATT>) {
  #Get elapsed and CPU time, resident and virtual memory, and page faults.
  #filesbyjob is the relevant hash.
  #ycrane 10424603     raTahexsp+    brown-a  COMPLETED 6-06:00:00 2021-03-10T18:35:02 2021-03-10T19:23:09   00:48:07                                  03:27:41        1         24      brown-a256
  #          10424603.ba+      batch             COMPLETED            2021-03-10T18:35:02 2021-03-10T19:23:09   00:48:07   2500656K   9337184K      250   03:27:41        1         24      brown-a256
  #                    10424603.ex+     extern             COMPLETED            2021-03-10T18:35:02 2021-03-10T19:23:14   00:48:12       928K    140376K        0  00:00.002        1         24      brown-a256
  if ($line =~ m/$user/) {
    chomp $line; @vars = split(/\s+/, $line); 
    if (exists($filesbyjob{$vars[2]})) {
      $parent = $filesbyjob{$vars[2]};
      $pine = <ATT>; chomp $pine;
      @pars = split(/\s+/, $pine);
      $eltime = $pars[6];
      $cotime = $pars[10];
      $resmem = $pars[7];
      $virmem = $pars[8];
      $pfault = $pars[9];
      $nobs{$codes{$parent}}++;
      if ($eltime =~ m/\-/) {($days, $time) = split(/-/, $eltime);}
      else {$days = 0; $time = $eltime;}
      ($hours, $minutes, $seconds) = split(/:/, $time);
      $neweltime = 3600 * (24 * $days + $hours) + 60 * $minutes + $seconds;
      push (@{$eltimes{$codes{$parent}}}, $neweltime);
      $totaleltime{$codes{$parent}} += $neweltime;
      if ($cotime =~ m/\-/) {($days, $time) = split(/-/, $cotime);}
      else {$days = 0; $time = $cotime;}
      ($hours, $minutes, $seconds) = split(/:/, $time);
      $newcotime = 3600 * (24 * $days + $hours) + 60 * $minutes + $seconds;
      push(@{$cotimes{$codes{$parent}}}, $newcotime);
      $totalcotime{$codes{$parent}} += $newcotime;
      if ($resmem =~ m/K/) {$resmem =~ s/K//; $resmem *= 1000;}
      elsif ($resmem =~ m/M/) {$resmem =~ s/M//; $resmem *= 1000000;}
      push(@{$resmems{$codes{$parent}}}, $resmem);
      $totalresmem{$codes{$parent}} += $resmem;
      if ($virmem =~ m/K/) {$virmem =~ s/K//; $virmem *= 1000;}
      elsif ($virmem =~ m/M/) {$virmem =~ s/M//; $virmem *= 1000000;}
      push(@{$virmems{$codes{$parent}}}, $virmem);
      $totalvirmem{$codes{$parent}} += $virmem;
      push(@{$pfaults{$codes{$parent}}}, $pfault);
      $totalpfault{$codes{$parent}} += $pfault;
      print OUT "job = $vars[2] eltime = $neweltime cotime = $newcotime resmem = $resmem virmem = $virmem pfault = $pfault code = $codes{$parent}\n";
    }
  }
}
for $key (sort(keys(%nobs))) {
  $meaneltimepercycle{$key} = $totaleltime{$key} / $totalcyclestried{$key};
  $meancotimepercycle{$key} = $totalcotime{$key} / $totalcyclestried{$key};
  $meanpfaultpercycle{$key} = $totalpfault{$key} / $totalcyclestried{$key};
  $meaneltime{$key} = $totaleltime{$key} / $nobs{$key};
  $meancotime{$key} = $totalcotime{$key} / $nobs{$key};
  $meanresmem{$key} = $totalresmem{$key} / $nobs{$key};
  $meanvirmem{$key} = $totalvirmem{$key} / $nobs{$key};
  $meanpfault{$key} = $totalpfault{$key} / $nobs{$key};
  print OUT "$key $nobs{$key} $meaneltime{$key} $meancotime{$key} $meanresmem{$key} $meanvirmem{$key} $meanpfault{$key}\n";
  $seeltimepercycle{$key} = 0;
  $secotimepercycle{$key} = 0;
  $sepfpercycle{$key} = 0;
  $seeltime{$key} = 0;
  $secotime{$key} = 0;
  $seresmem{$key} = 0;
  $sevirmem{$key} = 0;
  $sepfault{$key} = 0;
}
for $key (sort(keys(%nobs))) {
  for $i (0 .. $#{$eltimes{$key}}) {$seeltime{$key} += ($eltimes{$key}[$i] - $meaneltime{$key}) ** 2;}
  for $i (0 .. $#{$cotimes{$key}}) {$secotime{$key} += ($cotimes{$key}[$i] - $meancotime{$key}) ** 2;}
  for $i (0 .. $#{$resmems{$key}}) {$seresmem{$key} += ($resmems{$key}[$i] - $meanresmem{$key}) ** 2;}
  for $i (0 .. $#{$virmems{$key}}) {$sevirmem{$key} += ($virmems{$key}[$i] - $meanvirmem{$key}) ** 2;}
  for $i (0 .. $#{$pfaults{$key}}) {$sepfault{$key} += ($pfaults{$key}[$i] - $meanpfault{$key}) ** 2;}
  if ($nobs{$key} > 1) {
    $seeltime{$key} /= ($nobs{$key} - 1);
    $secotime{$key} /= ($nobs{$key} - 1);
    $seresmem{$key} /= ($nobs{$key} - 1);
    $sevirmem{$key} /= ($nobs{$key} - 1);
    $sepfault{$key} /= ($nobs{$key} - 1);
    $seeltime{$key} = sqrt($seeltime{$key}) / sqrt($nobs{$key});
    $secotime{$key} = sqrt($secotime{$key}) / sqrt($nobs{$key});
    $seresmem{$key} = sqrt($seresmem{$key}) / sqrt($nobs{$key});
    $sevirmem{$key} = sqrt($sevirmem{$key}) / sqrt($nobs{$key});
    $sepfault{$key} = sqrt($sepfault{$key}) / sqrt($nobs{$key});
    $gbmeanresmem{$key} = $meanresmem{$key} / 1000000000;
    $gbmeanvirmem{$key} = $meanvirmem{$key} / 1000000000;
    $gbseresmem{$key} = $seresmem{$key} / 1000000000;
    $gbsevirmem{$key} = $sevirmem{$key} / 1000000000;
  }
  else {
    $seeltime{$key} = "-"; $secotime{$key} = "-"; $seresmem{$key} = "-"; $sevirmem{$key} = "-"; $sepfault{$key} = "-";
    $gbseresmem{$key} = "-"; $gbsevirmem{$key} = "-";
  }
}
for $key (sort(keys(%nobs))) {
  print OUT "Group $key means and standard errors: elapsed time per run = $meaneltime{$key} +- $seeltime{$key}\n";
  print OUT "Group $key means and standard errors: CPU time per run = $meancotime{$key} +- $secotime{$key}\n";
  print OUT "Group $key means and standard errors: elapsed time per cycle = $meaneltimepercycle{$key} +- NA\n";
  print OUT "Group $key means and standard errors: CPU time per cycle = $meancotimepercycle{$key} +- NA\n";
  print OUT "Group $key means and standard errors: maximum resident size = $meanresmem{$key} +- $seresmem{$key}\n";
  print OUT "Group $key means and standard errors: maximum virtual size = $meanvirmem{$key} +- $sevirmem{$key}\n";
  print OUT "Group $key means and standard errors: maximum resident size in gb = $gbmeanresmem{$key} +- $gbseresmem{$key}\n";
  print OUT "Group $key means and standard errors: maximum virtual size in gb = $gbmeanvirmem{$key} +- $gbsevirmem{$key}\n";
  print OUT "Group $key means and standard errors: page faults per run = $meanpfault{$key} +- $sepfault{$key}\n";
  print OUT "Group $key means and standard errors: page faults per cycle = $meanpfaultpercycle{$key} +- NA\n";
}
