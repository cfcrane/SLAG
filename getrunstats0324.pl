#!/usr/bin/env perl
use warnings;
#This script sought runtime, number of completed cycles, maximum resident size, maximum virtual memory, and page faults for benchmarking jobs.  
open (INP, "< $ARGV[0]"); #attemptatjobstatsfrom01Febto24Mar.txt 
open (LOL, "< $ARGV[1]"); #listofcontigfilelists.txt 
open (JFN, "< $ARGV[2]"); #jobsfinished0326.txt
open (OUT, "> $ARGV[3]"); 
$user = $ARGV[4]; 
$firststring = $ARGV[5]; #"raTa raZm rsTa rsZm" 
$secondstring = $ARGV[6]; #"ca+ sp+"
$tol = $ARGV[7]; #5
$firststring =~ s/\"//g;
$secondstring =~ s/\"//g;
@firstcodes = split(/\s+/, $firststring);
@secondcodes = split(/\s+/, $secondstring);
while ($line = <JFN>) { #Code should be the same as in LOL.
  chomp $line;
  ($code, $count) = split(/\t/, $line);
  $jobsfinished{$code} = $count;
}
while ($line = <LOL>) {
  chomp $line;
  ($listfile, $listcode) = split(/ /, $line);
  push(@listfiles, $listfile);
  $parentlistcodes{$listfile} = $listcode;
  if (!defined($listcode)) {die "Listcode was not found in $ARGV[1].\n";}
  if (exists($listcodesums{$listcode})) {$listcodesums{$listcode}++;}
  else {$listcodesums{$listcode} = 1;}
}
for $key (sort(keys(%listcodesums))) {
  print OUT "listcodesums{$key} = $listcodesums{$key}\n";
  $cyclesums{$key} = 0;
  $eltimesums{$key} = 0;
  $cotimesums{$key} = 0;
  $resmemsums{$key} = 0;
  $virmemsums{$key} = 0;
  $pagefaultsums{$key} = 0;
}
$modtimes{"NULL"} = 0;
for ($k = 0; $k < scalar(@listfiles); $k++) {
  print OUT "listfiles[$k] = $listfiles[$k]\n";
  open (LIS, "< $listfiles[$k]") || die "$listfiles[$k] was not found.\n";
  while ($file = <LIS>) {
    chomp $file;
    if (-z $file) {print "$file is EMPTY.\n";}
    if ($file eq "NULL") {print "NULL file from $listfiles[$k], k = $k\n";}
    @vars = stat($file);
    $modtimes{$file} = $vars[9];
    $slurmouts{$file} = "NULL";
    $listcodes{$file} = $parentlistcodes{$listfiles[$k]};
    print OUT "file in $listfiles[$k] = $file last modified at $modtimes{$file} list code = $listcodes{$file}\n";
  }
  close(LIS);
}
while ($line = <INP>) {
  chomp $line;
#sacct output formatted as --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,MaxPages,TotalCPU,nnodes,ncpus,nodelist
#     ycrane 10403624     rsTahexca+    brown-a     FAILED 8-08:00:00 2021-03-07T20:19:46 2021-03-07T22:44:55   02:25:09                                  06:45:25        1         24      brown-a104
#          10403624.ba+      batch                FAILED            2021-03-07T20:19:46 2021-03-07T22:44:55   02:25:09  25261560K  63194636K       34   06:45:25        1         24      brown-a104 
#          10403624.ex+     extern             COMPLETED            2021-03-07T20:19:46 2021-03-07T22:44:58   02:25:12       920K    140376K        0  00:00.002        1         24      brown-a104 
  @vars = split(/\s+/, $line);
  if ($vars[1] eq $user) {
    #print OUT "diag: vars[3] = $vars[3]\n";
    $firstfound = 0; $secondfound = 0;
    if ($vars[3] eq "raTatrasp+") {print OUT "Code raTatrasp+ was seen in $line\n"}
    for ($i = 0; $i < scalar(@firstcodes); $i++) {
      if ($vars[3] =~ m/$firstcodes[$i]/) {$firstfound = 1;}
    }
    for ($i = 0; $i < scalar(@secondcodes); $i++) {
      if ($vars[3] =~ m/$secondcodes[$i]/) {$secondfound = 1;}
    }
    if ($firstfound == 1 && $secondfound == 1) {
      if ($vars[3] eq "raTatrasp+") {print OUT "Code raTatrasp+ passed line 72.\n";}
      if ($vars[5] !~ m/CANCELLED/) { #If WE did not cancel the job.  $pine can match CANCELLED if slurm cancelled the job.
        if ($vars[3] eq "raTatrasp+") {print OUT "Code raTatrasp+ passed line 74. Job number is $vars[2]\n";}
        $jobnumber = $vars[2];
        $currslurm = "slurm-".$jobnumber.".out";
        @stats = stat($currslurm);
        if ($vars[3] eq "raTatrasp+") {print OUT "Current slurm file is $currslurm and stats[9] = $stats[9] at line 79.\n";}
        #print OUT "line 53: currslurm = $currslurm stats[9] = $stats[9]\n";
        for $key (keys(%modtimes)) {
          if ($vars[3] eq "raTatrasp+") {print OUT "Value of key = $key for raTatrasp+ at line 82\n";}
          if (defined($modtimes{$key})) {
            if ($vars[3] eq "raTatrasp+") {print OUT "modtimes{$key} = $modtimes{$key} at line 84\n";}
            if (abs($stats[9] - $modtimes{$key}) < $tol) {
              #print OUT "line 56: key = $key modtime = $modtimes{$key} currslurm = $currslurm\n";
              $slurmouts{$key} = $currslurm;
              $jobnumbers{$key} = $jobnumber;
              #Get memory and time usage here.
              $pine = <INP>; chomp $pine;
              @pars = split(/\s+/, $pine);
              #          10403624.ba+      batch                FAILED            2021-03-07T20:19:46 2021-03-07T22:44:55   02:25:09  25261560K  63194636K       34   06:45:25        1         24      brown-a104 
              print OUT "$listcodes{$key} $key $pars[1] pars[6] = $pars[6] pars[7] = $pars[7] pars[8] = $pars[8] pars[9] = $pars[9] pars[10] = $pars[10]\n";
              if ($pars[6] =~ m/-/) {($days, $hourtime) = split(/-/, $pars[6]);}
              else {
                $days = 0;
                $hourtime = $pars[6];
              }
              ($hours, $minutes, $seconds) = split(/:/, $hourtime, 3);
              $eltime = 24 * 3600 * $days + 3600 * $hours + 60 * $minutes + $seconds;
              print OUT "diag: elapsed time for $key in $slurmouts{$key} is $eltime\n";
              $eltimesums{$listcodes{$key}} += $eltime;
              if ($pars[10] =~ m/-/) {($days, $hourtime) = split(/-/, $pars[10]);}
              else {
                $days = 0;
                $hourtime = $pars[10];
              }
              ($hours, $minutes, $seconds) = split(/:/, $hourtime, 3);
              $cotime = 24 * 3600 * $days + 3600 * $hours + 60 * $minutes + $seconds;
              $cotimesums{$listcodes{$key}} += $cotime;
              if ($pars[7] =~ m/K/) {
                $pars[7] =~ s/K//;
                $pars[7] *= 1000;
              }
              if ($ pars[7] =~ m/M/) {
                $pars[7] =~ s/M//;
                $pars[7] *= 1000000;
              }
              if ($pars[8] =~ m/K/) {
                $pars[8] =~ s/K//;
                $pars[8] *= 1000;
              }
              if ($ pars[8] =~ m/M/) {
                $pars[8] =~ s/M//;
                $pars[8] *= 1000000;
              }
              $resmemsums{$listcodes{$key}} += $pars[7];
              $virmemsums{$listcodes{$key}} += $pars[8];
              $pagefaultsums{$listcodes{$key}} += $pars[9];
            }
          }
        }
      }
    }
  }
}
for $key (sort(keys(%modtimes))) {
  if ($slurmouts{$key} ne "NULL") {
    print OUT "file = $key modtime = $modtimes{$key} job number = $jobnumbers{$key} slurm file = $slurmouts{$key}\n";
    if ($key =~ m/atram/) { #Get cycles completed from the slurm*out file.
      open (SLR, "< $slurmouts{$key}");
      $errorline = "";
      while ($line = <SLR>) {
        chomp $line;
        if ($line =~ m/Creating/) {$savedline = $line;}
        if ($line =~ m/ERROR/) {$errorline = $line;}
        if ($line =~ m/No new contigs/) {$errorline = $line;}
      }
      #2021-02-09 13:03:34 INFO: Creating new query files: iteration 21
      @vars = split(/\s+/, $savedline);
      $cyclescompleted{$key} = $vars[-1];
      $cyclesums{$listcodes{$key}} += $vars[-1];
      if (length($errorline) > 2) {
        if ($errorline =~ m/assembler failed/) {$errorcode{$key} = "assembly failed";}
        elsif ($errorline =~ m/TIME/) {$errorcode{$key} = "time expired";}
        elsif ($errorline =~ m/database is locked/) {$errorcode{$key} = "database locked";}
        elsif ($errorline =~ m/No new contig/) {$errorcode{$key} = "no new contigs";}
        else {$errorcode{$key} = "unspecified error in $slurmouts{$key}";}
      }
      else {$errorcode{$key} = "no error";}
    }
    elsif ($key =~ m/subset.contigs/) { #Get cycle number from the file name.
      #StanleyferredoxinsSRcap3extracted20.fasta.cap_20.subset.contigs
      @vars = split(/\./, $key);
      if ($vars[-3] =~ m/_/) {($filler, $number) = split(/_/, $vars[-3]);}
      else {print OUT "$key does not conform to the expected SLAG naming convention.\n"}
      while ($number =~ m/\D+/) {$number =~ s/\D+//g;}
      $number++; #Because SLAG starts with cycle 0.
      $cyclescompleted{$key} = $number;
      $errorcode{$key} = "not applicable";
      $cyclesums{$listcodes{$key}} += $number;
    }
    else {print OUT "$key does not conform to either convention.\n"; $cyclescompleted{$key} = "NULL";}
  }
}
for $key (sort(keys(%cyclescompleted))) {
  print OUT "$key\t$listcodes{$key}\tcycles completed = $cyclescompleted{$key}\terror code = $errorcode{$key}\n";
}
for $key (sort(keys(%cyclesums))) {print OUT "cycle sum for $key is $cyclesums{$key}\n";}
for $key (sort(keys(%eltimesums))) {print OUT "total elapsed time for $key is $eltimesums{$key}\n";}
for $key (sort(keys(%cotimesums))) {print OUT "total CPU time for $key is $cotimesums{$key}\n";}
for $key (sort(keys(%resmemsums))) {print OUT "total resident memory for $key is $resmemsums{$key}\n";}
for $key (sort(keys(%virmemsums))) {print OUT "total virtual memory for $key is $virmemsums{$key}\n";}
for $key (sort(keys(%pagefaultsums))) {print OUT "total page faults for $key are $pagefaultsums{$key}\n";}
for $key (sort(keys(%cyclesums))) {
  $eltimepercycle = $eltimesums{$key} / $cyclesums{$key};
  $cotimepercycle = $cotimesums{$key} / $cyclesums{$key};
  $coelratio = $cotimepercycle / $eltimepercycle;
  $pagefaultspercycle = $pagefaultsums{$key} / $cyclesums{$key};
  print OUT "Elapsed time per cycle for group $key = $eltimepercycle\nCPU time per cycle for group $key = $cotimepercycle\n";
  print OUT "CPU to elapsed time ratio for group $key = $coelratio\n";
  print OUT "Page faults per cycle for group $key = $pagefaultspercycle\n";
}
for $key (sort(keys(%eltimesums))) {
  $timeperrun = $eltimesums{$key} / $jobsfinished{$key};
  print OUT "Elapsed time per run for group $key = $timeperrun\n";
  $cputimeperrun = $cotimesums{$key} / $jobsfinished{$key};
  print OUT "CPU time per run for group $key = $cputimeperrun\n";
}
#        print OUT "$line\n";
#        $pine = <INP>; print OUT $pine;
#        chomp $pine;
#        @pars = split(/\s+/, $pine);
#        $maxresmem{$jobnumber} = $pars[7];
#        $maxvirmem{$jobnumber} = $pars[8];
#        $pagefaults{$jobnumber} = $pars[9];
#        $cputime{$jobnumber} = $pars[10];
#        $vine = <INP>; print OUT $vine;
#        chomp $vine;
#        @bars = split(/\s+/, $vine);
#        $walltime{$jobnumber} = $bars[6];
#        print OUT "$jobnumber $walltime{$jobnumber} $cputime{$jobnumber} $maxresmem{$jobnumber} $maxvirmem{$jobnumber} $pagefaults{$jobnumber}\n";
#        $slurmfile = "slurm-".$jobnumber.".out";
#        open (SLR, "< $slurmfile");
#        while ($dine = <SLR>) {
#          chomp $dine;
#          if ($dine =~ m/^atram.py /) {
#            $fine = <SLR>;
#            chomp $fine;
#            @fars = split(/\s+/, $fine);
#            $startmonth{$jobnumber} = $fars[1];
#            $startday{$jobnumber} = $fars[2];
#            $starttime{$jobnumber} = $fars[3];
#            $startyear{$jobnumber} = $fars[5];
#            print OUT "job = $jobnumber year = $startyear{$jobnumber} month = $startmonth{$jobnumber} day = $startday{$jobnumber} time = $starttime{$jobnumber}\n";
#          }
#          #Thu Mar 11 11:39:49 EST 2021
#          #2021-03-12 08:56:34 INFO: Blasting query against shards: iteration 12 //This is 
#          #2021-03-12 09:33:19 INFO: All 783 blast results completed
#          #2021-03-12 09:33:26 INFO: 3607279 blast hits in iteration 12
#          #2021-03-12 09:33:26 INFO: Writing assembler input files: iteration 12
#          #2021-03-12 12:06:53 INFO: Assembling shards with spades: iteration 12
#          #2021-03-12 13:07:24 INFO: Saving assembled contigs: iteration 12
#          #2021-03-12 16:25:55 INFO: Creating new query files: iteration 13 //This is the last line of the cycle.
#        }
#      }
#    }
#  }
#}
