#!/usr/bin/env perl
use warnings;
#This script collects cpu time, maximum memory usage, page faults, and processors used from sacct output.
$starttime = $ARGV[0]; #2021-01-20
open (OUT, "> $ARGV[1]");
$glop = `sacct --format="JobID,JobName,CPUTime,cputimeraw,MaxRSS,MaxPages,NCPUS,STATE" --starttime $starttime`;
#9655370      runenv012+   00:35:34       2134                              1  COMPLETED
#9655370.bat+      batch   00:35:34       2134   1389512K        0          1  COMPLETED 
#9655370.ext+     extern   00:35:35       2135       928K        0          1  COMPLETED
@lines = split(/\n/, $glop);
$getflag = 0;
for ($i = 2; $i < scalar(@lines); $i++) {
  @vars = split(/\s+/, $lines[$i]);
  if ($vars[0] eq "9655370") {printf "vars for 9655370 has %d fields\n", scalar(@vars);}
  if (!exists($jobids{$vars[0]}) && scalar(@vars) == 6) {
    if ($vars[-1] eq "COMPLETED") {
      $jobids{$vars[0]} = 1;
      $currjob = $vars[0];
      $jobname{$currjob} = $vars[1];
      $maxmem{$currjob} = 0;
      $maxtime{$currjob} = $vars[3];
      $maxpages{$currjob} = 0;
      $getflag = 1;
    }
    else {$getflag = 0;}
  }
  elsif ($getflag == 1 && scalar(@vars) == 8) {
    if ($vars[3] > $maxtime{$currjob}) {$maxtime{$currjob} = $vars[3];}
    $vars[4] =~ s/K/000/;
    $vars[4] =~ s/M/000000/;
    $vars[4] =~ s/G/000000000/;
    if ($vars[4] > $maxmem{$currjob}) {$maxmem{$currjob} = $vars[4];}
    if ($vars[5] > $maxpages{$currjob}) {$maxpages{$currjob} = $vars[5];}
  }
}
for $key (sort {$a <=> $b} (keys(%jobids))) {
  print OUT "jobid = $key jobname = $jobname{$key} maxtime = $maxtime{$key} maxmem = $maxmem{$key} maxpages = $maxpages{$key}\n";
}
