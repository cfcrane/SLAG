#!/usr/bin/env perl
use warnings;
#This script quickly calculated seconds per cycle for SLAG runs.
open (INP, "< $ARGV[0]"); #slagZeafiles0530a.txt, slagStanleyfiles0530a.txt, slaghemiStanleyfiles0530a.txt
open (OUT, "> $ARGV[1]");
$daysec = 24 * 60 * 60;
while ($line = <INP>) {
  chomp $line;
  #Zeahexokinasesbm_assembly19.subset.contigs      Zeahexokinasesbm_assembly20.subset.contigs slurm-10156004.out
  ($penultimatecontigfile, $lastcontigfile, $slurmfile) = split(/\t/, $line);
  @pstats = stat($penultimatecontigfile);
  @cstats = stat($lastcontigfile);
  $lastcyclelength = $cstats[9] - $pstats[9];
  open (SLR, "< $slurmfile");
  while ($line = <SLR>) {
    chomp $line;
    #Cycle 0 started at 15:28:50
    if ($line =~ m/Cycle/ && $line =~ m/started at/) {
      @vars = split(/\s+/, $line);
      @tars = split(/:/, $vars[-1]);
      $hours{$vars[1]} = $tars[0];
      $minutes{$vars[1]} = $tars[1];
      $seconds{$vars[1]} = $tars[2];
    }
  }
  $days = 0;
  $prevhour = $hours{"0"};
  for $key (sort {$a <=> $b} (keys(%hours))) {
    if ($hours{$key} < $prevhour) {$days++;}
    $prevhour = $hours{$key};
    $lastcycle = $key;
  }
  print OUT "hours{0) = $hours{0} minutes{0} = $minutes{0} seconds{0} = $seconds{0}\n";
  print OUT "days = $days last cycle = $lastcycle hours{$lastcycle} = $hours{$lastcycle} minutes{$lastcycle} = $minutes{$lastcycle} seconds{$lastcycle} = $seconds{$lastcycle} last cycle went $lastcyclelength seconds\n";
  $done = $lastcyclelength + $days * $daysec + 3600 * $hours{$lastcycle} + 60 * $minutes{$lastcycle} + $seconds{$lastcycle};
  $start = 3600 * $hours{0} + 60 * $minutes{0} + $seconds{0};
  $duration = $done - $start;
  $ncycles = $lastcycle + 1;
  $secpercycle{$lastcontigfile} = $duration / $ncycles;
}
for $key (sort(keys(%secpercycle))) {print OUT "$key\t$secpercycle{$key}\n";}
