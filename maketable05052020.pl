#!/usr/bin/env perl
use warnings;
#This script built the Zea results table for the SLAG paper.
open (LIA, "< $ARGV[0]"); #programlist05052020a.txt
open (LIB, "< $ARGV[1]"); #filelist05052020a.txt
open (GRP, "< $ARGV[2]"); #Zeagroups.txt
open (OUT, "> $ARGV[3]"); #slagtable2.txt
while ($program = <LIA>) {
  chomp $program;
  push (@programs, $program);
}
while ($file = <LIB>) {
  chomp $file;
  push (@files, $file);
  for ($i = 0; $i < scalar(@programs); $i++) {
    if ($file =~ m/$programs[$i]/) {$code{$file} = $programs[$i];}
  }
}
while ($group = <GRP>) { #The following assumes that no group name is a substring of another group name.
  chomp $group;
  push (@groups, $group);
}
for ($i = 0; $i < scalar(@programs); $i++) {
  print OUT $programs[$i];
  for ($j = 0; $j < scalar(@files); $j++) {
    if ($files[$j] =~ m/$programs[$i]/) {
      print OUT " $files[$j]\n";
      open (INP, "< $files[$j]");
      while ($line = <INP>) {
        if ($line =~ m/diag:/) {next;}
        chomp $line;
        @vars = split(/\t/, $line);
        #sucrosesynthases        72015   0.941198688603847
        for ($k = 0; $k < scalar(@groups); $k++) {
          if ($vars[0] eq $groups[$k]) {
            $longests{$groups[$k]}{$programs[$i]} = $vars[1];
            $identities{$groups[$k]}{$programs[$i]} = $vars[2];
          }
        }
        @tars = split(/ /, $line);
        if ($tars[0] eq "hitting") {
          #hitting contigs = 85 unhitting contigs = 0
          $hittingcontigs{$programs[$i]} = $tars[3];
          $nonhittingcontigs{$programs[$i]} = $tars[7]; 
        }
      }
    }
  }
}
print OUT "Group";
for ($i = 0; $i < scalar(@programs); $i++) {print OUT "\t$programs[$i]";}
print OUT "\n";
for ($i = 0; $i < scalar(@groups); $i++) {
  print OUT $groups[$i];
  for ($j = 0; $j < scalar(@programs); $j++) {
    if (exists($longests{$groups[$i]}{$programs[$j]})) {printf OUT "\t%d", $longests{$groups[$i]}{$programs[$j]};}
    else {print OUT "\t-";}
  }
  print OUT "\n";
  for ($j = 0; $j < scalar(@programs); $j++) {
    if (exists($identities{$groups[$i]}{$programs[$j]})) {printf OUT "\t%.4f", $identities{$groups[$i]}{$programs[$j]};}
    else {print OUT "\t-";}
  }
  print OUT "\n";
}
print OUT "Matching contigs";
for ($i = 0; $i < scalar(@programs); $i++) {print OUT "\t$hittingcontigs{$programs[$i]}";}
print OUT "\nNonmatching contigs";
for ($i = 0; $i < scalar(@programs); $i++) {
  if (exists($nonhittingcontigs{$programs[$i]})) {print OUT "\t$nonhittingcontigs{$programs[$i]}";}
  else {print OUT "\t0";}
}
print OUT "\n";
