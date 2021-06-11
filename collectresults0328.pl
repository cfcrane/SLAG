#!/usr/bin/env perl
use warnings;
#This script collected results from local-assembly files in its input list.
open (LIS, "< $ARGV[0]"); #filelist03282021a.txt
open (KEY, "< $ARGV[1]"); #searchtermtovariable328.txt
open (OUT, "> $ARGV[2]");
while ($line = <KEY>) {
  chomp $line;
  ($searchterm, $name) = split(/\t/, $line);
  push(@searchterms, $searchterm);
  push(@names, $name);
}
print OUT "Code\tFile";
for ($i = 0; $i < scalar(@searchterms); $i++) {
  print OUT "\t$names[$i]";
}
print OUT "\n";
while ($record = <LIS>) {
  chomp $record;
  ($code, $file) = split(/\s+/, $record);
  if (exists($codecounts{$code})) {$codecounts{$code}++;}
  else {
    $codecounts{$code} = 1;
    for ($i = 0; $i < scalar(@searchterms); $i++) {
      $codemeans{$code}[$i] = 0;
      $codeses{$code}[$i] = 0;
    }
  }
  open (INP, "< $file");
  while ($line = <INP>) {
    chomp $line;
    @vars = split(/\s+/, $line);
    for ($i = 0; $i < scalar(@searchterms); $i++) {
      if ($line =~ m/$searchterms[$i]/) {
        #print OUT "diag: $searchterms[$i] is matched by $line\n";
        for ($j = 0; $j < scalar(@vars) - 1; $j++) {
          if ($vars[$j] eq "=") {
            #print OUT "j = $j vars[$j+1] = $vars[$j+1]\n";
            if ($vars[$j+1] =~ m/\d+/) {
              $variables[$i] = $vars[$j+1];
              $codemeans{$code}[$i] += $vars[$j+1];
            }
            else {
              #We have a missing datum.
              $variables[$i] = "NA";
            }
          }
        }
      }
    }
  }
  print OUT "$code\t$file";
  for ($i = 0; $i < scalar(@variables); $i++) {print OUT "\t$variables[$i]";}
  print OUT "\n";
}
print OUT "Means:\n";
for ($i = 0; $i < scalar(@searchterms); $i++) {
  print OUT $searchterms[$i];
  for $key (sort(keys(%codemeans))) {
    $codemeans{$key}[$i] /= $codecounts{$key};
    print OUT "\t$codemeans{$key}[$i]";
  }
  print OUT "\n";
}
seek(LIS, 0, 0);
while ($record = <LIS>) {
  chomp $record;
  ($code, $file) = split(/\s+/, $record);
  open (INP, "< $file");
  while ($line = <INP>) {
    chomp $line;
    for ($i = 0; $i < scalar(@searchterms); $i++) {
      if ($line =~ m/$searchterms[$i]/) {
        #print OUT "diag: $searchterms[$i] is matched by $line\n";
        for ($j = 0; $j < scalar(@vars) - 1; $j++) {
          if ($vars[$j] eq "=") {
            #print OUT "j = $j vars[$j+1] = $vars[$j+1]\n";
            if ($vars[$j+1] =~ m/\d+/) {
              $variables[$i] = $vars[$j+1];
              $codeses{$code}[$i] += ($vars[$j+1] - $codemeans{$code}[$i]) ** 2;
            }
          }
        }
      }
    }
  }
}
for $key (keys(@codemeans)) {
  for $i (0 .. $#{$codemeans{$key}) {
    if ($codecounts{$key} > 1) {
      $codeses{$key}[$i] /= ($codecounts{$key} - 1);
      $codeses{$key}[$i] = sqrt($codeses{$key}[$i]) / sqrt($codecounts{$key});
      
    }
    else {$codeses{$key}[$i] = "-";}
  }
}
for ($i = 0; $i < scalar(@searchterms); $i++) {
  print OUT $searchterms[$i];

}
