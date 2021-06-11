#!/usr/bin/env perl
use warnings;
#This script collects counts and lengths of contigs from SLAG (localassembly1115.pl) runs.
open (LIS, "< $ARGV[0]"); #ZtUNcollectedfiles.txt
open (OUT, "> $ARGV[1]");
$searchstringa = $ARGV[2]; #individual accession, ALP48286.1_unicycler_increment_intact
$searchstringb = $ARGV[3]; #subset.contigs
$searchstringc = $ARGV[4]; #contigs, assembly.fasta
$searchstringd = $ARGV[5]; #retrieved
@files = <LIS>;
for ($i = 0; $i < scalar(@files); $i++) {
  chomp $files[$i];
  @stats = stat($files[$i]);
  $times{$files[$i]} = $stats[9];
}
$nsisters = 0; $nkeys = 0;
for ($k = 0; $k < scalar(@files); $k++) {
  if ($files[$k] =~ m/$searchstringa/ && $files[$k] =~ m/$searchstringb/) {
    #print OUT "DIAG: $files[$k] matches $searchstringa and $searchstringb\n";
    $nkeys++;
    @vars = split(/\//, $files[$k]);
    $stem = "$vars[0]";
    $idstem = "$vars[0]";
    for ($i = 1; $i < scalar(@vars) - 1; $i++) {$stem .= "/$vars[$i]";}
    #print OUT "Diag: stem = $stem\n";
    $mintimediff = 1000000000;
    for ($m = 0; $m < scalar(@files); $m++) {
      if ($files[$m] =~ m/$stem/) {
        #print OUT "diag: $files[$m] matches stem $stem\n";
        if ($files[$m] =~ m/$searchstringc/) {
          $timediff = abs($times{$files[$m]} - $times{$files[$k]});
          if ($timediff <= 1) {$sister{$files[$k]} = $files[$m]; $nsisters++;}
        }
      }
    }
    for ($i = 1; $i < scalar(@vars) - 2; $i++) {$idstem .= "/$vars[$i]";}
    $mintimediff = 1000000000;
    for ($m = 0; $m < scalar(@files); $m++) {
      if ($files[$m] =~ m/$idstem/) {
        if ($files[$m] =~ m/$searchstringd/) {
          $timediff = abs($times{$files[$m]} - $times{$files[$k]});
          if ($timediff < $mintimediff) {
            $mintimediff = $timediff;
            if ($timediff <= 1) {$brother{$files[$k]} = $files[$m];}
          }
        }
      }
    }
  }
}
$nbrothers = 0;
for $key (keys(%brother)) {$nbrothers++;}
print OUT "There are $nkeys files that match $searchstringb and $nsisters sister files that match $searchstringc\n";
print OUT "There are $nbrothers files that match $searchstringd\n";
if ($nkeys != $nsisters) {
  die "The file naming convention does not fit this script, which requires otherwise identical names with different suffixes, so some files are missing.\n";
}
if ($nkeys != $nbrothers) {
  print "Brother files were not necessarily found\; nkeys = $nkeys nbrothers = $nbrothers\n";
}
$totncontigs = 0; $totaltncontigs = 0;
for $key (sort {$times{$a} <=> $times{$b}} (keys(%sister))) {
  open (INP, "< $key");
  $ncontigs = 0;
  $seqlength = 0;
  while ($line = <INP>) {
    chomp $line;
    if ($line =~ m/>/) {
      if ($seqlength > 0) {push(@{$seqlengths{$key}}, $seqlength);}
      $ncontigs++;
      $seqlength = 0;
      #>Contig1
    }
    else {$seqlength += length($line);}
  }
  $totncontigs += $ncontigs;
  if ($seqlength > 0) {push(@{$seqlengths{$key}}, $seqlength);}
  close (INP);
  $altfile = $sister{$key};
  $altncontigs = 0;
  open (INP, "< $altfile") || die "$altfile could not be opened\n";
  $seqlength = 0;
  while ($line = <INP>) {
    if ($line =~ m/>/) {
      if ($seqlength > 0) {push(@{$altseqlengths{$key}}, $seqlength);}
      $altncontigs++;
      $seqlength = 0;
    }
    else {
      chomp $line;
      $seqlength += length($line);
    }
  }
  $totaltncontigs += $altncontigs;
  if ($seqlength > 0) {push(@{$altseqlengths{$key}}, $seqlength);}
  close (INP);
  if ($brother{$key} ne "NULL") {
    $retfile = $brother{$key};
    $nidentities = 0;
    open (INP, "< $retfile");
    while ($line = <INP>) { #tabular blast output, specific to localassembly1115.pl
      chomp $line;
      $nidentities++;
      #TraesCS1D01G398900.1    NODE_1_length_1549_cov_5.739130 8.87e-108       TraesCS1D01G398900.1    374     202     226     218     218     8       96.460  96.46   0       0       NODE_1_length_1549_cov_5.739130 1/-1    1       1       226     -1      891     666
      @vars = split(/\s+/, $line);
      push(@{$identities{$key}}, $vars[10]);
    }
    close (INP);
    #/scratch/halstead/c/ccrane/slagtests/TraesCS1A01G408800_cap3_increment_intactextracted3.fasta.cap_0_1_2_3.subset.contigs matching       1       1530    total   1530    identities      93.033
    @tars = split(/\//, $key);
    @pars = split(/_/, $tars[-1]);
    if (exists($longestmatchingcontiglength{$pars[0]})) {
      for $i (0 .. $#{$seqlengths{$key}}) {
        if ($seqlengths{$key}[$i] > $longestmatchingcontiglength{$pars[0]}) {
          $longestmatchingcontiglength{$pars[0]} = $seqlengths{$key}[$i];
          #$savedfile = $key;
        }
      }
    }
    else {$longestmatchingcontiglength{$pars[0]} = $seqlengths{$key}[0];}
    if (exists($longestaltcontiglength{$pars[0]})) {
      for $i (0 .. $#{$altseqlengths{$key}}) {
        if ($altseqlengths{$key}[$i] > $longestaltcontiglength{$pars[0]}) {
          $longestaltcontiglength{$pars[0]} = $altseqlengths{$key}[$i];
        }
      }
    }
    else {$longestaltcontiglength{$pars[0]} = $altseqlengths{$key}[0];}
    if (exists($highestidentity{$pars[0]})) {
      for $i (0 .. $#{$identities{$key}}) {
        if ($identities{$key}[$i] > $highestidentity{$pars[0]}) {
          $highestidentity{$pars[0]} = $identities{$key}[$i];
        }
      }
    }
    else {$highestidentity{$pars[0]} = $identities{$key}[0];}
  }
  else {
    $message = "$key has no brother blast output.\n";
    push (@{$identities{$key}}, $message);
  }
  print OUT "$key matching\t$ncontigs";
  for $i (0 .. $#{$seqlengths{$key}}) {print OUT "\t$seqlengths{$key}[$i]";}
  print OUT "\ttotal";
  for $i (0 .. $#{$altseqlengths{$key}}) {print OUT "\t$altseqlengths{$key}[$i]";}
  print OUT "\tidentities";
  for $i (0 .. $#{$identities{$key}}) {print OUT "\t$identities{$key}[$i]";}
  print OUT "\n";
}
$lkeys = 0;
$totlongestlengths = 0;
$longestlongestlength = -1;
$shortestlongestlength = 1000000000;
for $key (sort(keys(%longestmatchingcontiglength))) {
  print OUT "$key longest matching target is $longestmatchingcontiglength{$key}\n";
  $lkeys++;
  $totlongestlengths += $longestmatchingcontiglength{$key};
  if ($longestmatchingcontiglength{$key} > $longestlongestlength) {$longestlongestlength = $longestmatchingcontiglength{$key};}
  if ($longestmatchingcontiglength{$key} < $shortestlongestlength) {
    $shortestlongestlength = $longestmatchingcontiglength{$key};
  }
}
if ($lkeys > 1) {
  @temparray = ();
  for $key (sort {$longestmatchingcontiglength{$a} <=> $longestmatchingcontiglength{$b}} (keys(%longestmatchingcontiglength))) {
    push (@temparray, $longestmatchingcontiglength{$key});
  }
  $v = int($lkeys / 2);
  if ($lkeys % 2 == 0) {$medianlongestmatchingcontiglength = ($temparray[$v-1] + $temparray[$v]) / 2;}
  else {$medianlongestmatchingcontiglength = $temparray[$v];}
}
else {
  $medianlongestmatchingcontiglength = "mean";
}
$akeys = 0;
$totaltlongestlengths = 0;
$longestaltlongestlength = -1;
$shortestaltlongestlength = 1000000000;
for $key (keys(%longestaltcontiglength)) {
  $akeys++; $totaltlongestlengths += $longestaltcontiglength{$key};
  if ($longestaltcontiglength{$key} > $longestaltlongestlength) {$longestaltlongestlength = $longestaltcontiglength{$key};}
  if ($longestaltcontiglength{$key} < $shortestaltlongestlength) {$shortestaltlongestlength = $longestaltcontiglength{$key};}
}
if ($akeys > 1) {
  @temparray = ();
  for $key (sort {$longestaltcontiglength{$a} <=> $longestaltcontiglength{$b}} (keys(%longestaltcontiglength))) {
    push (@temparray, $longestaltcontiglength{$key});
  }
  $v = int($akeys / 2);
  if ($akeys % 2 == 0) {$medianlongestaltcontiglength = ($temparray[$v-1] + $temparray[$v]) / 2;}
  else {$medianlongestaltcontiglength = $temparray[$v];}
}
else {$medianlongestaltcontiglength = "mean";}
$ikeys = 0;
$tothighestidentities = 0;
$highesthighestidentity = -1;
$lowesthighestidentity = 200;
for $key (keys(%highestidentity)) {
  $ikeys++; $tothighestidentities += $highestidentity{$key};
  if ($highestidentity{$key} < $lowesthighestidentity) {$lowesthighestidentity = $highestidentity{$key};}
  if ($highestidentity{$key} > $highesthighestidentity) {$highesthighestidentity = $highestidentity{$key};}
}
if ($ikeys > 1) {
  @temparray = ();
  for $key (sort {$highestidentity{$a} <=> $highestidentity{$b}} (keys(%highestidentity))) {push (@temparray, $highestidentity{$key});}
  $v = int($ikeys / 2);
  if ($ikeys % 2 == 0) {$medianhighestidentity = ($temparray[$v-1] + $temparray[$v]) / 2;}
  else {$medianhighestidentity = $temparray[$v];}
}
else {$medianhighestidentity = "mean";}
if ($lkeys > 0) {$meanlongestlength = $totlongestlengths / $lkeys;}
else {$meanlongestlength = "undef";}
if ($ikeys > 0) {$meanhighestidentity = $tothighestidentities / $ikeys;}
else {$meanhighestidentity = "undef";}
if ($akeys > 0) {$meanaltlongestlength = $totaltlongestlengths / $akeys;}
else {$meanaltlongestlength = "undef";}
$stdvlongestmatchingcontiglength = 0;
for $key (keys(%longestmatchingcontiglength)) {
  $stdvlongestmatchingcontiglength += ($longestmatchingcontiglength{$key} - $meanlongestlength)**2;
}
if ($lkeys > 1) {
  $stdvlongestmatchingcontiglength /= ($lkeys - 1);
  $stdvlongestmatchingcontiglength = sqrt($stdvlongestmatchingcontiglength);
  $sterrlongestmatchingcontiglength = $stdvlongestmatchingcontiglength / sqrt($lkeys);
}
else {
  $stdvlongestmatchingcontiglength = "undef";
  $sterrlongestmatchingcontiglength = "undef";
}
$stdvaltlongestmatchingcontiglength = 0;
for $key (keys(%longestaltcontiglength)) {
  $stdvaltlongestmatchingcontiglength += ($longestaltcontiglength{$key} - $meanaltlongestlength)**2;
}
if ($akeys > 1) {
  $stdvaltlongestmatchingcontiglength /= ($akeys - 1);
  $stdvaltlongestmatchingcontiglength = sqrt($stdvaltlongestmatchingcontiglength);
  $stderraltlongestmatchingcontiglength = $stdvaltlongestmatchingcontiglength / sqrt($akeys);
}
else {
  $stdvaltlongestmatchingcontiglength = "undef";
  $stderraltlongestmatchingcontiglength = "undef";
}
print OUT "Number of matching contigs = $lkeys number of all contigs = $akeys number of identities = $ikeys\n";
print OUT "Mean longest matching contig length for $searchstringa is $meanlongestlength\n";
print OUT "Standard deviation of longest matching contig length for $searchstringa is $stdvlongestmatchingcontiglength\n";
print OUT "Standard error of the mean for longest matching contig length for $searchstringa is $sterrlongestmatchingcontiglength\n";
print OUT "Median longest matching contig length for $searchstringa is $medianlongestmatchingcontiglength\n";
print OUT "Range of longest matching contig lengths is $shortestlongestlength to $longestlongestlength\n";
print OUT "Mean longest any contig length for $searchstringa is $meanaltlongestlength\n";
print OUT "Standard deviation of longest any contig length for $searchstringa is $stdvaltlongestmatchingcontiglength\n";
print OUT "Standard error of the mean for longest any contig length for $searchstringa is $stderraltlongestmatchingcontiglength\n";
print OUT "Median longest any-contig length for $searchstringa is $medianlongestaltcontiglength\n";
print OUT "Range of longest any contig length for $searchstringa is $shortestaltlongestlength to $longestaltlongestlength\n";
print OUT "Mean highest identity for $searchstringa is $meanhighestidentity\n";
print OUT "Median highest identity for $searchstringa is $medianhighestidentity\n";
print OUT "Range of highest identities is $lowesthighestidentity to $highesthighestidentity\n";
print OUT "There are $nkeys .subset.contigs files, $nsisters parent .contigs files, and $nidentities blast vs retrieved files.\n";
$meannmatching = $totncontigs / $nkeys;
$meannall = $totaltncontigs / $nsisters;
print OUT "diag: totncontigs = $totncontigs nkeys = $nkeys totaltncontigs = $totaltncontigs nsisters = $nsisters\n";
print OUT "Mean number of contigs produced = $meannall mean that match the target = $meannmatching\n";
print OUT "Longest matching contig lengths ordered by target identifier for comparison with other assemblies:\n";
for $key (sort {$a cmp $b} (keys(%longestmatchingcontiglength))) {print OUT "$key\t$longestmatchingcontiglength{$key}\n";}
