#!/usr/bin/env perl
use warnings;
#This script collects counts and lengths of contigs from SLAG (localassembly1115.pl) runs.
$dir = $ARGV[0]; #/scratch/halstead/c/ccrane/slagtests
if ($dir =~ m/\/$/) {chop $dir;}
if (scalar(@ARGV) != 13) {die "This script needs 13 arguments.\n";}
open (OUT, "> $ARGV[1]");
$searchstringa = $ARGV[2]; #unicycler_increment_intact_workstem
$searchstringb = $ARGV[3]; #subset.contigs
$searchstringc = $ARGV[4]; #assembly.fasta
$searchstringd = $ARGV[5]; #retrieved
$searchstringe = $ARGV[6]; #vs
$searchstringf = $ARGV[7]; #unicycler
$searchstringg = $ARGV[8]; #TraesCS
$searchstringh = $ARGV[9]; #unicycler.log
$searchstringi = $ARGV[10]; #HU
$searchstringj = $ARGV[11]; #unicretrieved
$searchstringk = $ARGV[12]; #extracted
@extfiles = glob("$dir/*$searchstringf*$searchstringk*");
$scalarextfiles = scalar(@extfiles);
print OUT "There are $scalarextfiles files in array extfiles.\n";
#TraesCS1B01G397300.1_unicycler_increment_intactHUextracted1.fasta
for ($i = 0; $i < scalar(@extfiles); $i++) {
  if ($extfiles[$i] =~ m/$searchstringg/) {
    $currfile = $extfiles[$i];
    @stats = stat($currfile);
    $exttimes{$currfile} = $stats[9];
  }
}
$nexttimes = keys(%exttimes);
print OUT "The exttimes hash has $nexttimes keys.\n";
@subdirs = grep {-d} glob("$dir/*$searchstringa*"); #from https://stackoverflow.com/questions/1692492/getting-names-of-directories-under-given-path
$ndirs = 0;
for ($i = 0; $i < scalar(@subdirs); $i++) {
  #print OUT "diag: $subdirs[$i]\n";
  $subdir = $subdirs[$i];
  $ndirs++;
  opendir (DIR, $subdir);
  while ($filename = readdir(DIR)) {
    $currfile = $subdir."/".$filename;
    $curracc = "NULL";
    @vars = split(/\/|_/, $currfile); #This ASSUMES that the accession is set off with a slash or underscore.
    for ($j = 0; $j < scalar(@vars); $j++) {
      if ($vars[$j] =~ m/$searchstringg/) {
        $curracc = $vars[$j];
        if ($curracc =~ m/$searchstringi/) {$curracc =~ s/$searchstringi//g;} #expedient but limiting
        if ($curracc =~ m/$searchstringj/) {$curracc =~ s/$searchstringj//g;} #expedient but limiting
        if (!exists($children{$curracc})) {$children{$curracc} = ();}
      }
    }
    if ($curracc ne "NULL") {
      if ($currfile =~ m/$searchstringb$/) {
        @stats = stat($currfile);
        $times{$currfile} = $stats[9];
        push(@{$children{$curracc}}, $currfile);
      }
      if ($currfile =~ m/$searchstringc$/) {
        @stats = stat($currfile);
        $alttimes{$currfile} = $stats[9];
        push(@{$children{$curracc}}, $currfile);
      }
      if ($currfile =~ m/$searchstringh$/) {
        @stats = stat($currfile);
        $logtimes{$currfile} = $stats[9];
        push(@{$children{$curracc}}, $currfile);
      }
    }
  }
  closedir(DIR);
}
$nsecblastfiles = 0;
@allfiles = glob("$dir/*$searchstringf*$searchstringe*$searchstringd*");
for ($i = 0; $i < scalar(@allfiles); $i++) {
  if ($allfiles[$i] =~ m/$searchstringf/) {
    if ($allfiles[$i] =~ m/$searchstringe/) {
      if ($allfiles[$i] =~ m/$searchstringd/) {
        $nsecblastfiles++;
        $currfile = $allfiles[$i];
        @stats = stat($currfile);
        $rettimes{$currfile} = $stats[9];
        $curracc = "NULL";
        @vars = split(/\/|_/, $currfile);
        for ($j = 0; $j < scalar(@vars); $j++) {
          if ($vars[$j] =~ m/$searchstringg/) {
            $curracc = $vars[$j];
            if ($curracc =~ m/$searchstringi/) {$curracc =~ s/$searchstringi//g;} #expedient but limiting
            if ($curracc =~ m/$searchstringj/) {$curracc =~ s/$searchstringj//g;} #expedient but limiting
          }
        }
        if ($curracc ne "NULL") {push(@{$children{$curracc}}, $currfile);}
      }
    }
  }
}
print OUT "\nThere were $nsecblastfiles blast outputs involving $searchstringf and $searchstringd\n\n";
print OUT "List of accessions:\n";
for $key (sort(keys(%children))) {print OUT " $key";}
print OUT "\n\n";
print OUT "For the times hash:\n";
for $key (sort {$times{$a} <=> $times{$b}} (keys(%times))) {print OUT "$key\t$times{$key}\n";}
print OUT "\nFor the alttimes hash:\n";
for $key (sort {$alttimes{$a} <=> $alttimes{$b}} (keys(%alttimes))) {print OUT "$key\t$alttimes{$key}\n";}
print OUT "\nFor the rettimes hash:\n";
for $key (sort {$rettimes{$a} <=> $rettimes{$b}} (keys(%rettimes))) {print OUT "$key\t$rettimes{$key}\n";}
print OUT "\nFor the logtimes hash:\n";
for $key (sort {$logtimes{$a} <=> $logtimes{$b}} (keys(%logtimes))) {print OUT "$key\t$logtimes{$key}\n";}
print OUT "\n";
$nfailedassemblies = 0;
for $key (sort(keys(%logtimes))) {
  $curracc = "NULL";
  @vars = split(/\/|_/, $key); #This ASSUMES that the accession is set off with a slash or underscore.
  for ($i = 0; $i < scalar(@vars); $i++) {
    if ($vars[$i] =~ m/$searchstringg/) {
      $curracc = $vars[$i];
    }
  }
  if ($curracc !~ m/$searchstringg/) {next;}
  open (LOG, "< $key");
  while ($line = <LOG>) {
    if ($line =~ m/Error: miniasm assembly failed/) {
      print OUT "Assembly failed for $key\n";
      $nfailedassemblies++;
    }
  }
}
print OUT "There were $nfailedassemblies targets where assembly failed.\n\n";
for $key (sort(keys(%children))) {
  $nsubsetcontigfiles{$key} = 0;
  $nassemblyfastafiles{$key} = 0;
  $nretrievedblastfiles{$key} = 0;
  print OUT " $key:\n";
  for $i (0 .. $#{$children{$key}}) {
    print OUT "\t$children{$key}[$i]\n";
    if ($children{$key}[$i] =~ m/$searchstringb/) {$nsubsetcontigfiles{$key}++;}
    elsif ($children{$key}[$i] =~ m/$searchstringc/) {$nassemblyfastafiles{$key}++;}
    elsif ($children{$key}[$i] =~ m/$searchstringd/) {$nretrievedblastfiles{$key}++;}
  }
  print OUT "subset.contigs = $nsubsetcontigfiles{$key} assembly.fasta = $nassemblyfastafiles{$key} retrieved blast = $nretrievedblastfiles{$key}\n";
  $longestsubsetcontiglength{$key} = 0;
  $longestsubsetcontig{$key} = "NULL";
  $longestsubsetcontigfile{$key} = "NULL";
  $nsstcontigs{$key} = 0;
  $nsstfiles{$key} = 0;
  for $i (0 .. $#{$children{$key}}) {
    if ($children{$key}[$i] =~ m/$searchstringb/) {
      $currfile = $children{$key}[$i];
      $nsstfiles{$key}++;
      open (INP, "< $currfile");
      $currlength = 0;
      while ($line = <INP>) {
        chomp $line;
        if ($line =~ m/>/) {
          #>1
          $nsstcontigs{$key}++;
          if ($currlength > $longestsubsetcontiglength{$key}) {
            $longestsubsetcontiglength{$key} = $currlength;
            $longestsubsetcontigfile{$key} = $currfile;
            ($contignumber = $line) =~ s/>//;
            $longestsubsetcontig{$key} = $contignumber;
          }
          $currlength = 0;
        }
        else {$currlength += length($line);}
      }
    }
  }
  close(INP);
  for $j (0 .. $#{$children{$key}}) {
    if ($children{$key}[$j] =~ m/$searchstringd/) {
      $retfile = $children{$key}[$j];
      if (abs($times{$currfile} - $rettimes{$retfile}) < 2) {
        open (INP, "< $retfile");
        while ($pine = <INP>) {
          chomp $pine;
          #TraesCS1B01G451600.1    1       7.95e-96        TraesCS1B01G451600.1    339     183     233     217     217     12      93.133  93.13   4       3       1       1/-1    1       1       232     -1      22870   22641
          @tars = split(/\t/, $pine);
          if ($tars[1] eq $longestsubsetcontig{$key}) {$matchedidentity{$key} = $tars[10];}
        }
        close(INP);
      }
    }
  }
  $longestassembcontiglength{$key} = 0;
  $longestassembcontig{$key} = "NULL";
  $longestassembcontigfile{$key} = "NULL";
  $naltcontigs{$key} = 0;
  $naltfiles{$key} = 0;
  for $i (0 .. $#{$children{$key}}) {
    if ($children{$key}[$i] =~ m/$searchstringc/) {
      $altfile = $children{$key}[$i];
      $naltfiles{$key}++;
      open (INP, "< $altfile");
      $currlength = 0;
      while ($line = <INP>) {
        chomp $line;
        if ($line =~ m/>/) {
          #>1
          $naltcontigs{$key}++;
          ($contignumber = $line) =~ s/>//;
          if ($currlength > $longestassembcontiglength{$key}) {
            $longestassembcontiglength{$key} = $currlength;
            $longestassembcontigfile{$key} = $altfile;
            $longestassembcontig{$key} = $contignumber;
          }
          $currlength = 0;
        }
        else {$currlength += length($line);}
      }
      if ($currlength > $longestassembcontiglength{$key}) {
        $longestassembcontiglength{$key} = $currlength;
        $longestassembcontigfile{$key} = $altfile;
        $longestassembcontig{$key} = $contignumber;
      }
    }
  }
  $highestidentity{$key} = 0;
  $highestidentitycontig{$key} = "NULL";
  $highestidentityfile{$key} = "NULL";
  for $i (0 .. $#{$children{$key}}) {
    if ($children{$key}[$i] =~ m/$searchstringd/) {
      $retfile = $children{$key}[$i];
      open (INP, "< $retfile");
      while ($line = <INP>) {
        chomp $line;
        @bars = split(/\t/, $line);
        if ($bars[10] > $highestidentity{$key}) {
          $highestidentity{$key} = $bars[10];
          $highestidentitycontig{$key} = $bars[1];
          $highestidentityfile{$key} = $retfile;
        }
      }
      close(INP);
    }
  }
}
$meanlongestsstcontiglength = 0;
$meanlongestpinelength = 0;
$npinelengths = 0;
for $key (keys(%longestsubsetcontigfile)) {
  if ($longestsubsetcontigfile{$key} ne "NULL") {
    print OUT "DAAG: key = $key longest contig = $longestsubsetcontiglength{$key} longestsubsetcontigfile{key} = $longestsubsetcontigfile{$key}\n";
    $meanlongestsstcontiglength += $longestsubsetcontiglength{$key};
    $currfile = $longestsubsetcontigfile{$key};
    for $i (0 .. $#{$children{$key}}) {
      $altfile = $children{$key}[$i];
      if ($altfile =~ m/$searchstringc/) {
        if (abs($times{$currfile} - $alttimes{$altfile}) < 2) {
          $altmatchedtolongestsubsetcontigfile{$key} = $altfile;
          open (INP, "< $altfile");
          $npinelengths++;
          $maxtemplength = 0;
          $templength = 0;
          while ($pine = <INP>) {
            chomp $pine;
            if ($pine =~ m/^>/) {
              if ($templength > $maxtemplength) {$maxtemplength = $templength;}
              $templength = 0;
            }
            else {$templength += length($pine);}
          }
          if ($templength > $maxtemplength) {$maxtemplength = $templength;}
          close(INP);
          $longestpinelength{$key} = $maxtemplength;
          $meanlongestpinelength += $maxtemplength;
          print OUT "DUGG: $key $maxtemplength $altfile\n";
        }
      }
    }
  }
}
$meanlongestpinelength /= $npinelengths;
print OUT "meanlongestpinelength = $meanlongestpinelength for $npinelengths longest assembly.fasta contigs matched to subset.contigs contigs\n";
$nlongestsubsetcontigfiles = 0;
for $key (keys(%longestsubsetcontigfile)) {
  if ($longestsubsetcontigfile{$key} ne "NULL") {$nlongestsubsetcontigfiles++;}
}
$naltmatchedtolongestsubsetcontigfile = 0;
for $key (keys(%altmatchedtolongestsubsetcontigfile)) {
  if ($altmatchedtolongestsubsetcontigfile{$key} ne "NULL") {$naltmatchedtolongestsubsetcontigfile++;}
}
print OUT "There are $nlongestsubsetcontigfiles accessions that have a longest subset.contigs contig over all cycles.\n";
$meanlongestsstcontiglength /= $nlongestsubsetcontigfiles;
print OUT "Mean length of longest subset.contigs contig from any cycle over all accessions that produced subset.contigs = $meanlongestsstcontiglength\n";
print OUT "There are $naltmatchedtolongestsubsetcontigfile accessions that have a longest assembly.fasta contig over all cycles.\n";
$longestsstcontiglengthssd = 0; #These are for MATCHED subset.contigs and assembly.fasta files, i.e., from the same cycle.
$meanlongestpinelengthssd = 0;
for $key (keys(%longestsubsetcontigfile)) {
  if ($longestsubsetcontigfile{$key} ne "NULL") {
    $longestsstcontiglengthssd += ($meanlongestsstcontiglength - $longestsubsetcontiglength{$key})**2;
    $meanlongestpinelengthssd += ($meanlongestpinelength - $longestpinelength{$key})**2;
  }
}
$longestsstcontiglengthstddev = sqrt($longestsstcontiglengthssd / ($nlongestsubsetcontigfiles - 1));
$longestsstcontiglengthstderr = $longestsstcontiglengthstddev / sqrt($nlongestsubsetcontigfiles);
print OUT "Number of assembly.fasta files matched to longest subset.contigs files = $naltmatchedtolongestsubsetcontigfile\n";
$longestpinelengthstddev = sqrt($meanlongestpinelengthssd / ($naltmatchedtolongestsubsetcontigfile - 1));
$longestpinelengthstderr = $longestpinelengthstddev / sqrt($naltmatchedtolongestsubsetcontigfile);
print OUT "Standard deviation of longest subset.contigs contig length = $longestsstcontiglengthstddev\n";
print OUT "Standard error of mean of longest subset.contigs contig length = $longestsstcontiglengthstderr\n";
$naltmatchedtolongestsubsetcontigfile = keys(%altmatchedtolongestsubsetcontigfile);
print OUT "Mean length of longest assembly.fasta contig from cycles matched to longest subset.contigs contig = $meanlongestpinelength\n";
print OUT "Standard deviation of longest matched assembly.fasta contig length = $longestpinelengthstddev\n";
print OUT "Standard error of mean of longest matched assembly.fasta contig length = $longestpinelengthstderr\n";
$nhighestidentities = 0;
$meanhighestidentity = 0;
for $key (keys(%highestidentity)) {
  if ($highestidentity{$key} != 0) {$nhighestidentities++; $meanhighestidentity += $highestidentity{$key};}
}
$meanhighestidentity /= $nhighestidentities;
print OUT "Mean highest identity = $meanhighestidentity based on $nhighestidentities values.\n";
$totalnsstcontigs = 0;
$totalnaltcontigs = 0;
$totallongestsubsetcontiglength = 0;
$totallongestassembcontiglength = 0;
for $key (sort(keys(%naltcontigs))) {
  print OUT "For $key:\n";
  if ($nsstfiles{$key} > 0) {$meansstcontigsinkey = $nsstcontigs{$key} / $nsstfiles{$key};}
  else {$meansstcontigsinkey = "undef";}
  print OUT "nsubsetcontigs = $nsstcontigs{$key}\nnsstfiles = $nsstfiles{$key}\nmean = $meansstcontigsinkey\n";
  if ($naltfiles{$key} > 0) {$meanaltcontigsinkey = $naltcontigs{$key} / $naltfiles{$key};}
  else {$meanaltcontigsinkey = "undef";}
  print OUT "naltcontigs = $naltcontigs{$key}\nnaltfiles = $naltfiles{$key}\nmean = $meanaltcontigsinkey\n";
  print OUT "highest identity = $highestidentity{$key} its contig = $highestidentitycontig{$key} in file $highestidentityfile{$key}\n";
  if (exists($matchedidentity{$key})) {print OUT "percent identity for longest matching contig = $matchedidentity{$key}\n";}
  print OUT "longest matching contig length = $longestsubsetcontiglength{$key}\n";
  print OUT "longest matching contig is $longestsubsetcontig{$key} in file $longestsubsetcontigfile{$key}\n";
  print OUT "longest any contig length = $longestassembcontiglength{$key}\n";
  print OUT "longest any contig is $longestassembcontig{$key} in file $longestassembcontigfile{$key}\n";
  if ($longestsubsetcontiglength{$key} ne "NULL") {
    $totallongestsubsetcontiglength += $longestsubsetcontiglength{$key};
    $totalnsstcontigs += $nsstcontigs{$key};
  }
  if ($longestassembcontiglength{$key} ne "NULL") {
    $totallongestassembcontiglength += $longestassembcontiglength{$key};
    $totalnaltcontigs += $naltcontigs{$key};
  }
}
print OUT "There are $totalnsstcontigs subset.contigs contigs and $totalnaltcontigs assembly.fasta contigs.\n";
$nsubsetcontigfiles = keys(%times);
$nassembfastafiles = keys(%alttimes);
print OUT "There are $nsubsetcontigfiles subset.contigs files and $nassembfastafiles assembly.fasta files.\n";
$meansstcount = $totalnsstcontigs / $nsubsetcontigfiles;
$meanaltcount = $totalnaltcontigs / $nassembfastafiles;
print OUT "Mean contig count = $meansstcount per subset.contigs file and $meanaltcount per assembly.fasta file.\n";
#print OUT "Number of accessions = $naccs\n";
$meanlongestsstcontiglength = $totallongestsubsetcontiglength / keys(%longestsubsetcontigfile);
$meanlongestaltcontiglength = $totallongestassembcontiglength / keys(%alttimes);
print OUT "Mean longest subset.contigs contig length over all subset.contigs files = $meanlongestsstcontiglength\n";
print OUT "Mean longest assembly.fasta contig length over all assembly.fasta files = $meanlongestaltcontiglength\n";
$npairedassembfiles = keys(%times);
$meanpairedassemblongestcontiglength = 0;
for $key (keys(%longestpinelength)) {$meanpairedassemblongestcontiglength += $longestpinelength{$key};}
$meanpairedassemblongestcontiglength /= keys(%times);
print OUT "Number of assembly.fasta files that matched longest-contig subset.contigs files = $npairedassembfiles\n";
print OUT "Mean length of longest contig matched to a longest subset.contig = $meanpairedassemblongestcontiglength\n";
$meanlongestanycontig = 0;
$nnonzerolongestanycontig = 0;
for $key (keys(%longestassembcontiglength)) {
  if ($longestassembcontiglength{$key} > 0) {
    $nnonzerolongestanycontig++;
    $meanlongestanycontig += $longestassembcontiglength{$key};
  }
}
$meanlongestanycontig /= $nnonzerolongestanycontig;
$stddevlongestanycontig = 0;
for $key (keys(%longestassembcontiglength)) {
  if ($longestassembcontiglength{$key} > 0) {
    $stddevlongestanycontig += ($meanlongestanycontig - $longestassembcontiglength{$key})**2;
  }
}
$stddevlongestanycontig = sqrt($stddevlongestanycontig / ($nnonzerolongestanycontig - 1));
$stderrlongestanycontig = $stddevlongestanycontig / sqrt($nnonzerolongestanycontig);
print OUT "Mean length of any assembly.fasta longest contig = $meanlongestanycontig\n";
print OUT "Standard deviation of any assembly.fasta longest contig = $stddevlongestanycontig\n";
print OUT "Standard error of the mean of any assembly.fasta longest contig = $stderrlongestanycontig\n";
die "bailing out\n";


$lostprimarymatches = 0;
for $key (sort(keys(%alttimes))) {
  #/scratch/halstead/c/ccrane/slagtests/TraesCS1A01G399600.1_unicycler_increment_intact_workstem1/assembly.fasta
  $curracc = "NULL";
  @vars = split(/\/|_/, $key); #This ASSUMES that the accession is set off with a slash or underscore.
  for ($i = 0; $i < scalar(@vars); $i++) {
    if ($vars[$i] =~ m/$searchstringg/) {
      $curracc = $vars[$i];
      $nassembledcontigs{$curracc} = 0;
      $longestcontig{$curracc} = 0;
    }
  }
  if ($curracc eq "NULL") {die "$key gave a null value of curracc for searchstring $searchstringg!\n";}
  $brother{$key} = "NULL";
  $sister{$key} = "NULL";
  for $index (keys(%times)) {
    if ($index =~ m/$curracc/) {
      if (abs($alttimes{$key} - $times{$index}) < 2) {$brother{$key} = $index;} #brother = *subset.contigs key = *assembly.fasta
    }
  }
  for $index (keys(%rettimes)) {
    if ($index =~ m/$curracc/) {
      if (abs($alttimes{$key} - $rettimes{$index}) < 2) {$sister{$key} = $index;} #sister = 2nd*retrieved* key = *assembly.fasta
    }
  }
  if ($brother{$key} eq "NULL") {
    print OUT "$key does not have a brother file that matches $searchstringb .\n";
    $lostprimarymatches++;
  }
  if ($sister{$key} eq "NULL") {print OUT "$key does not have a sister file that matches $searchstringd .\n";}
}
print OUT "\nThere are $ndirs relevant directories under $dir\n";
print OUT "There were $lostprimarymatches cycles that lost the target sequence.\n";
for $key (sort(keys(%alttimes))) {
  $curracc = "NULL";
  @vars = split(/\/|_/, $key); #This ASSUMES that the accession is set off with a slash or underscore.
  for ($i = 0; $i < scalar(@vars); $i++) {
    if ($vars[$i] =~ m/$searchstringg/) {$curracc = $vars[$i];}
  }
  if ($curracc eq "NULL") {die "$key gave a null value of curracc for searchstring $searchstringg at the second invocation!\n";}
  if (!exists($seqlengths{$curracc})) { #This is not yet right.  We want to record all data per file.  These will be highest.
    $seqlengths{$curracc} = ();
    $altseqlengths{$curracc} = ();
    $identities{$curracc} = ();
  }
  open (SEQ, "< $key");
  
}
die "bailing out\n";

#TraesCS1B01G473100.1    1       1.66e-59        TraesCS1B01G473100.1    219     118     144     136     136     5       94.444  94.44   3       3       1       1/1     1       1       144     1       11058   11198
#Something is wrong with fields [0] and [3] here; they are identical for all retrieved reads. Check the main reads database for underscores and parse_seqids.






while ($filename = readdir(DIR)) {
  if ($filename =~ m/$searchstringa/) {
    if ($filename =~ m/$searchstringb$/) {
      $currfile = $dir."/".$filename;
      @stats = stat($currfile);
      $times{$currfile} = $stats[9];
    }
    elsif ($filename =~ m/$searchstringc$/) {
      $currfile = $dir."/".$filename;
      @stats = stat($currfile);
      $alttimes{$currfile} = $stats[9];
    }
    elsif ($filename =~ m/$searchstringd/ && $filename =~ m/txt$/) {
      $currfile = $dir."/".$filename;
      @stats = stat($currfile);
      $rettimes{$currfile} = $stats[9];
    }
  }
}
for $key (keys(%times)) {
  for $index (keys(%times)) {
    if ($index eq $key) {next;}
    elsif ($times{$key} == $times{$index}) {print OUT "$key and $index share time $times{$key}\n";}
  }
}
for $key (keys(%alttimes)) {
  for $index (keys(%alttimes)) {
    if ($index eq $key) {next;}
    elsif ($alttimes{$key} == $alttimes{$index}) {print OUT "$key and $index share time $alttimes{$key}\n";}
  }
}
for $key (keys(%rettimes)) {
  for $index (keys(%rettimes)) {
    if ($index eq $key) {next;}
    elsif ($rettimes{$key} == $rettimes{$index}) {print OUT "$key and $index share time $rettimes{$key}\n";}
  }
}
$nkeys = 0; $nsisters = 0; $nbrothers = 0;
for $key (keys(%times)) {
  $nkeys++;
  ($trialfile = $key) =~ s/$searchstringb/$searchstringc/;
  if (exists($alttimes{$trialfile})) {$sister{$key} = $trialfile; $nsisters++;}
  #TraesCS1D01G398900.1_cap3_increment_intactextracted2.fasta.cap_0_1_2.subset.contigs
  #TraesCS1D01G398900.1_cap3_increment_intact_vs_TraesCS1D01G398900.1retrieved_2_1e-20.txt
  $brother{$key} = "NULL";
  for $index (keys(%rettimes)) {
    if ($rettimes{$index} == $times{$key}) {$brother{$key} = $index; $nbrothers++;}
  }
}
print OUT "There are $nkeys files that match $searchstringb and $nsisters sister files that match $searchstringc\n";
print OUT "There are $nbrothers files that match $searchstringd\n";
if ($nkeys != $nsisters) {
  die "The file naming convention does not fit this script, which requires otherwise identical names with different suffixes, o some files are missing.\n";
}
if ($nkeys != $nbrothers) {
  print "Brother files were not necessarily found\; nkeys = $nkeys nbrothers = $nbrothers\n";
}
for $key (sort {$times{$a} <=> $times{$b}} (keys(%times))) {
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
      ($currcontig = $line) =~ s/>//;
    }
    else {$seqlength += length($line);}
  }
  if ($seqlength > 0) {push(@{$seqlengths{$key}}, $seqlength);}
  close (INP);
  $altfile = $sister{$key};
  $altncontigs = 0;
  open (INP, "< $altfile");
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
for $key (sort(keys(%longestmatchingcontiglength))) {
  print OUT "$key longest matching target is $longestmatchingcontiglength{$key}\n";
  $lkeys++;
  $totlongestlengths += $longestmatchingcontiglength{$key};
}
$akeys = 0;
$totaltlongestlengths = 0;
for $key (keys(%longestaltcontiglength)) {
  $akeys++; $totaltlongestlengths += $longestaltcontiglength{$key};
}
$ikeys = 0;
$tothighestidentities = 0;
for $key (keys(%highestidentity)) {
  $ikeys++; $tothighestidentities += $highestidentity{$key};
}
if ($lkeys > 0) {$meanlongestlength = $totlongestlengths / $lkeys;}
else {$meanlongestlength = "undef";}
if ($ikeys > 0) {$meanhighestidentity = $tothighestidentities / $ikeys;}
else {$meanhighestidentity = "undef";}
if ($akeys > 0) {$meanaltlongestlength = $totaltlongestlengths / $akeys;}
else {$meanaltlongestlength = "undef";}
print OUT "Mean longest matching contig length for $searchstringa is $meanlongestlength\n";
print OUT "Mean longest any contig length for $searchstringa is $meanaltlongestlength\n";
print OUT "Mean highest identity for $searchstringa is $meanhighestidentity\n";
