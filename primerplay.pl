#!/usr/bin/env perl
use warnings;
#This script deals with fraction of working primer pairs from the longest contig in a SLAG or aTRAM2 run.
#Find the longest contig.
#Find primers.
#Align primers to reference genome and check expected amplicon size.
#Tally and present fraction of primer pairs that would work in reference genome.
do $ARGV[0];
open (LIS, "< $listfile");
open (OUT, "> $mainoutfile");
$ulimit = $maxampliconlength + $tol;
$llimit = $minampliconlength - $tol;
if ($llimit < 100) {die "Bad tolerance and minimum amplicon size: tol = $tol and minampliconlength = $minampliconlength\n";}
while ($file = <LIS>) {
  chomp $file;
  open (SEQ, "< $file");
  $currseq = "";
  $maxlength = 0;
  while ($line = <SEQ>) {
    chomp $line;
    if ($line =~ m/>/) {
      if (length($currseq) > $maxlength) {
        $maxlength = length($currseq);
        ($saveddefline = $currdefline) =~ s/>//g;
        $savedseq = $currseq;
      }
      $currdefline = $line;
      $currseq = "";
    }
    else {$currseq .= $line;}
  }
}
if ($maxlength < $minampliconlength) {
  print OUT "Longest contig length, $maxlength, is less than needed minimum amplicon length, $minampliconlength.\n";
  die "Insufficient contig length to produce amplicons at least $minampliconlength long.\n";
}
open (BLD, "> $boulderfile");
print BLD "SEQUENCE_ID=$enzymename", "_$saveddefline\n";
print BLD "SEQUENCE_TEMPLATE=$savedseq\n";
print BLD "PRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_RIGHT_PRIMER=1\n";
print BLD "PRIMER_MIN_SIZE=$minprimerlength\nPRIMER_MAX_SIZE=$maxprimerlength\nPRIMER_MAX_NS_ACCEPTED=0\n";
if (length($savedseq) > $maxampliconlength) {
  print BLD "PRIMER_PRODUCT_SIZE_RANGE=$minampliconlength-$maxampliconlength\nPRIMER_NUM_RETURN=$numreturn\n";
}
else {
  $vs = length($savedseq);
  print BLD "PRIMER_PRODUCT_SIZE_RANGE=$minampliconlength-$vs\nPRIMER_NUM_RETURN=$numreturn\n";
}
print BLD "PRIMER_EXPLAIN_FLAG=1\n=\n";
close (BLD);
`$primer3exe -output=$primer3outfile -error=$primer3errfile $boulderfile`;
open (POU, "< $primer3outfile");
while ($line = <POU>) {
  chomp $line;
  if ($line =~ m/PRIMER_LEFT/ && $line =~ m/_SEQUENCE/) {
    ($filler, $leftprimer) = split(/=/, $line);
    push (@leftprimers, $leftprimer);
  }
  if ($line =~ m/PRIMER_RIGHT/ && $line =~ m/_SEQUENCE/) {
    ($filler, $rightprimer) = split(/=/, $line);
    push (@rightprimers, $rightprimer);
  }
}
$vl = scalar(@leftprimers); $vr = scalar(@rightprimers);
if ($vl != $vr) {print OUT "We have unequal numbers of left and right primers. Left = $vl and rigth = $vr\n";}
if ($vl != $numreturn) {print OUT "Number of primer pairs = $vl but number requested is $numreturn.\n";}
open (QUE, "> $blastinfile");
for ($i = 0; $i < $vl; $i++) {print QUE ">left_primer_$i\n$leftprimers[$i]\n>right_primer_$i\n$rightprimers[$i]\n";}
close (QUE);
`$blastexe -task blastn-short -query $blastinfile -db $refgenomedb -out $blastoutfile -evalue $evalue -dust no -outfmt "6 qseqid qlen qstart qend sseqid sstart send pident sstrand"`;
open (BLA, "< $blastoutfile");
while ($line = <BLA>) {
  chomp $line;
  ($qname,$qlen,$qstart,$qend,$sname,$sstart,$send,$hitpct,$sstrand) = split(/\t/, $line);
#left_primer_11  20      1       20      ref|NC_024461.2|        6252632 6252651 100.000 plus
#right_primer_11 20      1       20      ref|NC_024463.2|        217596225       217596244       100.000 plus
#right_primer_11 20      3       20      ref|NC_024463.2|        160086227       160086210       100.000 minus
#right_primer_11 20      1       20      ref|NC_024461.2|        6254351 6254332 100.000 minus
  @pars = split(/_/, $qname);
  $basename = $pars[1]."_".$pars[2];
  if ($qname =~ m/left/) {
    #This will count but overwrite multiple hits on the same chromosome.
    if ($qend - $qstart + 1 == $qlen && $hitpct == 100) {
      push(@{$leftstarts{$basename}{$sname}}, $sstart);
      push(@{$leftends{$basename}{$sname}}, $send);
      push(@{$leftoriens{$basename}{$sname}}, $sstrand);
      if (exists($lefthits{$basename})) {$lefthits{$basename}++;}
      else {$lefthits{$basename} = 1;}
    }
  }
  if ($qname =~ m/right/) {
    if ($qend - $qstart + 1 == $qlen && $hitpct == 100) {
      push(@{$rightstarts{$basename}{$sname}}, $sstart);
      push(@{$rightends{$basename}{$sname}}, $send);
      push(@{$rightoriens{$basename}{$sname}}, $sstrand);
      if (exists($righthits{$basename})) {$righthits{$basename}++;}
      else {$righthits{$basename} = 1;}
    }
  }
}
for $key (keys(%leftstarts)) {
  $k = 0;
  for $index (keys(%{$leftstarts{$key}})) {$k++;}
  if ($lefthits{$key} != $k) {
    print OUT "$key has more than one left hit on at least one chromosome.\n";
  }
}
for $key (keys(%rightstarts)) {
  $k = 0;
  for $index (keys(%{$rightstarts{$key}})) {$k++;}
  if ($righthits{$key} != $k) {print OUT "$key has more than one right hit on at least one chromosome.\n";}
}
#for $key (sort(keys(%leftstarts))) {
#  for $index (sort(keys(%{$leftstarts{$key}}))) {
#    for $i (0 .. $#{$leftstarts{$key}{$index}}) {
#      print OUT "left $key $index $i $leftstarts{$key}{$index}[$i] $leftends{$key}{$index}[$i] $leftoriens{$key}{$index}[$i]\n";
#    }
#  }
#}
#for $key (sort(keys(%rightstarts))) {
#  for $index (sort(keys(%{$rightstarts{$key}}))) {
#    for $i (0 .. $#{$rightstarts{$key}{$index}}) {
#      print OUT "right $key $index $i $rightstarts{$key}{$index}[$i] $rightends{$key}{$index}[$i] $rightoriens{$key}{$index}[$i]\n";
#    }
#  }
#}
$workingpairs = 0;
$nloci = 0;
for $key (sort(keys(%leftstarts))) { #key is a primer name, like xx
  $foundforpair = 0;
  #print OUT "key = $key\n";
  for $index (sort(keys(%{$leftstarts{$key}}))) { #index is a chromosome
    #print OUT "	index = $index\n";
    for $i (0 .. $#{$leftstarts{$key}{$index}}) {
      #print OUT "		i = $i leftstarts{$key}{$index}[$i] = $leftstarts{$key}{$index}[$i]\n";
      $found = 0;
      if (exists($rightstarts{$key}{$index})) {
        for $j (0 .. $#{$rightstarts{$key}{$index}}) {
          if ($leftoriens{$key}{$index}[$i] ne $rightoriens{$key}{$index}[$j]) {
            $maxdist = abs($rightends{$key}{$index}[$j] - $leftstarts{$key}{$index}[$i]);
            $dist = abs($rightends{$key}{$index}[$j] - $leftends{$key}{$index}[$i]);
            if ($dist > $maxdist) {$maxdist = $dist;}
            $dist = abs($rightstarts{$key}{$index}[$j] - $leftstarts{$key}{$index}[$i]);
            if ($dist > $maxdist) {$maxdist = $dist;}
            $dist = abs($rightstarts{$key}{$index}[$j] - $leftends{$key}{$index}[$i]);
            if ($dist > $maxdist) {$maxdist = $dist;}
            if ($maxdist <= $ulimit && $maxdist >= $llimit) {$found++; $savedj = $j;}
          }
        }
      }
      if ($found > 0) {
        $nloci++; $foundforpair = 1; 
        print OUT "found = $found for $key $index i = $i j = $savedj maxdist = $maxdist\n";
      }
    }
  }
  if ($foundforpair > 0) {$workingpairs++;}
}
print OUT "There are $workingpairs potentially working primer pairs out of $vl pairs generated.\n";
print OUT "Number of potentially distinct amplicons is $nloci\n";
