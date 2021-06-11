#!/usr/bin/env perl
use warnings;
#This script collects count, length of longest contig, and percentage matching of a longest contig to a reference genome.
open (LIS, "< $ARGV[0]"); #Stanleyfullslaghexokinasecontigfiles.txt
open (OUT, "> $ARGV[1]");
$blastexe = $ARGV[2]; #/group/bioinfo/apps/apps/blast-2.10.0+/bin/blastn
$blastdb = $ARGV[3];
$evalue = $ARGV[4]; #1e-30
$longestlength = 0;
$greatestmeanlength = 0;
$greatestcount = 0;
$savedfile = "NULL";
while ($file = <LIS>) {
  chomp $file;
  open(INP, "< $file");
  $sequence = "";
  $count = 0;
  $totallength = 0;
  while ($line = <INP>) {
    chomp $line;
    if ($line =~ m/>/) {
      $totallength += length($sequence);
      if (length($sequence) > $longestlength) {
        $saveddefline = $currdefline;
        $longestlength = length($sequence);
        $savedfileforlength = $file;
        $savedsequence = $sequence;
      }
      $count++;
      $sequence = "";
    }
    else {$sequence .= $line;}
  }
  #print OUT "diag: file = $file longestlength = $longestlength contig count = $contigcount{$file}\n";
  $count++;
  $totallength += length($sequence);
  if (length($sequence) > $longestlength) { #Pick up a straggler.
    $saveddefline = $currdefline;
    $longestlength = length($sequence);
    $savedsequence = $sequence;
    $savedfileforlength = $file;
  }
  $meanlength = $totallength / $count;
  if ($count > $greatestcount) {$greatestcount = $count;}
  if ($meanlength > $greatestmeanlength) {$greatestmeanlength = $meanlength;}
  #Longest length has been handled as individual contigs are processed.
}
#
#$maxround = 0;
#for $key (keys(%meanlength)) {
#  #hemiStanleyhexokinasesSRcap3extracted20.fasta.cap_20.subset.contigs
#  @sars = split(/\./, $key);
#  ($parta, $partb) = split(/_/, $sars[-3]);
#  if ($partb > $maxround) {
#    $maxround = $partb;
#    $savedkey = $key;
#  }
#  if ($contigcount{$key} > 0) {$meanlength{$key} /= $contigcount{$key};}
#  else {print "$key produced no contigs.\n"; $meanlength{$key} = 0;}
#  #print OUT "diag: meanlength key = $key mean length = $meanlength{$key}\n";
#}
#$maxcount = 0;
#for $key (keys(%contigcount)) {
#  #print OUT "diag: key = $key contigcount{$key} = $contigcount{$key} maxcount = $maxcount\n";
#  if ($contigcount{$key} > $maxcount) {
#    $maxcount = $contigcount{$key};
#    $savedfileforcount = $key;
#  }
#}
#
@bars = split(/\./, $savedfileforlength);
$blastoutfile = $bars[0];
for ($i = 1; $i < scalar(@bars) - 1; $i++) {$blastoutfile .= ".".$bars[$i];}
$blastoutfile .= ".accblastout.txt";
`$blastexe -query $savedfileforlength -db $blastdb -out $blastoutfile -evalue $evalue -outfmt "6 qseqid qlen qstart qend sseqid sstart send nident length pident sstrand"`;
open (BLA, "< $blastoutfile");
$lastquery = "NULL";
$sum7 = 0;
$sum8 = 0;
while ($line = <BLA>) {
  chomp $line;
  @vars = split(/\t/, $line);
  if ($vars[0] ne $lastquery) {
    $sum7 += $vars[7];
    $sum8 += $vars[8];
    #$match = $vars[7] / $vars[8];
    # "diag: vars[7] = $vars[7] vars[8] = $vars[8] pident = $vars[9]\n";
  }
  $lastquery = $vars[0];
}
if ($sum8 > 0) {$match = $sum7 / $sum8;}
else {$match = "undefined";}
print OUT "maximum number of contigs = $greatestcount\n";
print OUT "maximum mean length of contigs = $greatestmeanlength\n";
print OUT "longest contig length = $longestlength in $savedfileforlength\n";
print OUT "sum of vars[7] / sum of vars[8] = $match\n";
