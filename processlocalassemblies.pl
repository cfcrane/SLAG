#!/usr/bin/env perl
use warnings;
#This script filled out the SRspades, SRcap3, and SRphrap lines for Zea in masterSLAGtable.docx.
open (LIA, "< $ARGV[0]"); #ZeaSRcap3basecontigsfiles.txt
open (LIB, "< $ARGV[1]"); #ZeaSRcap3subsetcontigsfiles.txt
open (DIC, "< $ARGV[2]"); #Zeagroups.txt
open (OUT, "> $ARGV[3]");
$maxfilecount = $ARGV[4]; #21
#-rw-r--r-- 1 ccrane student 18086 Apr 25 20:34 ZeaSRspadescontigsfastafiles.txt
#-rw-r--r-- 1 ccrane student 19346 Apr 25 20:35 ZeaSRspadessubsetcontigfiles.txt
#-rw-r--r-- 1 ccrane student 18926 Apr 25 20:40 ZeaSRphrapfastacontigsfiles.txt
#-rw-r--r-- 1 ccrane student 20396 Apr 27 11:09 ZeaSRphrapsubsetcontigsfiles.txt
#-rw-r--r-- 1 ccrane student 20086 Apr 25 20:42 ZeaSRcap3basecontigsfiles.txt
#-rw-r--r-- 1 ccrane student 21556 Apr 27 11:07 ZeaSRcap3subsetcontigsfiles.txt
while ($line = <DIC>) {chomp $line; push(@dictionary, $line);}
$nmaxfiles = 0;
$grandmeanbasecontigcount = 0;
$grandmeansubsetcontigcount = 0;
$grandmeanlongestbasecontiglength = 0;
$grandmeanlongestsubsetcontiglength = 0;
$grandmeanbasecontiglength = 0;
$grandmeansubsetcontiglength = 0;
$grandmeansetnbasecontigs = 0;
$grandmeannsubsetcontigs = 0;
$nproducedcontigs = 0;
for ($i = 0; $i < scalar(@dictionary); $i++) {
  seek (LIA, 0, 0);
  @parentfiles = ();
  while ($file = <LIA>) {
    if ($file =~ m/$dictionary[$i]/) {chomp $file; push(@parentfiles, $file);}
  }
  if (scalar(@parentfiles) == 0) {
    print OUT "$dictionary[$i] produced no contigs.\n";
    next;
  }
  $nproducedcontigs++;
  if (scalar(@parentfiles) == $maxfilecount) {$nmaxfiles++;}
  seek (LIB, 0, 0);
  @subsetcontigfiles = ();
  while ($file = <LIB>) {
    if ($file =~ m/$dictionary[$i]/) {chomp $file; push(@subsetcontigfiles, $file);}
  }
  if (scalar(@subsetcontigfiles) == 0) {next;} #Exclude a group whose initial assembly failed to produce a matching contig.
#  print OUT "term = $dictionary[$i]:\nparent files:\n";
#  for ($j = 0; $j < scalar(@parentfiles); $j++) {print OUT "$parentfiles[$j]\n";}
#  print OUT "\nsubset.contigs files:\n";
#  for ($j = 0; $j < scalar(@subsetcontigfiles); $j++) {print OUT "$subsetcontigfiles[$j]\n";}
  $longestbasecontiglength = 0;
  $longestbasecontigfile = "";
  $meanbasecontigcount = 0;
  for ($j = 0; $j < scalar(@parentfiles); $j++) {
    open (INP, "< $parentfiles[$j]");
    $currlength = 0;
    while ($line = <INP>) { #
      chomp $line;
      if ($line =~ />/) {
        if ($currlength > $longestbasecontiglength) {
          $longestbasecontiglength = $currlength;
          $longestbasecontigfile = $parentfiles[$j];
        }
        $currlength = 0;
        $meanbasecontigcount++;
      }
      else {$currlength += length($line);}
    }
    if ($currlength > $longestbasecontiglength) {
      $longestbasecontiglength = $currlength;
      $longestbasecontigfile = $parentfiles[$j];
    }
  }
  $meanbasecontigcount /= scalar(@parentfiles);
  $grandmeanbasecontigcount += $meanbasecontigcount;
  $grandmeanlongestbasecontiglength += $longestbasecontiglength;
  push (@longestbasecontiglengths, $longestbasecontiglength);
  print OUT "term = $dictionary[$i] mean count = $meanbasecontigcount longest base contig = $longestbasecontiglength in file $longestbasecontigfile\n";
  $longestsubsetcontiglength = 0;
  $longestsubsetcontigfile = "";
  $meansubsetcontigcount = 0;
  for ($j = 0; $j < scalar(@subsetcontigfiles); $j++) {
    open (INP, "< $subsetcontigfiles[$j]");
    $currlength = 0;
    while ($line = <INP>) {
      chomp $line;
      if ($line =~ m/>/) {
        if ($currlength > $longestsubsetcontiglength) {
          $longestsubsetcontiglength = $currlength;
          $longestsubsetcontigfile = $subsetcontigfiles[$j];
        }
        $currlength = 0;
        $meansubsetcontigcount++;
      }
      else {$currlength += length($line);}
    }
    if ($currlength > $longestsubsetcontiglength) {
      $longestsubsetcontiglength = $currlength;
      $longestsubsetcontigfile = $subsetcontigfiles[$j];
    }
  }
  $meansubsetcontigcount /= scalar(@subsetcontigfiles);
  $grandmeansubsetcontigcount += $meansubsetcontigcount;
  $grandmeanlongestsubsetcontiglength += $longestsubsetcontiglength;
  push (@longestsubsetcontiglengths, $longestsubsetcontiglength);
  print OUT "term = $dictionary[$i] mean count = $meansubsetcontigcount longest matching contig = $longestsubsetcontiglength in file $longestsubsetcontigfile\n";
  open (INP, "< $longestsubsetcontigfile");
  $meansubsetcontiglength = 0;
  $nsubsetcontigs = 0;
  $currlength = 0;
  while ($line = <INP>) {
    chomp $line;
    if ($line =~ m/>/) {
      $meansubsetcontiglength += $currlength;
      $nsubsetcontigs++;
      $currlength = 0;
    }
    else {$currlength += length($line);}
  }
  $meansubsetcontiglength += $currlength;
  $meansubsetcontiglength /= $nsubsetcontigs;
  $grandmeansubsetcontiglength += $meansubsetcontiglength;
  push (@meansubsetcontiglengths, $meansubsetcontiglength);
  push (@setnsubsetcontigs, $nsubsetcontigs);
  $grandmeannsubsetcontigs += $nsubsetcontigs;
  close(INP);
  open (INP, "< $longestbasecontigfile");
  $meanbasecontiglength = 0;
  $nbasecontigs = 0;
  $currlength = 0;
  while ($line = <INP>) {
    chomp $line;
    if ($line =~ m/>/) {
      $meanbasecontiglength += $currlength; #Precede this with $currlength = 0 before the while loop!
      $nbasecontigs++;
      $currlength = 0;
    }
    else {$currlength += length($line);}
  }
  $meanbasecontiglength += $currlength;
  $meanbasecontiglength /= $nbasecontigs;
  $grandmeanbasecontiglength += $meanbasecontiglength;
  push (@meanbasecontiglengths, $meanbasecontiglength);
  push (@setnbasecontigs, $nbasecontigs);
  $grandmeansetnbasecontigs += $nbasecontigs;
}
print OUT "Number of groups that produced at least one contig = $nproducedcontigs\n";
print OUT "Number of runs that went all $maxfilecount cycles = $nmaxfiles\n";
$grandmeanbasecontigcount /= $nproducedcontigs;
$grandmeansubsetcontigcount /= $nproducedcontigs;
$grandmeanlongestbasecontiglength /= $nproducedcontigs;
$grandmeanlongestsubsetcontiglength /= $nproducedcontigs;
$grandmeansubsetcontiglength /= $nproducedcontigs;
$grandmeanbasecontiglength /= $nproducedcontigs;
$grandmeansetnbasecontigs /= $nproducedcontigs;
$grandmeannsubsetcontigs /= $nproducedcontigs;
print OUT "For files that contain longest contigs, mean subset contig length = $grandmeansubsetcontiglength\n";
print OUT "For files that contain longest contigs, mean base contig count = $grandmeansetnbasecontigs\n";
print OUT "For files that contain longest contigs, mean base contig length = $grandmeanbasecontiglength\n";
print OUT "For files that contain longest contigs, mean subset contig count = $grandmeannsubsetcontigs\n";
#Put in code to print out grand mean contig counts for base and subset files.
print OUT "Mean matching contig count = $grandmeansubsetcontigcount overall mean base contig count = $grandmeanbasecontigcount\n";
print OUT "Mean longest matching contig length = $grandmeanlongestsubsetcontiglength\n";
print OUT "Mean longest base contig length = $grandmeanlongestbasecontiglength\n";
$stddevlongestbasecontiglength = 0;
$stddevlongestsubsetcontiglength = 0;
for ($i = 0; $i < scalar(@longestbasecontiglengths); $i++) {
  $stddevlongestbasecontiglength += ($grandmeanlongestbasecontiglength - $longestbasecontiglengths[$i]) ** 2;
}
for ($i = 0; $i < scalar(@longestsubsetcontiglengths); $i++) {
  $stddevlongestsubsetcontiglength += ($grandmeanlongestsubsetcontiglength - $longestsubsetcontiglengths[$i]) ** 2;
}
$stddevlongestbasecontiglength = sqrt($stddevlongestbasecontiglength / (-1 + scalar(@longestbasecontiglengths)));
$stddevlongestsubsetcontiglength = sqrt($stddevlongestsubsetcontiglength / (-1 + scalar(@longestsubsetcontiglengths)));
$stderrlongestbasecontiglength = $stddevlongestbasecontiglength / sqrt(scalar(@longestbasecontiglengths));
$stderrlongestsubsetcontiglength = $stddevlongestsubsetcontiglength / sqrt(scalar(@longestsubsetcontiglengths));
print OUT "Standard error of mean of longest base contig length = $stderrlongestbasecontiglength\n";
print OUT "Standard error of mean of longest matching contig length = $stderrlongestsubsetcontiglength\n";
$stddevmeansubsetcontiglengths = 0;
$stddevsetnsubsetcontig = 0;
for ($i = 0; $i < scalar(@meansubsetcontiglengths); $i++) {
  $stddevmeansubsetcontiglengths += ($grandmeansubsetcontiglength - $meansubsetcontiglengths[$i]) ** 2;
  $stddevsetnsubsetcontig += ($grandmeannsubsetcontigs - $setnsubsetcontigs[$i]) ** 2;
}
$stddevmeansubsetcontiglengths = sqrt($stddevmeansubsetcontiglengths / (-1 + scalar(@meansubsetcontiglengths)));
$stddevsetnsubsetcontig = sqrt($stddevsetnsubsetcontig / (-1 + scalar(@setnsubsetcontigs)));
$stderrmeansubsetcontiglengths = $stddevmeansubsetcontiglengths / sqrt(scalar(@meansubsetcontiglengths));
$stderrsetnsubsetcontigs = $stddevsetnsubsetcontig / sqrt(scalar(@setnsubsetcontigs));
print OUT "Standard error of mean of subset contig lengths in files with longest contig = $stderrmeansubsetcontiglengths\n";
print OUT "Standard error of mean of subset contig counts in files with longest contig = $stderrsetnsubsetcontigs\n";
$stddevmeanbasecontiglength = 0;
$stddevsetnbasecontigs = 0;
for ($i = 0; $i < scalar(@meanbasecontiglengths); $i++) {
  $stddevmeanbasecontiglength += ($grandmeanbasecontiglength - $meanbasecontiglengths[$i]) ** 2;
  $stddevsetnbasecontigs += ($grandmeansetnbasecontigs - $setnbasecontigs[$i]) ** 2;
}
$stddevmeanbasecontiglength = sqrt($stddevmeanbasecontiglength / (-1 + scalar(@meanbasecontiglengths)));
$stddevsetnbasecontigs = sqrt($stddevsetnbasecontigs / (-1 + scalar(@setnbasecontigs)));
$stderrmeanbasecontiglength = $stddevmeanbasecontiglength / sqrt(scalar(@meanbasecontiglengths));
$stderrsetnbasecontigs = $stddevsetnbasecontigs / sqrt(scalar(@setnbasecontigs));
print OUT "Standard error of mean of base contig lengths in files with longest contig = $stderrmeanbasecontiglength\n";
print OUT "Standard error of mean of base contig counts in files with longest contig = $stderrsetnbasecontigs\n";
