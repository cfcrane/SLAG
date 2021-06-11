#!/usr/bin/env perl
use warnings;
#This script set up runs of localassembly1115.pl for the SLAG paper.  Some naming conventions have been hard-coded.
if (scalar(@ARGV) != 4) {die "This script has four command-line arguments.\n";}
open (SEQ, "< $ARGV[0]"); #iwgsc_refseqv1.0_homeologous_group_1_chosen_groups_1107.fasta
open (TPL, "< $ARGV[1]"); #/scratch/halstead/c/ccrane/slagtests/SRspadesincrementtemplate.cfg
$basedir = $ARGV[2]; #/scratch/halstead/c/ccrane/slagtests
$differentiator = $ARGV[3]; #SR, MR, LR, LR1, LR2
$name = "NULL";
$nseqs = 0;
while ($line = <SEQ>) {
  chomp $line;
  if ($line =~ m/>/) {
    #>TraesCS1A01G397600.1 exon from chr1A position 563186055 to 563185774
    $line =~ s/>//;
    @vars = split(/\s+/, $line);
    $name = $vars[0];
    $sequence{$name} = "";
    $nseqs++;
  }
  else {
    $sequence{$name} .= $line;
  }
}
$nkeys = 0;
for $key (keys(%sequence)) {$nkeys++;}
print "There are $nseqs sequences in $ARGV[0] and $nkeys keys in the sequence hash.\n";
for $key (keys(%sequence)) {
  $fastafile = $basedir."/".$key.".fasta";
  open (OUT, "> $fastafile");
  print OUT ">$key\n$sequence{$key}\n";
  close (OUT);
  seek(TPL, 0, 0);
  while ($line = <TPL>) {
    chomp $line;
    if ($line =~ m/\$assembler/) {
      #$assembler = "spades";
      $line =~ s/\;//; $line =~ s/\"//g;
      @sars = split(/\s+/, $line);
      $assembler = $sars[2];
    }
    elsif ($line =~ m/\$extractionoption/) {
      #$extractionoption = "increment";
      $line =~ s/\;//; $line =~ s/\"//g;
      @cars = split(/\s+/, $line);
      $extractionoption = $cars[2];
    }
    elsif ($line =~ m/\$longread/) {
      #$longread = 2;
      $line =~ s/\;//; $line =~ s/\"//g;
      @bars = split(/\s+/, $line);
      if ($bars[2] == 0) {$fragcode = "intact";}
      elsif ($bars[2] == 1) {$fragcode = "frag";}
      elsif ($bars[2] == 2) {$fragcode = "doublefrag";}
    }
  }
  $configfile = $basedir."/".$key."_".$differentiator."_".$assembler."_".$extractionoption."_".$fragcode.".cfg";
  print "`./localassembly1115.pl $configfile`\;\n";
  seek (TPL, 0, 0);
  open (CON, "> $configfile");
  while ($line = <TPL>) {
    if ($line =~ m/\$seqfile/) {print CON "\$seqfile = \"$basedir", "/$key", ".fasta\"\;\n"; next;}
    elsif ($line =~ m/MMMM/) {
      $line =~ s/MMMM/$key/g;
    }
    if ($line =~ m/NNNN/) {
      $replacement = $basedir."/".$key."_".$assembler."_".$extractionoption."_".$fragcode;
      $line =~ s/NNNN/$replacement/g;
      print CON $line;
    }
    elsif ($line =~ m/workstem/) {
      chomp $line;
      @pars = split(/\s+/, $line);
      $workstem = $basedir."/".$key."_".$assembler."_".$extractionoption."_".$fragcode."_workstem";
      print CON "$pars[0] = \"$workstem\"\;\n";
    }
    elsif ($line !~ m/\$seqfile/) {print CON $line;}
  }
  close(CON);
}



#file section
$seqfile = "NNNN.fasta";
$spadesoutstem = "NNNN_spades_assembly";
$spadesworkstem = "/scratch/halstead/c/ccrane/spadeswork";
#general section
$restartflag = 0;
$maxcycle = 30;
$extractionoption = "increment";
$extincrement = 30;
$longread = 0;
$pairedend = 1;
#blast section
$querytype = "nucleotide";
$blastdir = "/group/bioinfo/apps/apps/blast-2.7.1+/bin";
$forwardblastdb = "/scratch/halstead/c/ccrane/iwgsc_refseqv2.0_cutouts_fake_paired_end_reads_60x_R1_db";
$reverseblastdb = "/scratch/halstead/c/ccrane/iwgsc_refseqv2.0_cutouts_fake_paired_end_reads_60x_R2_db";
$nthreads = 1;
$evalue = 1e-20;
$secevalue = 1e-30;
$runalign = 10000;
$carryforwardevalue = 1e-40;
$stem = "NNNN";
$tempdbname = "NNNN_tempdb";
$tempdboutstem = "NNNN_vs_temp_retrieved";
#assembler section
$assembler = "spades";
$spadesexe = "/group/bioinfo/apps/apps/spades-3.13.0/bin/spades.py";
$phredoffset = 33;
