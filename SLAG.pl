#!/usr/bin/perl
use warnings;
#This script is set up for blastn with tabular output.
use Text::Wrap;
$configfile = $ARGV[0];
$restartflag = 0;
$longread = 0; #Override in $configfile for phrap or cap3 on fragmented long reads.
$pairedend = 0; #Override in $configfile for paired-end reads.
$useracon = 0; #Override in $configfile to use Racon to polish contigs at each cycle.
$readnamedelimiter = '\.'; #Override in $configfile to split paired-end read names on another character or string.
if (-e $configfile) {do $configfile;}
else {die "$configfile was not found, we cannot proceed.\n";}
if (!-e $seqfile) {
  $seqfile =~ s/fasta$/fa/;
  if (!-e $seqfile) {
    $seqfile =~ s/fa$/fna/;
    if (!-e $seqfile) {
      die "$seqfile was not found.  It must end in .fasta, .fa, or .fna.  Please check its path and spelling.\n";
    }
  }
}
$found = 0;
if ($assembler eq "phrap") {$found = 1;}
elsif ($assembler eq "cap3") {$found = 1;}
elsif ($assembler eq "canu") {$found = 1;}
elsif ($assembler eq "shasta") {$found = 1;}
elsif ($assembler eq "unicycler") {
  $found = 1;
  #`module load bioinfo`;
  #`module load Unicycler`;
}
elsif ($assembler eq "spades") {$found = 1;}
elsif ($assembler eq "masurca") {$found = 1;}
if ($found == 0) {die "Assembler $assembler is not supported in this program.\n";}
if ($longread == 1) {
  if (!defined($fragseqfilestem)) {die "For fragmented long reads and phrap or cap3 you need to specify the fragseqfile stem name.\n";}
  if (!defined($chunksize)) {die "For fragmented long reads and phrap or cap3 you need to specify the chunksize.\n";}
}
if ($longread == 2) {
  if (!defined($fragseqfilestem)) {die "For fragmented long reads and phrap or cap3 you need to specify the fragseqfile stem name.\n";}
  if (!defined($chunksizea)) {die "For doubly fragmented long reads and phrap or cap3 you need to specify chunksizea and chunksizeb.\n";}
}
if ($assembler eq "phrap") {
  if (!defined($phrapexe)) {die "Please specify the phrap executable with any necessary path in the configuration file.\n";}
  if (!defined($phrapsettings)) {die "Please specify the phrap parameter settings as a quoted string (possibly empty) in the configuration file.\n";}
}
elsif ($assembler eq "cap3") {
  if (!defined($cap3exe)) {die "Please specify the cap3 executable with any necessary path in the configuration file.\n";}
  if (!defined($cap3options)) {die "Please specify the cap3 options as a quoted string (possibly empty) in the configuration file.\n";}
  if (defined($cap3options)) {
    if ($cap3options =~ m/" -x "/) {
      @tars = split(/\s+/, $cap3options);
      for ($i = 0; $i < scalar(@tars); $i++) {
        #If tars[$i+1] is another option flag, this will give an error from cap3.
        if ($tars[$i] eq "-x") {$particlebase = $tars[$i+1];}
      }
    }
    else {$particlebase = "cap"; $cap3options .= " -x cap";}
  }
  else {$particlebase = "cap"; $cap3options = "-x cap";}
}
elsif ($assembler eq "canu") {
  if (!defined($canuexe)) {die "Please specify the canu executable with any necessary path in the configuration file.\n";}
  if (!defined($canusettings)) {die "Please specify the canu command line options as a quoted string (possibly empty) in the configuration file.\n";}
}
elsif ($assembler eq "spades") {
  if (!defined($spadesexe)) {die "Please specify the spades executable with any necessary path in the configuration file.\n";}
}
elsif ($assembler eq "masurca") {
  if (!defined($masurcaexe)) {die "Please specify the masurca executable with any necessary path in the configuration file.\n";}
}
if ($restartflag > 0) {
  if (!defined($foundingfile)) {die "Restarts require a distinct founding file.\n";}
  if (!defined($cycle)) {die "Restarts require a cycle greater than one.\n";}
}
else {$cycle = 0;} #We are starting at the beginning.
$logfile = $stem."log.txt";
open (LOG, "> $logfile");
$lastnacc = 0;
if ($restartflag == 0) {$foundingfile = $seqfile;}
$cycle = 0;
while (1) { #We leave the loop somewhere below.
  $currhr = (localtime)[2];
  $currmin = (localtime)[1];
  $currsec = (localtime)[0];
  print LOG "Cycle $cycle started at $currhr", ":$currmin", ":$currsec\n";
  print "Cycle $cycle started at $currhr", ":$currmin", ":$currsec\n";
  if ($pairedend == 1) {
    $forwardoutfile = $stem."forwardblastout".$cycle.".txt";
    $reverseoutfile = $stem."reverseblastout".$cycle.".txt";
  }
  elsif ($pairedend == 0) {$outfile = $stem."blastout".$cycle.".txt";}
  elsif ($pairedend == 2) {
    $forwardoutfile = $stem."forwardblastout".$cycle.".txt";
    $reverseoutfile = $stem."reverseblastout".$cycle.".txt";
    $lroutfile = $stem."blastout".$cycle.".txt";
  }
  if ($pairedend == 1) {
    if ($cycle == 0) {
      if ($querytype =~ m/protein/) {
        $forwardoutfile = $stem."initialforwardblastout.txt";
        print LOG "$blastdir/tblastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/tblastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        $reverseoutfile = $stem."initialreverseblastout.txt";
        print LOG "$blastdir/tblastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/tblastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      }
      else { #for nucleotide seeds
        print LOG "$blastdir/blastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/blastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        print LOG "$blastdir/blastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/blastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        print "Blast is done for cycle $cycle.\n";
      }
    }
    else { #for subsequent cycles
      print LOG "$blastdir/blastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
      `$blastdir/blastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      print LOG "$blastdir/blastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
      `$blastdir/blastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      print "Blast is done for cycle $cycle.\n";
    }
  }
  elsif ($pairedend == 0) { #Reads are singles.
    if ($cycle == 0) {
      if ($querytype =~ m/protein/) {
        $outfile = $stem."initialblastout.txt";
        #ppsomaticembryogenesisprotein_1209bcextracted8.fasta.Contig96   GINDECP02D4RAY  3.6e-64 1       245.66  1597    359     337     337     22      93.87   93.87   0       0       0       0       +1      1       359     +1      331     689
        print LOG "$blastdir/tblastn -db $blastdb -query $seqfile -out $outfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/tblastn -db $blastdb -query $seqfile -out $outfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        print "Blast is done for cycle $cycle.\n";
      }
      else {
        print LOG "$blastdir/blastn -db $blastdb -query $seqfile -out $outfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/blastn -db $blastdb -query $seqfile -out $outfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        print "Blast is done for cycle $cycle.\n";
      }
    }
    else { #This is for all subsequent cycles.
      print LOG "$blastdir/blastn -db $blastdb -query $seqfile -out $outfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
      `$blastdir/blastn -db $blastdb -query $seqfile -out $outfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      #`$blastdir/blastn -db $blastdb -query $seqfile -out $outfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt "6 0:qseqid 1:sseqid 2:evalue 3:qacc 4:bitscore 5:score 6:length 7:nident 8:positive 9:mismatch 10:pident 11:ppos 12:gaps 13:gapopen 14:sacc 15:frames 16:qframe 17:qstart 18:qend 19:sframe 20:sstart 21:send"`;
      print "Blast is done for cycle $cycle.\n";
    }
  }
  elsif ($pairedend == 2) {
    #Third option for MaSuRCA or for Racon polishing.
    if ($cycle == 0) {
      if ($querytype =~ m/protein/) {
        $forwardoutfile = $stem."initialforwardblastout.txt";
        print LOG "$blastdir/tblastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/tblastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        $reverseoutfile = $stem."initialreverseblastout.txt";
        print LOG "$blastdir/tblastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/tblastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        $lroutfile = $stem."initiallrblastout.txt";
        #ppsomaticembryogenesisprotein_1209bcextracted8.fasta.Contig96   GINDECP02D4RAY  3.6e-64 1       245.66  1597    359     337     337     22      93.87   93.87   0       0       0       0       +1      1       359     +1      331     689
        print LOG "$blastdir/tblastn -db $blastdb -query $seqfile -out $lroutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/tblastn -db $blastdb -query $seqfile -out $lroutfile -num_threads $nthreads -evalue $evalue -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      }
      else { #for nucleotide seeds
        print LOG "$blastdir/blastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/blastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        print LOG "$blastdir/blastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/blastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
        $lroutfile = $stem."initiallrblastout.txt";
        print LOG "$blastdir/blastn -db $blastdb -query $seqfile -out $lroutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
        `$blastdir/blastn -db $blastdb -query $seqfile -out $lroutfile -num_threads $nthreads -evalue $evalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      }
    }
    else { #for subsequent cycles
      print LOG "$blastdir/blastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
      `$blastdir/blastn -db $forwardblastdb -query $seqfile -out $forwardoutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      print LOG "$blastdir/blastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
      `$blastdir/blastn -db $reverseblastdb -query $seqfile -out $reverseoutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      print LOG "$blastdir/blastn -db $blastdb -query $seqfile -out $lroutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
      `$blastdir/blastn -db $blastdb -query $seqfile -out $lroutfile -num_threads $nthreads -evalue $secevalue -dust no -max_target_seqs $runalign -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
      print "Blast is done for cycle $cycle.\n";
    }
  }
  if ($pairedend == 1) {
    @blastats = stat($forwardoutfile);
    $hitstoblastouts = $blastats[7];
    @blastats = stat($reverseoutfile);
    $hitstoblastouts += $blastats[7];
    if ($hitstoblastouts == 0) {
      die "There were no blast hits to the founding sequence in $forwardoutfile or $reverseoutfile at cycle $cycle.\n";
    }
    @forwardaccarray = ();
    %forwardaccbitscores = ();
    $lastacc = "NULL";
    open (BLA, "< $forwardoutfile");
    while ($line = <BLA>) {
      chomp $line;
      @vars = split(/\t/, $line);
      if ($vars[7] + $vars[9] == 0) {next;} #Skip exceptional zero-valued lines in the blast output.
      if ($vars[1] ne $lastacc) {push(@forwardaccarray, $vars[1]); $forwardaccbitscores{$vars[1]} = $vars[4];}
      $lastacc = $vars[1];
    }
    close (BLA);
    @reverseaccarray = ();
    %reversebitscores = ();
    $lastacc = "NULL";
    open (BLA, "< $reverseoutfile");
    while ($line = <BLA>) {
      chomp $line;
      @vars = split(/\t/, $line);
      if ($vars[7] + $vars[9] == 0) {next;} #Skip exceptional zero-valued lines in the blast output.
      if ($vars[1] ne $lastacc) {push(@reverseaccarray, $vars[1]); $reverseaccbitscores{$vars[1]} = $vars[4];}
      $lastacc = $vars[1];
    }
    close (BLA);
    @combinedaccarray = ();
    %combinedaccbitscores = ();
    for $key (keys(%forwardaccbitscores)) {
      if (exists($reverseaccbitscores{$key})) {
        if ($forwardaccbitscores{$key} > $reverseaccbitscores{$key}) {$combinedaccbitscores{$key} = $forwardaccbitscores{$key};}
        else {$combinedaccbitscores{$key} = $reverseaccbitscores{$key};}
      }
      else {$combinedaccbitscores{$key} = $forwardaccbitscores{$key};}
    }
    for $key (keys(%reverseaccbitscores)) {
      if (!exists($forwardaccbitscores{$key})) {$combinedaccbitscores{$key} = $reverseaccbitscores{$key};}
    }
    for $key (sort {$combinedaccbitscores{$b} <=> $combinedaccbitscores{$a}} (keys(%combinedaccbitscores))) {
      push(@accarray, $key);
      push(@accbitscores, $combinedaccbitscores{$key});
    }
  }
  elsif ($pairedend == 0) { #Reads are singles.
    @blastats = stat($outfile);
    if ($blastats[7] == 0) {die "There were no blast hits to the founding sequence in $outfile at cycle $cycle.\n";}
    @accarray = ();
    @accbitscores = ();
    $lastacc = "NULL";
    open (BLA, "< $outfile");
    while ($line = <BLA>) {
      chomp $line;
      @vars = split(/\t/, $line);
      if ($vars[7] + $vars[9] == 0) {next;} #Skip exceptional zero-valued lines in the blast output.
      if ($vars[1] ne $lastacc) {push(@accarray, $vars[1]); push(@accbitscores, $vars[4]);}
      $lastacc = $vars[1];
    }
    close (BLA);
  }
  elsif ($pairedend == 2) { #Get both long and short reads for MaSuRCA.
    @blastats = stat($forwardoutfile);
    $hitstoblastouts = $blastats[7];
    @blastats = stat($reverseoutfile);
    $hitstoblastouts += $blastats[7];
    if ($hitstoblastouts == 0) {
      die "There were no blast hits to the founding sequence in $forwardoutfile or $reverseoutfile at cycle $cycle.\n";
    }
    @forwardaccarray = ();
    %forwardaccbitscores = ();
    $lastacc = "NULL";
    open (BLA, "< $forwardoutfile");
    while ($line = <BLA>) {
      chomp $line;
      @vars = split(/\t/, $line);
      if ($vars[7] + $vars[9] == 0) {next;} #Skip exceptional zero-valued lines in the blast output.
      if ($vars[1] ne $lastacc) {push(@forwardaccarray, $vars[1]); $forwardaccbitscores{$vars[1]} = $vars[4];}
      $lastacc = $vars[1];
    }
    close (BLA);
    @reverseaccarray = ();
    %reversebitscores = ();
    $lastacc = "NULL";
    open (BLA, "< $reverseoutfile");
    while ($line = <BLA>) {
      chomp $line;
      @vars = split(/\t/, $line);
      if ($vars[7] + $vars[9] == 0) {next;} #Skip exceptional zero-valued lines in the blast output.
      if ($vars[1] ne $lastacc) {push(@reverseaccarray, $vars[1]); $reverseaccbitscores{$vars[1]} = $vars[4];}
      $lastacc = $vars[1];
    }
    close (BLA);
    @combinedaccarray = ();
    %combinedaccbitscores = ();
    for $key (keys(%forwardaccbitscores)) {
      if (exists($reverseaccbitscores{$key})) {
        if ($forwardaccbitscores{$key} > $reverseaccbitscores{$key}) {$combinedaccbitscores{$key} = $forwardaccbitscores{$key};}
        else {$combinedaccbitscores{$key} = $reverseaccbitscores{$key};}
      }
      else {$combinedaccbitscores{$key} = $forwardaccbitscores{$key};}
    }
    for $key (keys(%reverseaccbitscores)) {
      if (!exists($forwardaccbitscores{$key})) {$combinedaccbitscores{$key} = $reverseaccbitscores{$key};}
    }
    for $key (sort {$combinedaccbitscores{$b} <=> $combinedaccbitscores{$a}} (keys(%combinedaccbitscores))) {
      push(@accarray, $key);
      push(@accbitscores, $combinedaccbitscores{$key});
    }
    @blastats = stat($lroutfile);
    if ($blastats[7] == 0) {die "There were no blast hits to the founding sequence in $lroutfile at cycle $cycle.\n";}
    @lraccarray = ();
    @lraccbitscores = ();
    $lastacc = "NULL";
    open (BLA, "< $lroutfile");
    while ($line = <BLA>) {
      chomp $line;
      @vars = split(/\t/, $line);
      if ($vars[7] + $vars[9] == 0) {next;} #Skip exceptional zero-valued lines in the blast output.
      if ($vars[1] ne $lastacc) {push(@lraccarray, $vars[1]); push(@lraccbitscores, $vars[4]);}
      $lastacc = $vars[1];
    }
    close (BLA);
  }
  if ($extractionoption eq "population") {
    if ($cycle == 0) {$currlimit = scalar(@accarray);}
    else {
      $nprevused = 0;
      for $key (keys(%previouslyused)) {$nprevused++; $prevusedfoundincurrent{$key} = 0;}
      $usedsofar = 0;
      for ($i = 0; $i < scalar(@accarray); $i++) {
        if (exists($previouslyused{$accarray[$i]})) {$usedsofar++; $prevusedfoundincurrent{$accarray[$i]} = 1;}
        if ($usedsofar == $nprevused) {$currlimit = $i; last;}
      }
      if ($usedsofar < $nprevused) {
        for $key (keys(%prevusedfoundincurrent)) {
          if ($prevusedfoundincurrent{$key} == 0) {print "$key did not appear in the reads of this cycle $cycle.\n";}
        }
        print "Not all previously used reads appear in this cycle $cycle .\n"; 
        print "Total previously used = $nprevused and $usedsofar of them have appeared now.\n";
        die "Consider increasing runalign from its current value of $runalign .\n";
      }
    }
  }
  elsif ($extractionoption eq "increment") {
    if ($cycle == 0) {$baselimit = scalar(@accarray);}
    $currlimit = $baselimit + $cycle * $extincrement;
  }
  elsif ($extractionoption eq "manual") { #@manarray is set directly in the configuration file.
    if ($cycle < scalar(@manarray)) {$currlimit = $manarray[$cycle];}
  }
  elsif ($extractionoption eq "all") {$currlimit = scalar(@accarray);}
  elsif ($extractionoption eq "bitscore") {
    $currlimit = scalar(@accbitscores);
    print LOG "DIAGNOSTICS: Preliminary value of currlimit = $currlimit\n";
    for ($j = 0; $j < 40; $j++) {print LOG "cycle = $cycle accbitscores[$j] = $accbitscores[$j]\n";}
    for ($i = 0; $i < scalar(@accbitscores); $i++) {
      if ($cycle > 0) { 
        if ($accbitscores[$i] < $minbitscore) {$currlimit = $i - 1; last;}
      }
    }
    $v = scalar(@accbitscores);
    print LOG "DIAGNOSTICS: outfile = $outfile currlimit = $currlimit in cycle $cycle , minimum bitscore = $minbitscore, scalar(accbitscores) = $v\n";
  }
  else {die "The value of the extraction option is not recognized, exiting.\n";}
  print LOG "The current number of reads to assemble in cycle $cycle is $currlimit.\n";
  $accessionfile = $stem."accessions".$cycle.".txt"; #The extracted accession names go here.
  open (ACC, "> $accessionfile");
  $nacc = 0;
  %previouslyused = ();
  for ($i = 0; $i < scalar(@accarray); $i++) { #@accarray is in descending order of bitscore.
    if (!exists($previouslyused{$accarray[$i]})) { #Also prevent duplicate names in ACC.
      print ACC "$accarray[$i]\n";
      $previouslyused{$accarray[$i]} = 1;
      $nacc++;
    }
    if ($nacc >= $currlimit) {last;}
  }
  close(ACC);
  print LOG "There were $nacc reads that were used in cycle $cycle.\n";
  if ($assembler eq "phrap" && $nacc > $maxnacc) { #Stop if phrap will be overwhelmed.
    print LOG "Hit count exceeded $maxnacc, phrap could be overwhelmed.\n";
    last;
  }
  if ($pairedend == 1) {
    $accessionfilef = $accessionfile."forwards.txt";
    $accessionfiler = $accessionfile."reverses.txt";
    open (ACS, "< $accessionfile");
    open (ACF, "> $accessionfilef");
    open (ACR, "> $accessionfiler");
    while ($line = <ACS>) {
      chomp $line;
      @mars = split(/$readnamedelimiter/, $line);
      $readname = $mars[0];
      for ($i = 1; $i < scalar(@mars) - 1; $i++) {
        if ($readnamedelimiter eq '\.') {$readname .= ".".$mars[$i];}
        else {$readname .= $readnamedelimiter.$mars[$i];}
      }
      $readnames{$readname} = 1;
    }
    close(ACS);
    for $key (keys(%readnames)) {
      if ($readnamedelimiter eq '\.') {
        print ACF $key, ".1\n";
        print ACR $key, ".2\n";
      }
      else {
        print ACF $key, $readnamedelimiter, "1\n";
        print ACR $key, $readnamedelimiter, "2\n";
      }
    }
    close (ACF); close (ACR);
    #Put forwards of each pair in filef, reverses in filer.
    $forwarddboutfile = $stem."_R1_extracted".$cycle.".fasta";
    print LOG "$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfilef -out $forwarddboutfile -dbtype nucl -db $forwardblastdb\n";
    `$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfilef -out $forwarddboutfile -dbtype nucl -db $forwardblastdb`;
    $reversedboutfile = $stem."_R2_extracted".$cycle.".fasta";
    print LOG "$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfiler -out $reversedboutfile -dbtype nucl -db $reverseblastdb\n";
    `$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfiler -out $reversedboutfile -dbtype nucl -db $reverseblastdb`;
    #Merge files to $dboutfile here; make sure that R1 and R2 reads have different names.
    #While $dboutfile is not needed by spades, it is needed by phrap.
    $dboutfile = $stem."extracted".$cycle.".fasta";
    open (DBO, "> $dboutfile");
    open (INP, "< $forwarddboutfile");
    while ($line = <INP>) {
      if ($line =~ m/>/) {
        if ($line !~ m/R1/) {
          chomp $line;
          @pars = split(/\s+/, $line);
          $pars[0] .= "_R1";
          print DBO $pars[0];
          for ($i = 1; $i < scalar(@pars); $i++) {print DBO " $pars[$i]";}
          print DBO "\n";
        }
        else {print DBO $line;}
      }
      else {print DBO $line;}
    }
    close (INP);
    open (INP, "< $reversedboutfile");
    while ($line = <INP>) {
      if ($line =~ m/>/) {
        if ($line !~ m/R2/) {
          chomp $line;
          @pars = split(/\s+/, $line);
          $pars[0] .= "_R2";
          print DBO $pars[0];
          for ($i = 1; $i < scalar(@pars); $i++) {print DBO " $pars[$i]";}
          print DBO "\n";
        }
        else {print DBO $line;}
      }
      else {print DBO $line;}
    }
    close(DBO);
  }
  elsif ($pairedend == 0) {
    $dboutfile = $stem."extracted".$cycle.".fasta"; #The extracted sequences go here.
    print LOG "$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfile -out $dboutfile -dbtype nucl -db $blastdb\n";
    print "$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfile -out $dboutfile -dbtype nucl -db $blastdb\n";
    `$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfile -out $dboutfile -dbtype nucl -db $blastdb`;
  }
  elsif ($pairedend == 2) {
    #third option for running MaSuRCA. Extract reads from two databases --> ONT and Illumina.
    $forwarddboutfile = $stem."_R1_extracted".$cycle.".fasta";
    print LOG "$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfile -out $forwarddboutfile -dbtype nucl -db $forwardblastdb\n";
    `$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfile -out $forwarddboutfile -dbtype nucl -db $forwardblastdb`;
    $reversedboutfile = $stem."_R2_extracted".$cycle.".fasta";
    print LOG "$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfile -out $reversedboutfile -dbtype nucl -db $reverseblastdb\n";
    `$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $accessionfile -out $reversedboutfile -dbtype nucl -db $reverseblastdb`;
    $lraccessionfile = $stem."lraccessions".$cycle.".txt"; #The extracted long-read accession names go here.
    open (LCC, "> $lraccessionfile");
    $lrnacc = 0;
    %previouslyused = ();
    for ($i = 0; $i < scalar(@lraccarray); $i++) { #@lraccarray is in descending order of bitscore.
      print LCC "$lraccarray[$i]\n";
      $previouslyused{$lraccarray[$i]} = 1;
      $lrnacc++;
      if ($lrnacc >= $currlimit) {last;}
    }
    close(LCC);
    $lrdboutfile = $stem."extracted".$cycle.".fasta"; #The extracted sequences go here.
    print LOG "$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $lraccessionfile -out $lrdboutfile -dbtype nucl -db $blastdb\n";
    print "$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $lraccessionfile -out $lrdboutfile -dbtype nucl -db $blastdb\n";
    `$blastdir/blastdbcmd -get_dups -outfmt %f -entry_batch $lraccessionfile -out $lrdboutfile -dbtype nucl -db $blastdb`;
  }
  print "Sequence retrieval is done for cycle $cycle.\n";
  if ($assembler eq "phrap") {
    if ($longread == 1) { #Fragment the reads in $dboutfile.
      open (SEC, "< $dboutfile");
      $fragseqfile = $fragseqfilestem.$cycle.".fasta";
      open (FRG, "> $fragseqfile");
      $sequence = "";
      $lastdefline = "NULL";
      while ($line = <SEC>) {
        chomp $line;
        if ($line =~ m/^#/) {next;}
        if ($line =~ m/>/) {
          $seqlength = length($sequence);
          if ($seqlength > 4) {
            $div = int($seqlength / $chunksize);
            $divlength = 1 + int($seqlength / $div);
            for ($i = 0; $i < $div; $i++) {
              $lastdefline =~ s/ /_/g;
              $piececount{$lastdefline} = $div;
              print FRG "$lastdefline.$i\n";
              $part = substr($sequence, $i * $divlength, $divlength);
              print FRG wrap('', '', $part), "\n";
            }
          }
          $lastdefline = $line;
          $sequence = "";
        }
        else {$sequence .= $line;}
      }
      $seqlength = length($sequence);
      if ($seqlength > 4) {
        $div = int($seqlength / $chunksize);
        $divlength = 1 + int($seqlength / $div);
        for ($i = 0; $i < $div; $i++) {
          $lastdefline =~ s/ /_/g;
          $piececount{$lastdefline} = $div;
          print FRG "$lastdefline.$i\n";
          $part = substr($sequence, $i * $divlength, $divlength);
          print FRG wrap('', '', $part), "\n";
        }
      }
      $dboutfile = $fragseqfile;
      close(FRG);
    }
    if ($longread == 2) { #Doubly fragment the retrieved reads.
      if (!defined($chunksizea) || !defined($chunksizeb)) {die "Values of chunksizea and chunksizeb need to be given in the configuration file.\n";}
      open (SEC, "< $dboutfile");
      $fragseqfile = $fragseqfilestem.$cycle.".fasta";
      open (FRG, "> $fragseqfile");
      $sequence = "";
      $lastdefline = "NULL";
      while ($line = <SEC>) {
        chomp $line;
        if ($line =~ m/^#/) {next;}
        if ($line =~ m/>/) {
          $seqlength = length($sequence);
          if ($seqlength > 4) {
            $div = int($seqlength / $chunksizea);
            $divlength = 1 + int($seqlength / $div);
            for ($i = 0; $i < $div; $i++) {
              $lastdefline =~ s/ /_/g;
              $piececount{$lastdefline} = $div;
              print FRG $lastdefline, "a.$i\n";
              $part = substr($sequence, $i * $divlength, $divlength);
              print FRG wrap('', '', $part), "\n";
            }
          }
          $lastdefline = $line;
          $sequence = "";
        }
        else {$sequence .= $line;}
      }
      $seqlength = length($sequence);
      if ($seqlength > 4) {
        $div = int($seqlength / $chunksizea);
        $divlength = 1 + int($seqlength / $div);
        for ($i = 0; $i < $div; $i++) {
          $lastdefline =~ s/ /_/g;
          $piececount{$lastdefline} = $div;
          print FRG $lastdefline, "a.$i\n";
          $part = substr($sequence, $i * $divlength, $divlength);
          print FRG wrap('', '', $part), "\n";
        }
      }
      seek(SEC, 0, 0);
      while ($line = <SEC>) {
        chomp $line;
        if ($line =~ m/^#/) {next;}
        if ($line =~ m/>/) {
          $seqlength = length($sequence);
          if ($seqlength > 4) {
            $div = int($seqlength / $chunksizeb);
            $divlength = 1 + int($seqlength / $div);
            for ($i = 0; $i < $div; $i++) {
              $lastdefline =~ s/ /_/g;
              $piececount{$lastdefline} = $div;
              print FRG $lastdefline, "b.$i\n";
              $part = substr($sequence, $i * $divlength, $divlength);
              print FRG wrap('', '', $part), "\n";
            }
          }
          $lastdefline = $line;
          $sequence = "";
        }
        else {$sequence .= $line;}
      }
      $seqlength = length($sequence);
      if ($seqlength > 4) {
        $div = int($seqlength / $chunksizeb);
        $divlength = 1 + int($seqlength / $div);
        for ($i = 0; $i < $div; $i++) {
          $lastdefline =~ s/ /_/g;
          $piececount{$lastdefline} = $div;
          print FRG $lastdefline, "b.$i\n";
          $part = substr($sequence, $i * $divlength, $divlength);
          print FRG wrap('', '', $part), "\n";
        }
      }
      $dboutfile = $fragseqfile;
      close(FRG);
    }
    $phrapoutfile = $stem."phrapped".$cycle.".txt";
    print LOG "$phrapexe $phrapsettings $dboutfile > $phrapoutfile\n";
    print "$phrapexe $phrapsettings $dboutfile > $phrapoutfile\n";
    `$phrapexe $phrapsettings $dboutfile > $phrapoutfile`;
    $tigfile = $dboutfile.".contigs";
    @tigfilestats = stat($tigfile);
    if ($tigfilestats[7] < 2) {die "Phrap produced no contigs in cycle $cycle\n";}
    $tggfile = $dboutfile.".subset.contigs";
    $singletfile = $dboutfile.".singlets";
    open (SNG, "< $singletfile");
    while ($line = <SNG>) {
      if ($line =~ m/>/) {
        for $key (keys(%previouslyused)) {
          #This is not perfect, but it does not presuppose any particular structure to the read names.
          #If n1114 is to be deleted, but n11149 is also present in the list of used accessions, then BOTH will be deleted.
          if ($line =~ m/$key/) {delete($previouslyused{$key});}
        }
      }
    }
    close (SNG);
  }
  elsif ($assembler eq "cap3") {
    if ($longread == 1) {
      #Apply fragmentcontigs.pl to $dboutfile.
      open (SEC, "< $dboutfile") || die "File $dboutfile could not be opened.\n";
      print "Sequences in file $dboutfile will be fragmented for cap3.\n";
      $fragseqfile = $fragseqfilestem.$cycle.".fasta";
      open (FRG, "> $fragseqfile");
      $sequence = "";
      $lastdefline = "NULL";
      while ($line = <SEC>) {
        chomp $line;
        if ($line =~ m/^#/) {next;}
        if ($line =~ m/>/) {
          $seqlength = length($sequence);
          if ($seqlength > 4) {
            $div = int($seqlength / $chunksize);
            $divlength = 1 + int($seqlength / $div);
            for ($i = 0; $i < $div; $i++) {
              $lastdefline =~ s/ /_/g;
              print FRG "$lastdefline.$i\n";
              $part = substr($sequence, $i * $divlength, $divlength);
              print FRG wrap('', '', $part), "\n";
            }
          }
          $lastdefline = $line;
          $sequence = "";
        }
        else {$sequence .= $line;}
      }
      $seqlength = length($sequence);
      if ($seqlength > 4) {
        $div = int($seqlength / $chunksize);
        $divlength = 1 + int($seqlength / $div);
        for ($i = 0; $i < $div; $i++) {
          $lastdefline =~ s/ /_/g;
          print FRG "$lastdefline.$i\n";
          $part = substr($sequence, $i * $divlength, $divlength);
          print FRG wrap('', '', $part), "\n";
        }
      }
      $dboutfile = $fragseqfile;
      close(FRG);
    }
    if ($longread == 2) {
      if (!defined($chunksizea) || !defined($chunksizeb)) {die "Values of chunksizea and chunksizeb need to be given in the configuration file.\n";}
      open (SEC, "< $dboutfile");
      $fragseqfile = $fragseqfilestem.$cycle.".fasta";
      open (FRG, "> $fragseqfile");
      $sequence = "";
      $lastdefline = "NULL";
      while ($line = <SEC>) {
        chomp $line;
        if ($line =~ m/^#/) {next;}
        if ($line =~ m/>/) {
          $seqlength = length($sequence);
          if ($seqlength > 4) {
            $div = int($seqlength / $chunksizea);
            $divlength = 1 + int($seqlength / $div);
            for ($i = 0; $i < $div; $i++) {
              $lastdefline =~ s/ /_/g;
              $piececount{$lastdefline} = $div;
              print FRG $lastdefline, "a.$i\n";
              $part = substr($sequence, $i * $divlength, $divlength);
              print FRG wrap('', '', $part), "\n";
            }
          }
          $lastdefline = $line;
          $sequence = "";
        }
        else {$sequence .= $line;}
      }
      $seqlength = length($sequence);
      if ($seqlength > 4) {
        $div = int($seqlength / $chunksizea);
        $divlength = 1 + int($seqlength / $div);
        for ($i = 0; $i < $div; $i++) {
          $lastdefline =~ s/ /_/g;
          $piececount{$lastdefline} = $div;
          print FRG $lastdefline, "a.$i\n";
          $part = substr($sequence, $i * $divlength, $divlength);
          print FRG wrap('', '', $part), "\n";
        }
      }
      seek(SEC, 0, 0);
      while ($line = <SEC>) {
        chomp $line;
        if ($line =~ m/^#/) {next;}
        if ($line =~ m/>/) {
          $seqlength = length($sequence);
          if ($seqlength > 4) {
            $div = int($seqlength / $chunksizeb);
            $divlength = 1 + int($seqlength / $div);
            for ($i = 0; $i < $div; $i++) {
              $lastdefline =~ s/ /_/g;
              $piececount{$lastdefline} = $div;
              print FRG $lastdefline, "b.$i\n";
              $part = substr($sequence, $i * $divlength, $divlength);
              print FRG wrap('', '', $part), "\n";
            }
          }
          $lastdefline = $line;
          $sequence = "";
        }
        else {$sequence .= $line;}
      }
      $seqlength = length($sequence);
      if ($seqlength > 4) {
        $div = int($seqlength / $chunksizeb);
        $divlength = 1 + int($seqlength / $div);
        for ($i = 0; $i < $div; $i++) {
          $lastdefline =~ s/ /_/g;
          $piececount{$lastdefline} = $div;
          print FRG $lastdefline, "b.$i\n";
          $part = substr($sequence, $i * $divlength, $divlength);
          print FRG wrap('', '', $part), "\n";
        }
      }
      $dboutfile = $fragseqfile;
      close(FRG);
    }
    $cap3outfile = $cap3outstem."_$cycle".".txt";
    $particle = $particlebase."_".$cycle;
    if (!defined($cap3options)) {$currcap3options = "-x $particle";}
    else {($currcap3options = $cap3options) =~ s/-x $particlebase/-x $particle/;}
    `$cap3exe $dboutfile $currcap3options > $cap3outfile`;
    print LOG "$cap3exe $dboutfile $currcap3options > $cap3outfile\n";
    print LOG "Value of particle at line 492 is $particle\n";
    $tigfile = $dboutfile.".".$particle.".contigs";
    @tigfilestats = stat($tigfile);
    if ($tigfilestats[7] < 2) {die "Cap3 produced no contigs in cycle $cycle.\n";}
    $tggfile = $dboutfile.".".$particle.".subset.contigs";
    print LOG "tigfile = $tigfile and tggfile = $tggfile\n";
  }
  elsif ($assembler eq "canu") {
    #canu -p lepbiglobcanuassemb_102972_700_overends -d /scratch/halstead/c/ccrane/lepcanu16 genomeSize=10k usegrid=false correctedErrorRate=0.105 -nanopore-raw lepbiglobreads_102972_103672_0_700.fastq
    $canuworkdir = $canuworkstem.$cycle;
    if (-e $canuworkdir) {die "Please empty and remove $canuworkdir before proceeding.\n";}
    ($currcanusettings = $canusettings) =~ s/canuworkdirectory/$canuworkdir/;
    $currcanusettings =~ s/estgenomesize/$lastnuccount/;
    print LOG "$canuexe $currcanusettings $dboutfile\n";
    `$canuexe $currcanusettings $dboutfile`;
    $tigfile = $canuworkdir."/".$canuprefix.".contigs.fasta";
    $tggfile = $canuworkdir."/".$canuprefix.$cycle.".subset.contigs";
    print LOG "From cycle $cycle, tigfile = $tigfile and tggfile = $tggfile .\n";
  }
  elsif ($assembler eq "shasta") {
    $moddboutfile = "mod".$dboutfile;
    open (DBT, "< $dboutfile");
    open (MOD, "> $moddboutfile");
    $modcount = 0;
    while ($line = <DBT>) {
      if ($line =~ m/>/) {
        if ($modcount == 0) {print MOD $line;}
        else {print MOD "\n$line";}
      }
      else {chomp $line; print MOD $line;}
      $modcount++;
    }
    close(DBT); close(MOD);
    $shastaworkdir = $shastaworkstem.$cycle;
    if (-e $shastaworkdir) {die "Please clean out and delete $shastaworkdir before proceeding.\n";}
    ($currshastasettings = $shastasettings) =~ s/shastaworkdirectory/$shastaworkdir/;
    #/group/bioinfo/apps/apps/shasta-0.3.0/Shasta --input $inputfastafile --assemblyDirectory $workdirectory --Reads.minReadLength $minreadlength --threads $nthreads
    print LOG "$shastaexe --input $moddboutfile --assemblyDirectory $shastaworkdir $currshastasettings\n";
    `$shastaexe --input $moddboutfile --assemblyDirectory $shastaworkdir $currshastasettings`;
    $tigfile = $shastaworkdir."/Assembly.fasta";
    if (!-e $tigfile) {die "Post-assembly contig file $tigfile does not exist.\n";}
    $tggfile = $shastaworkdir."/".$shastaoutstem.$cycle.".subset.contigs";
    print LOG "From cycle $cycle, tigfile = $tigfile and tggfile = $tggfile .\n";
  }
  elsif ($assembler eq "unicycler") {
    $unicyclerworkdir = $unicyclerworkstem.$cycle;
    if (-e $unicyclerworkdir) {die "Please clean out and delete $unicyclerworkdir before proceeding.\n";}
#    if ($dboutfile =~ m/.fasta$/) {($dboutfastq = $dboutfile) =~ s/.fasta/.fastq/;}
#    elsif ($dboutfile =~ m/.fna$/) {($dboutfastq = $dboutfile) =~ s/.fna/.fastq/;}
#    elsif ($dboutfile =~ m/.fas$/) {($dboutfastq = $dboutfile) =~ s/.fas/.fastq/;}
#    else {$dboutfastq = $dboutfile.".fastq";}
#    open (INP, "< $dboutfile");
#    open (FSQ, "> $dboutfastq");
#    $sequence = "";
#    while ($line = <INP>) {
#      chomp $line;
#      if ($line =~ m/>/) {
#        if (length($sequence) > 2) {push(@sequences, $sequence);}
#        $line =~ s/>/@/;
#        push (@deflines, $line);
#        $sequence = "";
#      }
#      else {$sequence .= $line;}
#    }
#    if (length($sequence) > 2) {push(@sequences, $sequence);}
#    for ($i = 0; $i < scalar(@deflines); $i++) {
#      print FSQ "$deflines[$i]\n$sequences[$i]\n+\n";
#      for ($j = 0; $j < length($sequences[$i]); $j++) {print FSQ "B";}
#      print FSQ "\n";
#    }
#    close (FSQ);
    print LOG "$unicyclerexe -l $dboutfile -o $unicyclerworkdir -t $nthreads\n";
    `$unicyclerexe -l $dboutfile -o $unicyclerworkdir -t $nthreads --linear_seqs 1 --no_rotate`;
    $tigfile = $unicyclerworkdir."/assembly.fasta";
    if (!-e $tigfile) {die "Post-assembly contig file $tigfile does not exist.\n";}
    $tggfile = $unicyclerworkdir."/".$unicycleroutstem.".subset.contigs";
    print LOG "From cycle $cycle, tigfile = $tigfile and tggfile = $tggfile .\n";
  }
  elsif ($assembler eq "soap") {
    #/group/bioinfo/apps/apps/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-63mer
  }
  elsif ($assembler eq "spades") {
    if (!defined($maxmem)) {$maxmem = 250;} #This is spades's default.
    $spadesworkdir = $spadesworkstem.$cycle;
    unless(mkdir $spadesworkdir) {die "Please clean out and delete $spadesworkdir before proceeding.\n";}
    $forwardextractedsfastq = $spadesoutstem."forwardextracteds.fastq";
    $reverseextractedsfastq = $spadesoutstem."reverseextracteds.fastq";
    for $key (keys(%readnames)) {
      print LOG "readnames key = $key\n";
      if ($readnamedelimiter eq '\.') {
        $fname = $key.".1";
        $rname = $key.".2";
      }
      else {
        $fname = $key.$readnamedelimiter."1";
        $rname = $key.$readnamedelimiter."2";
      }
      $forwardnames{$fname} = 1;
      $reversenames{$rname} = 2;
    }
    open (FSQ, "< $forwardreadfile");
    open (FFQ, "> $forwardextractedsfastq");
    $printflag = 0;
    while ($line = <FSQ>) {
      if ($line =~ m/^@/) {
        #@ERR3288215.1536.1 NB501206:425:HN5T2BGX7:1:11101:18885:1114 length=151
        @mars = split(/\s+/, $line);
        $mars[0] =~ s/@//;
        if (exists($forwardnames{$mars[0]})) {$printflag = 1;}
        else {$printflag = 0;}
      }
      if ($printflag == 1) {print FFQ $line;}
    }
    close (FFQ); close (FSQ);
    open (FSQ, "< $reversereadfile");
    open (FRQ, "> $reverseextractedsfastq");
    $printflag = 0;
    while ($line = <FSQ>) {
      if ($line =~ m/^@/) {
        #@ERR3288215.1536.1 NB501206:425:HN5T2BGX7:1:11101:18885:1114 length=151
        @mars = split(/\s+/, $line);
        $mars[0] =~ s/@//;
        if (exists($reversenames{$mars[0]})) {$printflag = 1;}
        else {$printflag = 0;}
      }
      if ($printflag == 1) {print FRQ $line;}
    }
    close (FRQ); close (FSQ);
    if (-s $forwardextractedsfastq < 10) {die "We failed to load the matching forward reads into a fastq file.\n";}
    if (-s $reverseextractedsfastq < 10) {die "We failed to load the matching reverse reads into a fastq file.\n";}
    #spades.py --careful -t 16 -o /scratch/halstead/c/ccrane/spadesonSRR3185380 -1 /scratch/halstead/c/ccrane/SRR3185380_1trimmedP.fastq -2 /scratch/halstead/c/ccrane/SRR3185380_2trimmedP.fastq --s1 /scratch/halstead/c/ccrane/SRR3185380_1trimmedU.fastq --s2 /scratch/halstead/c/ccrane/SRR3185380_2trimmedU.fastq
    print LOG "$spadesexe -t $nthreads --phred-offset $phredoffset -m $maxmem -o $spadesworkdir -1 $forwardextractedsfastq -2 $reverseextractedsfastq\n";
    `$spadesexe -t $nthreads --phred-offset $phredoffset -m $maxmem -o $spadesworkdir -1 $forwardextractedsfastq -2 $reverseextractedsfastq`;
    $tigfile = $spadesworkdir."/contigs.fasta";
    if (!-e $tigfile) {die "Post-assembly contig file $tigfile does not exist.\n";}
    $tggfile = $spadesoutstem.$cycle.".subset.contigs";
    print LOG "From cycle $cycle, tigfile = $tigfile and tggfile = $tggfile .\n";
  }
  elsif ($assembler eq "masurca") {
    $masurcaworkdir = $masurcaworkstem.$cycle;
    unless(mkdir $masurcaworkdir) {die "Please clean out and delete $masurcaworkdir before proceeding.\n";}
    $forwardextractedsfastq = $masurcaoutstem."forwardextracteds.fastq";
    $reverseextractedsfastq = $masurcaoutstem."reverseextracteds.fastq";
    @seqlengths = ();
    open (INP, "< $forwarddboutfile");
    open (FSQ, "> $forwardextractedsfastq") || die "We could not open $forwardextractedsfastq.\n";
    $lastdefline = "NULL";
    while ($line = <INP>) {
      if ($line =~ m/>/) {
        if ($lastdefline ne "NULL") {
          $lastdefline =~ s/>/@/;
          print FSQ $lastdefline, $sequence, "\n+\n";
          for ($i = 0; $i < length($sequence); $i++) {print FSQ "i";}
          print FSQ "\n";
          push(@seqlengths, length($sequence));
        }
        $lastdefline = $line;
        $sequence = "";
      }
      else {
        chomp $line;
        $sequence .= $line;
      }
    }
    if ($lastdefline ne "NULL") {
      $lastdefline =~ s/>/@/;
      print FSQ $lastdefline, $sequence, "\n+\n";
      for ($i = 0; $i < length($sequence); $i++) {print FSQ "i";}
      print FSQ "\n";
      push(@seqlengths, length($sequence));
    }
    close(INP); close(FSQ);
    open (INP, "< $reversedboutfile");
    open (FSQ, "> $reverseextractedsfastq") || die "We could not open $reverseextractedsfastq.\n";
    $lastdefline = "NULL";
    while ($line = <INP>) {
      if ($line =~ m/>/) {
        if ($lastdefline ne "NULL") {
          $lastdefline =~ s/>/@/;
          print FSQ $lastdefline, $sequence, "\n+\n";
          for ($i = 0; $i < length($sequence); $i++) {print FSQ "i";}
          print FSQ "\n";
          push(@seqlengths, length($sequence));
        }
        $lastdefline = $line;
        $sequence = "";
      }
      else {
        chomp $line;
        $sequence .= $line;
      }
    }
    if ($lastdefline ne "NULL") {
      $lastdefline =~ s/>/@/;
      print FSQ $lastdefline, $sequence, "\n+\n";
      for ($i = 0; $i < length($sequence); $i++) {print FSQ "i";}
      print FSQ "\n";
      push(@seqlengths, length($sequence));
    }
    $meanlength = 0;
    for ($i = 0; $i < scalar(@seqlengths); $i++) {$meanlength += $seqlengths[$i];}
    $meanlength /= scalar(@seqlengths);
    $ss = 0;
    for ($i = 0; $i < scalar(@seqlengths); $i++) {$ss += ($seqlengths[$i] - $meanlength)**2;}
    $ss /= scalar(@seqlengths);
    $ss = sqrt($ss);
    close(INP); close(FSQ);
    open (IMM, "> $immcfg");
    print IMM "DATA\n";
    print IMM "PE= aa $meanlength $ss $forwardextractedsfastq $reverseextractedsfastq\n";
    #Generate a file of long reads.
    print IMM "NANOPORE=$lrdboutfile\n";
    print IMM "END\n\nPARAMETERS\n";
    print IMM "NUM_THREADS=$masurcanumthreads\n";
    print IMM "END\n";
    close(IMM);
    die "bailing out\n";
  }
  else {die "Assembler $assembler is not supported in this program.\n";}
  $maxdeflinelength = 0;
  open (TIG, "< $tigfile");
  while ($vine = <TIG>) {
    if ($vine =~ m/>/) {
      if (length($vine) > $maxdeflinelength) {$maxdeflinelength = length($vine);}
    }
  }
  if ($maxdeflinelength > 50) { #This limit is imposed by makeblastdb. Edit $tigfile to shorten the deflines.
    seek(TIG, 0, 0);
    $ranno = int(rand(10000));
    @bars = split(/\//, $tigfile);
    $temptigfile = $bars[0];
    for ($j = 1; $j < scalar(@bars) - 1; $j++) {$temptigfile .= "/".$bars[$j];}
    $temptigfile .= "t".$ranno.$bars[-1];
    open (TTG, "> $temptigfile");
    if ($assembler eq "canu") { #Distinguishing substrings are at the front of the line.
      while ($vine = <TIG>) {
        if ($vine =~ m/>/) {
          @gars = split(/\s+/, $vine);
          print TTG "$gars[0] $gars[1] $gars[2]\n";
        }
        else {print TTG $vine;}
      }
    }
    else {
      while ($vine = <TIG>) { #Distinguishing substrings are at the back of the line.
        if ($vine =~ m/>/) {
          if (length($vine) > 30) {$vinet = ">".substr($vine, -30); print TTG $vinet;}
          else {print TTG $vine;}
        }
        else {print TTG $vine;}
      }
    }
    close (TIG); close(TTG);
    `mv $temptigfile $tigfile`;
  }
  print LOG "$blastdir/makeblastdb -in $tigfile -out $tempdbname -title \"Temporary database for cycle $cycle\" -dbtype nucl -parse_seqids\n";
  `$blastdir/makeblastdb -in $tigfile -out $tempdbname -title "Temporary database for cycle $cycle" -dbtype nucl -parse_seqids`;
  $examplefile = $tempdbname.".nhr";
  if (!-e $examplefile) {die "makeblastdb failed with $tigfile at cycle $cycle.\n";}
  $curroutfile = $tempdboutstem."_".$cycle."_".$evalue.".txt";
  if ($querytype =~ m/protein/) {
    print LOG "$blastdir/tblastn -db $tempdbname -query $foundingfile -out $curroutfile -evalue $evalue -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
    `$blastdir/tblastn -db $tempdbname -query $foundingfile -out $curroutfile -evalue $evalue -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
  }
  else {
    print LOG "$blastdir/blastn -db $tempdbname -query $foundingfile -out $curroutfile -evalue $evalue -dust no -outfmt \"6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send\"\n";
    `$blastdir/blastn -db $tempdbname -query $foundingfile -out $curroutfile -evalue $evalue -dust no -outfmt "6 qseqid sseqid evalue qacc bitscore score length nident positive mismatch pident ppos gaps gapopen sacc frames qframe qstart qend sframe sstart send"`;
  }
  if (-e $curroutfile) {
    @curroutfilestats = stat($curroutfile);
    if ($curroutfilestats[7] < 2) {die "Blastn of founding sequence against $tempdbname returned no hits at cycle $cycle.\n";}
  }
  else {die "$curroutfile was not produced by blastn of founding sequence against $tempdbname at cycle $cycle.\n";}
  open (TMP, "< $curroutfile"); #blast output
  @tempaccs = ();
  %tempfounds = ();
  while ($line = <TMP>) { #Carry forward only contigs that match the founding sequence sufficiently.
    chomp $line;
    @vars = split(/\t/, $line);
    #The carry-forward e-value can be more stringent than the initial e-value, since the contig will have fewer sequencing errors.
    if (!exists($tempfounds{$vars[1]})) {
      if ($vars[2] <= $carryforwardevalue) {
        push(@tempaccs, $vars[1]);
        $tempfounds{$vars[1]} = 1;
      }
    }
  }
  for ($i = 0; $i < scalar(@tempaccs); $i++) {print LOG "DIAG1: tempaccs[$i] = $tempaccs[$i] in cycle $cycle .\n";}
  if (scalar(@tempaccs) == 0) {
    print LOG "No contigs matched the founding sequence in cycle $cycle.\n";
    die "No contigs matched the founding sequence in cycle $cycle.\n";
  }
  %carriedsequences = ();
  open (TIG, "< $tigfile"); #Is this the right file of contigs from the latest round?
  while ($line = <TIG>) { #The first line of this fasta file is assumed to start with '>'.
    chomp $line;
    if ($line =~ m/>/) {
      $line =~ s/>//;
      ($curraccession, $filler) = split(/\s+/, $line, 2);
      $carriedsequences{$curraccession} = "";
    }
    else {$carriedsequences{$curraccession} .= $line;}
  }
  close (TIG);
  @tigfilestats = stat($tigfile);
  if ($tigfilestats[7] < 2) {die "Assembler $assembler produced no contigs in cycle $cycle.\n";}
  $lastnuccount = $tigfilestats[7];
  for $key (keys(%carriedsequences)) {print LOG "DIAG2: carriedsequences has key $key.\n";}
  @tempacclengths = ();
  open (TGG, "> $tggfile");
  for ($i = 0; $i < scalar(@tempaccs); $i++) {
    print TGG ">$tempaccs[$i]\n", wrap('', '', $carriedsequences{$tempaccs[$i]}), "\n";
    push(@tempacclengths, length($carriedsequences{$tempaccs[$i]}));
  }
  close(TGG);
  if (exists($manarray[0])) { 
    if ($cycle >= scalar(@manarray)) {
      print LOG "All entries in manarray have been used, exiting at the end of cycle $cycle.\n";
      die "All entries in manarray have been used, exiting at the end of cycle $cycle.\n";
    }
  }
  elsif (defined($minprogress)) {
    @sortedtempacclengths = sort {$b <=> $a} @tempacclengths;
    if (exists($lastsortedtempacclengths[0])) {
      if ($sortedtempacclengths[0] < $lastsortedtempacclengths[0] + $minprogress) {
        print LOG "Longest contig did not grow enough, exiting at the end of cycle $cycle.\n";
        die "The longest contig did not grow enough, exiting after cycle $cycle.\n";
      }
    }
    @lastsortedtempacclengths = @sortedtempacclengths;
  }
  if ($useracon == 0) {$seqfile = $tggfile;}
  else { #$readsfile and $poblastdb are generated outside of this program.
    if (!-e $readsfile) {
      $readsfile =~ s/fastq$/fq/;
      if (!-e $readsfile) {die "$readsfile was not found.  Please check its path and spelling.\n";}
    }
    $currhr = (localtime)[2];
    $currmin = (localtime)[1];
    $currsec = (localtime)[0];
    print LOG "Racon began $nrounds of polishing at $currhr", ":$currmin", ":$currsec.\n";
    for ($i = 0; $i < $nrounds; $i++) {
      $reffile = $postem."_cyc".$cycle."_refseq_rd".$i.".fasta";
      if ($i == 0) {`cp $tggfile $reffile`;}
      else {`cp $polfile $reffile`;}
      $bowtie2idxstem = $postem."_cyc".$cycle."_forbowtie2_rd".$i;
      $samfile = $postem."_cyc".$cycle."_bowtie2out.sam";
      $polfile = $postem."_cyc".$cycle."_polished_rd".$i.".fasta";
      $blastout = $postem."_cyc".$cycle."_polished_vs_ref_rd".$i."_".$evalue.".txt";
      print LOG "$bowtie2dir/bowtie2-build $reffile $bowtie2idxstem\n";
      `$bowtie2dir/bowtie2-build $reffile $bowtie2idxstem`;
      print LOG "$bowtie2dir/bowtie2 --local --no-unal -x $bowtie2idxstem -U $readsfile -S $samfile\n";
      `$bowtie2dir/bowtie2 -p $nthreads --local --no-unal -x $bowtie2idxstem -U $readsfile -S $samfile`;
      if (!-e $samfile) {die "Previous bowtie2 step failed.\n";}
      elsif (-z $samfile) {die "Previous bowtie2 step failed.\n";}
      print LOG "$raconexe $readsfile $samfile $reffile > $polfile\n"; #Remember to module load Racon; there is a dependency.
      `$raconexe $readsfile $samfile $reffile > $polfile`;
      if (!-e $polfile) {die "Previous Racon step failed.\n";}
      elsif (-z $polfile) {die "Previous Racon step failed.\n";}
      #print LOG "$blastdir/blastn -query $polfile -db $blastdb -out $blastout -evalue $poevalue\n";
      #`$blastdir/blastn -query $polfile -db $poblastdb -out $blastout -evalue $poevalue`;
    }
    $seqfile = $polfile;
  }
  if (defined($maxcycle)) {
    if ($cycle >= $maxcycle) {
      print LOG "Maximum cycle is on. The program is exiting at the end of cycle $cycle.\n";
      die "The program is exiting at the end of maximum cycle $cycle\n";
    }
  }
  $lastnacc = $nacc;
  $cycle++;
  $currhr = (localtime)[2];
  $currmin = (localtime)[1];
  $currsec = (localtime)[0];
  print LOG "Cycle $cycle began at $currhr", ":$currmin", ":$currsec.\n";
}#end of the cycle loop
open (TGG, "< $tggfile");
$sequence = "";
while ($line = <TGG>) {
  chomp $line;
  if ($line =~ m/>/ || eof(TGG)) {
    $seqlength = length($sequence);
    if ($seqlength > 2) {print LOG "$header\tlength = $seqlength\n";}
    $header = $line;
    $sequence = "";
  }
  else {$sequence .= $line;}
}
$currhr = (localtime)[2];
$currmin = (localtime)[1];
$currsec = (localtime)[0];
print LOG "Program ended at $currhr", ":$currmin", ":$currsec\n";
