# SLAG
     SLAG (Seeded Local Assembly of Genes) is a local genome assembler.  It is stategically similar to aTRAM and aTRAM2, but it was independently derived and goes 
its own way with strategies to leverage low-copy read libraries.  SLAG iteratively collects reads that match one or more target sequences, assembles them, and feeds
target-matching contigs into the next iteration.  A manuscript has been submitted to Molecular Ecology Resources describing SLAG:
     SLAG: A Program for Seeded Local Assembly of Genes in Complex Genomes.  2020.  Charles F. Crane, Jill A. Nemacheck, Subhashree Subramanyam, Christie E. Williams, and Stephen B. Goodwin.
     SLAG's executable is named localassembly1115.pl, and it takes a single argument, the name of a configuration file.  Here is a typical invocation in a Unix-like
environment:
     ./localassembly1115.pl Zeahexokinasescanuincrement.cfg
The content of the configuration file depends strongly upon the choice of assembler and the length of available reads.  In some cases, the assembler depends upon a 
machine-specific, institutional installation-specific environment, which must be set up outside of SLAG.  Please peruse the provided examples of configuration files for
settings that have been used with your choice of assembler.  SLAG currently supports SPAdes for short reads (e.g., 2x150), cap3 and phrap for intermediate (Sanger) reads, and canu and Unicycler for long (Nanopore or PacBio) reads.
     Here the file Zeahexokinasescanuincrement.cfg is annotated to illustrate the roles of the variables.  Lines that begin with "##" are comments that are not part
of the configuration file itself.
#file section
$seqfile = "/scratch/halstead/c/ccrane/slagtests/Zeahexokinases.fasta";
##$seqfile is the fasta-formatted file of target sequences.
$canuoutstem = "/scratch/halstead/c/ccrane/slagtests/Zeahexokinases_assembout";
##Any variable named ..outstem designates an output directory and prefix to the names of SLAG's output files.
#general section
$restartflag = 0;
##Start this run from the target sequences, do not pick up after an earlier run.
$maxcycle = 20;
##Iterate 20 cycles after the initial collection of reads that match the target sequence(s), for a total of 21 cycles.
$extractionoption = "increment";
##Add a fixed number of reads to the assembly per cycle.  Other allowed values of $extractionoption are "manual", "all", "bitscore", and "population".
$extincrement = 30;
##The number of reads to add per cycle, which depends on sequencing depth.
$longread = 0;
##Use intact long reads; do not fragment them. $longread = 1 specifies to break the reads up into fragments of length specified by $chunksize.  $longread = 2
##specifies to fragment the reads first by one chunk size and then by another chunksize, then to combine the results for the assembler.
#blast section
$querytype = "protein";
##The alternative is "nucleotide".
$blastdir = "/group/bioinfo/apps/apps/blast-2.7.1+/bin";
##Where the blast executable resides on your installation.
$blastdb = "\"/scratch/halstead/c/ccrane/ERR3288290db /scratch/halstead/c/ccrane/ERR3288291db /scratch/halstead/c/ccrane/ERR3288292db /scratch/halstead/c/ccrane/ERR3288293db /scratch/halstead/c/ccrane/ERR3288294db /scratch/halstead/c/ccrane/ERR3288295db\"";
##The blast database(s) of reads.  Read sequences are retrieved from the database by blastdbcmd.  Note in this example how the inner quotes have been escaped.
$nthreads = 4;
##Number of threads for blastn/blastx to use.
$evalue = 1e-10;
##Blast e-value for initial read matches to target sequences, applied only in the first cycle.
$secevalue = 1e-20;
##Blast e-value for subsequent cycles.
$runalign = 10000;
##Number of blast alignments to get.  The number seems excessive, and possibly it is excessive, but we want to make sure that all the closest reads are examined
#and collected.
$carryforwardevalue = 1e-20;
##Blast e-value for contigs versus target sequences.  Only contigs that match the target at least this well proceed to the next cycle.
$stem = "/scratch/halstead/c/ccrane/slagtests/Zeahexokinases";
##Directory and prefix for more output files, including queries of contigs against the reads in subsequent cycles.
$tempdbname = "/scratch/halstead/c/ccrane/slagtests/Zeahexokinasescanuretrieveddb";
##Base name of the blastable database of contigs, to be queried by the target sequence(s) in $seqfile or prior-run contigs in $foundingfile. 
$tempdboutstem = "/scratch/halstead/c/ccrane/slagtests/Zeahexokinases_vs_Zeahexokinasesretrieved";
##Prefix for blast output versus $tempdbname.
#assembler section
$assembler = "canu";
##The name of the chosen assembler.  Other options include "phrap", "cap3", "spades", and "unicycler".
$canuexe = "/group/bioinfo/apps/apps/canu-1.8/Linux-amd64/bin/canu";
##The name of the assembler executable, including any necessary path.
$canuprefix = "Zeahexokinases";
##Different assemblers have different variables.
$canusettings = "-p $canuprefix -d canuworkdirectory genomeSize=estgenomesize usegrid=false correctedErrorRate=0.105 -nanopore-raw";
##Invoke canu with these settings.  In general, assembler settings are inserted in the assembler invocation as a single string.
$canuworkstem = "/scratch/halstead/c/ccrane/slagtests/Zeahexokinasescanuwork";
##Directory basename for canu's temporary files.  A cycle number is appended to $canuworkstem to make the directory for each cycle.
$lastnuccount = 10000;
##Estimated size of assembly; updated at each cycle by total contig length in the previous cycle.
##END OF CONFIGURATION FILE.
