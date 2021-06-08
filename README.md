	SLAG is an iterative targeted assembly pipeline that is strategically
similar to aTRAM2* but offers a different set of assemblers and can produce
assemblies at low read coverage. At the first cycle, one or more protein or
nucleotide seeding sequences are BLAST aligned to a database of reads. Matching
reads are retrieved and assembled, and the resulting contigs are screened by
BLAST against the seeding sequences. Matching contigs go on to the next cycle
of blast-retrieve-assemble-screen. SLAG currently can use CAP3, phrap, SPAdes,
canu, and Unicycler, and it can fragment long reads for assembly with CAP3 or
phrap when coverage is insufficient for canu or Unicycler.

	The current edition of SLAG is localassembly1115.pl, which takes a
single argument, a configuration file that specifies up to 40 parameter values.
Non-comment lines in this file have the form "$variable = $value" for example,
"$blastexe = "/bin/blastn";" without the outer quotes. Comment lines begin with
"#". The output of SLAG is a series of "*subset.contigs" fasta files, one per
completed cycle. These files can be parsed with adhoc03082021a.pl or
quickacciostats.pl to identify the longest contig, the number of contigs, and
the mean length of contigs, for example, by "./adhoc03082021a.pl
filelist03082021k.txt Zeasucrosesynthasesresults.txt" or "./quickacciostats.pl
attemptedaccs1023A_1119.txt accs1023AmatchinIWGSPacciostats1202.txt
accs1023AinIWGSPacciobrief1202.txt /scratch/lustreC/c/ccrane
1023A_extractedIWGSPseqs".

	Several configuration files are provided here as examples for
particular read lengths and assemblers:

SPAdes assembly of short reads: Zeahexokinasespureslagbm.cfg

CAP3 assembly of short reads: Stanleyhexokinasespureslagbm.cfg

CAP3 assembly of singly fragmented long reads:
ZeahexokinasesLR1cap3increment.cfg

Phrap assembly of singly fragmented long reads:
ZeahexokinasesLR1phrapincrement.cfg

CAP3 assembly of doubly fragmented long reads:
TraesCS1A01G421800.1_LR2_cap3_increment_doublefrag.cfg

Assembly of long reads with canu:
TraesCS1B01G473100.1_CANU_canu_increment_intact.cfg

Unicycler assembly of long reads: Zeahexokinases_unic.cfg

Manual read retrieval with phrap assembly:
configTraesCS1A01G397600formanual.cfg

Phrap assembly with read retrieval based on bitscore:
configTraesCS1A01G397600forbitscore.cfg

	Settable parameters and allowed values or types in configuration files
include:

$assembler   "unicycler", "canu", "cap3", "phrap", "spades"

$blastdb    name of the reads database for standalone BLAST 

$blastdir    path to the BLAST executable, e.g.,
"/group/bioinfo/apps/apps/blast-2.7.1+/bin"

$canuexe    location and name of the canu executable, e.g.,
"/group/bioinfo/apps/apps/canu-1.8/Linux-amd64/bin/canu"

$canuprefix   base name for canu output files, e.g., "Zeahistonedeacetylases"

$canusettings  settings passed to canu, e.g., " "-p $canuprefix -d
canuworkdirectory genomeSize=estgenomesize usegrid=false
correctedErrorRate=0.105 -nanopore-raw"

$canuworkstem	path to and base name of canu’s working directories, of which
there is one per cycle. Example:
"/scratch/halstead/c/ccrane/slagtests/Zeahistonedeacetylasescanuwork"

$cap3exe    location and name of the CAP3 executable, e.g.,
"/group/bioinfo/apps/apps/CAP3-12.21.07/cap3"

$cap3options  settings for CAP3, e.g., "-b 20 -m 2 -n -4 -g 5 -s 600 -p 83 -o
40 -y 150 -z 3 -h 25 -j 70"

$cap3outstem	path and base name of CAP3 output files, e.g., 
"/scratch/halstead/c/ccrane/slagtests/Zeahistonedeacetylases_LR1cap3_assembly"

$carryforwardevalue   maximum accepted e-value for contig hits against founding
sequences, e.g., 1e-20

$chunksize	fragment size of singly-fragmented long reads for phrap or CAP3
assembly, e.g., 600

$chunksizea	first fragment size of doubly fragmented long reads for phrap
or CAP3 assembly, e.g., 490

$chunksizeb	second fragment size of doubly fragmented long reads for phrap
or CAP3 assembly, e.g., 610

$evalue BLAST e-value for first-cycle alignment of reads to seeding sequences,
e.g.,1e-10

$extincrement	number of reads to increase in each subsequent cycle, e.g., 30.
These reads come from the top of tabular BLAST output and therefore are the
closest matching reads.

$extractionoption	how to use reads that pass the BLAST e-value. Allowed
values are "increment", "population", "manual", "all", or "bitscore". We
recommend "increment" with a value near the read depth to twice the read depth,
e.g., 20.

$forwardblastdb path and name of a forward-reads BLAST database, e.g.,
"/scratch/brown/ycrane/dbforslag/ERR328821xcleaned_1_db_R1"

$forwardreadfile	file of forward reads from paired-end reads, e.g.,
"/scratch/brown/ycrane/ERR328821xcleaned_1.fastq"

$fragseqfilestem	path to and base name of file of fragmented long reads,
e.g.,	 
"/scratch/halstead/c/ccrane/slagtests/ZeahistonedeacetylasesphrapLR1extractedfr
agments"

$lastnuccount	value to substitute for the string "estgenomesize" in
canusettings, e.g, 10000, size of the expected longest contig

$longread	Set this to 1 to singly fragment long reads for CAP3 or phrap,
2 for doubly fragmenting long reads for CAP3 or phrap, or 0 otherwise.

@manarray	Array of numbers of reads to retrieve for manual read
retrieval, e.g. for six cycles, (17, 24, 32, 40, 48, 57).

$maxcycle	number of cycles to run after the initial cycle with seeding
sequences, e.g., 20

$maxnacc	maximum number of retrieved reads to allow for phrap assembly,
e.g., 60000

$minbitscore	minimum BLAST bitscore to accept reads, e.g., 5000

$minprogress	minimum allowed length increase for the longest contig over the
previous cycle, e.g., 100. This variable usually is not set with the current
usage, which allows growth and shrinkage over consecutive cycles.

$nthreads	number of threads for BLAST searches, e.g., 4

$pairedend	Set to 1 for paired-end reads, 0 for singleton reads. 2 for
paired-end reads with Racon polishing.

$phrapexe	path and file name of the phrap executable, e.g.,
"/group/bioinfo/apps/apps/genome/bin/phrap"

$phrapoutstem	path and base name for phrap output files, e.g., 
"/scratch/halstead/c/ccrane/slagtests/Zeahistonedeacetylases_LR1phrap_assembly"

$phrapsettings	settings for phrap, e.g., "-minmatch 18 -maxgap 20
-repeat_stringency 0.85 -retain_duplicates" 

$phredoffset	value of phred offset for SPAdes assembly, e.g., 33

$querytype	type of seeding sequence, "protein" or "nucleotide"

$restartflag	1 to restart a run that ran out of time before completing, 0
otherwise

$reverseblastdb path and name of reverse-reads database if reads are paired,
e.g., "/scratch/brown/ycrane/dbforslag/ERR328821xcleaned_1_db_R2"

$reversereadfile	path and name of reverse reads fasttq file, e.g.,
"/scratch/brown/ycrane/ERR328821xcleaned_2.fastq"

$runalign	 value of BLAST mas_target_seqs, e.g., 10000

$secevalue		e-value for BLAST alignments of contigs to reads, e.g.,
1e-20

$seqfile	fasta file of query sequences, e.g.,
"/scratch/halstead/c/ccrane/slagtests/Zeahistonedeacetylases.fasta"

$spadesexe	SPADes executable with path, e.g.,
"/group/bioinfo/apps/apps/spades-3.14.1/bin/spades.py"

$spadesoutstem		prefix for names of SPAdes output files, e.g.,
"Zeahexokinasesbm_assembly"

$spadesworkstem prefix of SPAdes working directories, e.g.,
"/scratch/brown/ycrane/Zeahexokinasesbmspadeswork"

$stem	prefix of intermediary file names, e.g., "Zeahistonedeacetylases"

$tempdbname	name of contig database to be aligned with seeding sequences,
e.g., "Zeahistonedeacetylasesunicretrieveddb"

$tempdboutstem	prefix of BLAST output files for alignments to the contig
database, e.g., "Zeahistonedeacetylases_vs_Zeahistonedeacetylasesunicretrieved"

$unicyclerexe	name with path of the Unicycler executable, e.g., "unicycler"

$unicycleroutstem	prefix of Unicycler output files, e.g.,
"Zeahistonedeacetylases_unicycler_assembout"

$unicyclerworkstem	prefix of Unicycler working directories, e.g.,
"/scratch/halstead/c/ccrane/Zeahistonedeacetylasesunicyclerwork"


*Allen, J. M., LaFrance, R., Folk, R. A., Johnson, K. P., & Guralnick, R. P.
(2018). aTRAM 2.0: an improved, flexible locus assembler for NGS data.
Evolutionary Bioinformatics, 14, 1-4.

