<h1 align="center">SLAG</h1>
<p align="center">
</p>

SLAG is a local assembly program that is strategically similar to aTRAM2 (Allen et al., 2018, Evolutionary Bioinformatics 14: 1-4). However, SLAG uses a different set of assemblers ([`canu`](https://github.com/marbl/canu),
[`CAP3`](https://faculty.sites.iastate.edu/xqhuang/cap3-and-pcap-sequence-and-genome-assembly-programs), [`phrap`](http://www.phrap.org/phredphrapconsed.html), [`SPAdes`](https://github.com/ablab/spades),  [`Unicycler`](https://github.com/rrwick/Unicycler)) and can handle low read depths where aTRAM2 fails to yield assemblies. The usual use cases will be mining promoters using EST or heterologous protein sequences, and mining alleles or paralogs or homeologs using nucleotide or protein sequences. 

A manuscript describing SLAG and its local assemblies has been submitted to Molecular Ecology Resources: Charles F. Crane, Christie E. Williams, Jill Nemacheck, Subhashree Subramanyam, and Stephen B. Goodwin. [expectedly 2022] SLAG: A Program for Seeded Local Assembly of Genes in Complex Genomes. Molecular Ecology Resources xx: xx-xx.

### Features
SLAG:
* can assemble short, medium, or long reads
* has obtained assemblies from read coverage **as low as 5x**, which is an advantage when working with large, polyploid genomes
* is about as fast as aTRAM2, although either SLAG or aTRAM2 can be up to 15 times faster than the other for particular inputs
* is more robust than aTRAM2, especially at lower coverage and when SLAG uses incremental addition of reads to the set to be assembled

### Installation
Please see the INSTALL.md file in this distribution. Basically, one downloads the distribution file SLAGball.tgz from the releases section of this page, puts it where SLAG should go, uncompresses it, and untars it. One then has the option to add this path to the user's $PATH.

### Usage
This distribution supplies two SLAG executables: a "free-living" `SLAG.pl` that requires the user to know the file locations of various dependencies, and a Singularity-containerized `SLAG.sif` that carries all the dependencies within it.

The containerized form is more portable and likely to outlast obsolescence of SLAG's current dependencies. However, because of a licensing restriction, the containerized executable ⚠️ **does not** run the phrap assembler. 

Before you begin, SLAG requires some preparation. You are presumed to have a
* fasta or fastq file(s) of reads, and 
* a fasta file of one or more founding nucleotide or protein sequences. 
 
SLAG.pl requires access to `blast` and at least one of these programs: [`CAP3`](https://faculty.sites.iastate.edu/xqhuang/cap3-and-pcap-sequence-and-genome-assembly-programs), [`phrap`](http://www.phrap.org/phredphrapconsed.html), [`SPAdes`](https://github.com/ablab/spades), [`canu`](https://github.com/marbl/canu), [`Unicycler`](https://github.com/rrwick/Unicycler), and [`Racon`](https://github.com/lbcb-sci/racon). You will either need to know the fully qualified location of these programs in your filesystem, or have access to a working "module load" command that can modify your `$PATH` to include the program's location. 

With either `SLAG.pl` or `SLAG.sif`, you will need to use `makeblastdb` to make a blastable database of your reads. 

If you start with paired- end `fastq` files, you will need to copy and convert them to `fasta` before preparing separate forward- and reverse-read databases. 

Invoking SLAG is simple: `./SLAG.pl configuration_file` or, if your system supports Singularity: `./SLAG.sif configuration_file` .

Setting up the configuration file is not so simple. There is expectedly a different configuration file for every run of SLAG; the configuration file is the only way for the user to direct SLAG for each specific local assembly. Because SLAG runs <u>often last for hours</u>, SLAG will usually be invoked within a shell script or at least a Slurm sbatch command. 

⚠️ Note for using `SLAG.sif`: This container can see `/home` and `/scratch` directories in the external (host) filesystem. This includes all subdirectories under the external `/scratch`. If your blastable databases and reads are in some directory $outsidedir that is not rooted in `/scratch` or `/home` on the host, you will have to use singularity run -B: `singularity run -B $outsidedir:/scratch SLAG.sif $configurationfile` .

#### Structure of the configuration file
The configuration file is a series of statements of the form 
```
$parameter1 = $value1;
$parameter2 = $value2;
```
without any quotes, i.e., a series of assignment statements in perl. 

SLAG has 59 settable parameters. Typically the configuration file will define 25 to 30 parameters, depending greatly on which local assembler is to be used and whether an assembly of long reads is being attempted from shallow read depth. 

Users are **strongly encouraged** to use one of the 24 template configuration files provided with this release in subdirectory `slagbase/templates`. These templates cover the common combinations of read length, read coverage, and assembler.  An ancillary Perl script, `setupSLAGconfigurations.pl`, can substitute user-supplied file names into a copy of the template to make a run-specific configuration file. 

Other changes to the configuration file require hand-editing and some familiarity with the assemblers. Since there are 59 settable parameters, many combinations of values have never been tested, and some are internally inconsistent.

#### Usage of `setupSLAGconfigurations.pl`: 
`./setupSLAGconfigurations.pl $templatefile $newfile AAAA=$assemblername DDDD=$databasename ...`

where each argument after the second is a pairing of coding letters with a name or value. For example, if your institution has "module load" and a canu module, an argument could be as simple as `AAAA=canu`. 

Otherwise, the argument could be `AAAA=/wherever/the/sysadmin/put/it/canu` .

Because `setupSLAGconfigurations.pl` will substitute whatever is to the left of the equals sign with whatever is to the right of it, be careful if you create a new template file that you do not duplicate an unrelated string with the code letters.

#### Output
Regardless of the assembler, SLAG will iterate through some number of cycles of read retrieval and assembly, outputting a distinct [prefix][cycle].subset.contigs file at each successfully completed cycle. This fasta file contains the locally assembled contigs that match the founding, target sequence. The user sets the value of prefix in the configuration file. For canu and Unicycler, a .subset.contigs file resides in the working directory produced by each cycle.

### Variables

#### Common settable variables

| Variable | Explanation |
| ----------| ----------- |
| $assembler | The assembler to be used in this run of SLAG. Each run uses only one assembler. The choices are "canu", "cap3", "phrap", "spades", and "unicycler"; there is no default choice. Several variables are specific to the choice of assembler. They are listed below as "Variables specific to $assembler assembly". |
| $carryforwardevalue | The maximum e-value allowed for contigs to pass on to the next cycle. This e-value comes from blastn or tblastn of the latest contigs against the founding sequence. Because the alignments here are of whole contigs to the founding sequence, the value of $carryforwardevalue should be less than or equal to the value of $secevalue. This variable affects how many contigs are carried forward to seed the next cycle. Reasonable values range from 1e-50 to 1e-20. |
| $extincrement | The maximum number of additional matching reads to include at each cycle after the first, if $extractionoption is set to "increment". All matching reads are assembled at the first cycle. A good starting value for $extincrement is the read depth.|
| $extractionoption | The method to select reads to be assembled. The choices are "all", "bitscore", "increment", "manual", and "population". In each cycle of SLAG, a blast alignment returns a list of read identifiers in ascending order of e-value, i.e., the closest alignments are at the top. <br/><br/>Option "all" causes SLAG to assemble all matching reads, so that the number of reads used is determined by the blast e-value.<br/><br/>Option "bitscore" uses only matching reads that exceed a bitscore value specified by $minbitscore. <br/><br/>Option "increment" causes a maximum of $extincrement additional reads from the top of the list to be assembled at each cycle. <br/><br/>With option "manual", the number of reads to be assembled is set by the user, and this number can vary among cycles. The values are set as an array @manarray, for example by @manarray = (30, 40, 50, 60, 70), where the numbers within the array are chosen to try to find a maximum contig length. It makes sense that the values are in ascending order in @manarray, but SLAG will run even if they are not. <br/><br/>Option "population" takes reads from the list of ascending e-values until all reads used in the previous cycle have been chosen. Any additional reads to be assembled in the current cycle must have an e-value less than or equal to the largest e-value taken in the previous cycle. <br/><br/>Option "all" (along with setting $minprogress to a small positive number) makes SLAG behave like aTRAM2, while option "increment" allows the longest contig length to fluctuate and generally has produced the longest contigs with the deepest penetration of flanking repetitive sequence at some point during the run.  However, if $extincrement is too small, contig growth is noticeably slowed. |
| $cycle | If SLAG is continuing a previous run, this is the cycle number for this run to begin with. It should be <= $maxcycle. If SLAG is not continuing a previous run, do not set $cycle. |
| $foundingfile | If SLAG is continuing a previous run, this is the *subset.contigs file from the last completed cycle of the previous run. If SLAG is not continuing a previous run, do not set $foundingfile. |
| $longread | Choice of a read-fragmentation method to enable long reads to be assembled with cap3 or phrap when read depth is too shallow to run canu or unicycler. Meaningful values are 1 for single fragmentation and 2 for double, staggered fragmentation. Leaving $longread unset or specifying any other value does not result in any fragmentation. |
| $maxcycle  | SLAG performs this number of cycles after the first cycle. For example, setting $maxcycle = 20 causes SLAG to go through a maximum of 21 cycles. Because of a hard-coded print statement in SLAG, the maximum allowed value is 39. |
| $minbitscore | The minimum blast bitscore for accepting a read for assembly. This variable applies only if $extractionoption = "bitscore". |
| $minprogress | Minimum length increase of the longest contig from the previous cycle. Stop the run if the length increase is less than this value. If $minprogress is not set, the run can continue regardless of contig expansion or contraction in consecutive cycles. |
| $pairedend | A flag indicating if reads are from paired ends. Allowed values are 0, 1, and 2. Zero indicates that the reads are not from paired ends. One and two indicate that reads are from paired ends; two also specifies that Racon polishing will be performed. For values 1 or 2, there will be twin blast jobs against twin databases, where forward and reverse reads are handled separately. |
| $querytype | The type of founding, target sequence. Allowed values are "nucleotide" and "protein". |
| $readnamedelimiter | A character or short string on which to split read names that distinguish forward and reverse paired-end reads. The default value is "." if this parameter is not specified. Example: forward read is read_0.1 and reverse read of the same pair is read_0.2. A frequent alternative value is "/". |
| $restartflag | A flag indicating if this run is to continue a previous run of SLAG, for example, if the previous run ran out of time. Allowed values are zero and 1. Zero indicates a fresh run, and 1 directs continuation of a previous run. If continuing a previous run, $foundingfile and $cycle must be set. |
| $runalign | An integer that specifies the value of -max_target_seqs for blast searches. |
| $seqfile | The fasta file that contains the founding, target sequence, which can be nucleotide or protein as indicated by the value of $querytype. |
| $stem | Prefix for a number of output and database file names. It should be distinct for each input locus and sample to be assembled. |

#### blast-specific variables

| Variable | Explanation |
| ----------| ----------- |
| $blastdb | File name, without the three-letter extension, of a blastable database of reads, from which reads will be selected for assembly based on e-value of alignment to a founding, target sequence. If there is more than one database, the list of space-delimited names should be enclosed in double quotes. |
| $blastdir | Directory where blastn, tblastn, and blastdbcmd are located. On a personal filesystem, "find / -name "blastn"" should reveal this directory. On a cluster where "module load blast" works, a "which blastn" will reveal it. |
| $evalue | The e-value to be used for blastn or tblastn of the founding sequence against the reads database at the first cycle. |
| $nthreads | The number of threads for blast to use anywhere in the cycle. |
| $runalign | An integer that specifies the value of -max_target_seqs for blast searches, where it is the maximum number of alignments to keep. Blast's default is 500, but SLAG often needs a value of 10000. |
| $secevalue | The e-value to be used for blastn of carried-forward contigs against the reads database from the second cycle onward. This will generally be less than the primary $evalue for the founding sequence against the reads database. A typical value is 1e-20. |
| $tempdbname | Name of a blastable database rebuilt from all assembled contigs near the end of each cycle. A blastn or tblastn aligns the founding sequence to this database. On the basis of the alignment e-values, some or all contigs will be carried forward to the next cycle. |
| $tempdboutstem | Base name for files output by blastn or tblastn of the founding sequence against the database $tempdbname. The output file names are concatenations of $tempdboutstem, the cycle number, $evalue, and ".txt". |

#### canu assembly-specific variables

| Variable | Explanation |
| ----------| ----------- |
| $canuexe | Location and name of the canu executable. To find this name on a personal filesystem, run for example "find / -name "canu"". On a computing cluster with loadable modules, run "module avail" and search for canu, then run "module load canu", and finally run "which canu". Use the full path that is returned. |
| $canuprefix | Prefix for output files from canu, including path. |
| $canusettings | A double-quoted string of settings passed verbatim to canu, e.g., "-p $canuprefix -d canuworkdirectory genomeSize=estgenomesize usegrid=false correctedErrorRate=0.105 -nanopore-raw". SLAG replaces the words canuworkdirectory and estgenomesize with the names of the upcoming canu working directory (which SLAG composes using $canuworkstem) and the number of nucleotides in the .subset.contigs file output by the previous cycle. |
| $canuworkstem | Canu uses a separate working directory for each SLAG cycle. Its name is a concatenation of $canuworkstem and $cycle. The canu working directories should not exist before starting SLAG. |

#### cap3 assembly-specific variables

| Variable | Explanation |
| ----------| ----------- |
| $cap3exe | Location and name of the CAP3 executable. To find this name on a personal filesystem, run "find / -name "cap3"". On a computing cluster with loadable modules, run "module avail" and search for "CAP3" or "cap3", then run "module load cap3[or CAP3]", and finally run "which cap3" or "which CAP3". It is safest to use the full path that these commands return. |
| $cap3options | A double-quoted string of options to be passed verbatim to cap3, e.g., "-b 20 -m 2 -n -4 -g 5 -s 600 -p 85 -o 40 -y 150 -z 3 -h 25 -j 70". Look in the cap3 documentation for further explanation. |
| $cap3outstem | Prefix with path for cap3 output files. |

#### SPAdes assembly-specific variables

| Variable | Explanation |
| ----------| ----------- |
| $forwardblastdb | Full path and name of a blastable database of forward reads. <br/><br/>SLAG blasts the sequences in $seqfile against this database and retrieves reads from it. The user must prepare this database beforehand by converting the forward-read fastq file to fasta and using makeblastdb. |
| $forwardreadfile | Full path and name of the forward fastq file of Illumina reads, generally containing "R1" in the name. SLAG awkwardly but necessarily goes back and forth between fastq and fasta files as it sets up the next SPAdes assembly, and the fastq file is needed because blastdbcmd retrieves reads in fasta format. |
| $phredoffset | Phred offset used by SPAdes to interpret input quality scores. <br/><br/>Allowed values are 33 and 64. All of our tests have used 33. |
| $reverseblastdb | Full path and name of a blastable database of reverse reads.<br/><br/>SLAG also blasts the sequences in $seqfile against this database and retrieves reads from it. This database must be prepared by the user beforehand by converting the reverse-read fastq file to fasta and using makeblastdb. |
| $reversereadfile | Full path and name of the reverse fastq file of Illumina reads, generally containing "R2" in the name. SLAG awkwardly but necessarily goes back and forth between fastq and fasta files as it sets up the next SPAdes assembly, and the fastq file is needed because blastdbcmd retrieves reads in fasta format. |
| $spadesexe | Location and name of the SPAdes assembler. To find this name on a personal filesystem, run "find / -name "spades.py"". On a computing cluster with loadable modules, run "module avail" and search for "spades", run "module load spades", and run "which spades". It is safest to use the full path that these commands return. |
| $spadesoutstem | Prefix of the name of the output file from SPAdes. The contigs that still match the founding, target sequence are in a file named by concatenating $spadesoutstem, $cycle, and ".subset.contigs". |
| $spadesworkstem | SPAdes creates a new directory at each cycle, and $spadesworkstem is the prefix of its name. SLAG will complain and halt if these directories are not empty when the run begins. |

#### Unicycler assembly-specific variables

| Variable | Explanation |
| ----------| ----------- |
| $unicyclerexe | Location and name of the Unicycler assembler. To find this name on a personal filesystem, run "find / -name "unicycler" ". On a computing cluster with loadable modules, run "module avail" and search for "Unicycler" or "unicycler", run "module load Unicycler", and run "which unicycler". |
| $unicycleroutstem | This is the prefix of the output file name, which is a concatenation of $unicycleroutstem and ".subset.contigs". This file contains contigs that match the founding, target sequence. |
| $unicyclerworkstem | This is the prefix of a subdirectory name. SLAG makes a separate subdirectory for each cycle, and the subdirectory names contain the cycle number. For each successfully completed cycle, this subdirectory will contain a file simply ending in ".contigs" and another file ending in ".subset.contigs". The latter is a subset of the former and contains only contigs that match the target sequence. It is recommended to set a different value of $unicyclerworkstem for each run of SLAG beforehand. |

#### phrap assembly-specific variables

| Variable | Explanation |
| ----------| ----------- |
| $maxnacc | The maximum number of reads to submit to phrap assembly. Without the -manyreads option, phrap can handle a maximum of about 65000 reads, and values of $maxnacc betwee 30000 and 60000 have proven less likely to cause phrap to abort. |
| $phrapexe | Location and name of the phrap executable. To find this name on a personal filesystem, run "find / -name "phrap" ". On a computing cluster with loadable modules, run "module avail" and search for "phrap", then run "module load phrap", and finally run "which phrap". It is safest to use the full path that these commands return. |
| $phrapoutstem | Prefix for phrap output file names. This should include a path so that the results are not written to the same directory where phrap resides. |
| $phrapsettings | A double-quoted string of argument names and values that is inserted into the command that invokes phrap, such as "-minmatch 18 -maxgap 20 - repeat_stringency 0.85 -retain_duplicates". The value of repeat_stringency largely controls the stringency of phrap assembly, i.e., how much differing reads are collapsed into consensus contigs. The default value in phrap is 0.95, and 0.83 is as low as we have dared to go. A value above 0.97 risks splitting the assembly into a contig for each read, sequencing errors and all. |

#### Variables specific to fragmented reads for cap3 or phrap assembly

| Variable | Explanation |
| ----------| ----------- |
| $chunksize | Long reads will be fragmented to this size if $longread is set to 1. |
| $chunksizea | Long reads will be fragmented to this size in a first pass if $longread is set to 2. |
| $chunksizeb | Long reads will be fragmented to this size in a second pass if $longread is set to 2. $chunksizeb should not equal $chunksizea and should not be a multiple of $chunksizea. |
| $fragseqfilestem | The base name for fasta files of fragmented reads. SLAG names the individual files as a concatenation of $fragseqfilestem, the cycle number, and ".fasta". |
| $longread | This variable controls the choice of long-read fragmentation, a method to coax assemblies out of long reads when read coverage is too low to run canu. <br/><br/>Allowed values are 1 for single-fragmentation, 2 for double fragmentation to two different lengths, and any other integer for no fragmentation. The default is zero. |


#### Variables specific to Racon polishing (error correction) of contigs
Note: Racon polishing has been implemented only for cap3 and phrap assembly.
The results are generally disappointing in polyploid species because reads get matched across different genomes.

| Variable | Explanation |
| ----------| ----------- |
| $useracon | Flag to polish or not to polish contigs with Racon. Meaningful values are 1 (use Racon) and 0 (do not use Racon). The default is 0. |
| $raconexe | Location of Racon, including path. To find this name on a personal filesystem, run "find / -name "Racon"". On a computing cluster with loadable modules, run "module avail" and search for "Racon" or "racon", run "module load Racon", and run "which racon". |
| $bowtie2dir | Directory where bowtie2 resides. To find this name on a personal filesystem, run "find / -name "bowtie2"" and copy the path. On a computing cluster with loadable modules, run "module avail" and search for "bowtie2", run "module load bowtie2", run "which bowtie2", and copy the path. |
| $readsfile | Name of a file of unpaired short reads in fastq format, which are used to polish contigs with Racon. This file must be prepared before starting SLAG. |
| $nrounds | Number of iterations to use Racon to polish contigs in each cycle. |
| $postem | Prefix for several file names generated during Racon polishing. The output files from Racon are named by concatenating $postem, "_cyc", the current SLAG cycle, "_polished_rd", the current round of Racon, and ".fasta". |
| $poblastdb | Base name of a blastable database of sequences from a reference genome. This database is used to compare polished contigs to reference sequence, and it is how the effect of polishing on contig accuracy was measured for the manuscript. The alignment output against this database is not used further within SLAG, but SLAG can execute the two following lines for convenient assessment of polishing effectiveness: print LOG "$blastdir/blastn -query $polfile -db $blastdb -out $blastout -evalue $evalue\n"; $blastdir/blastn -query $polfile -db $poblastdb -out $blastout -evalue $poevalue; These are currently lines 1049 and 1050 in SLAG.pl. Uncomment them if you want to run them. |
| $poevalue | The e-value for blastn of $polfile against $poblastdb. This has been set high (1e-10) where there were very few contigs after polishing. |


### Acknowledgement
I thank Lev Gorenstein of Purdue University Research and Academic Computing for writing and building SLAG.sif, the Singularity container for SLAG.pl.  I thank Melinda Crane for help with formatting this README.md.
