Instructions to install SLAG

1. Download SLAGball.tgz from Github.

2. Run "md5sum SLAGball.tgz" and compare the result to the file 
SLAGballmd5sum.txt.  Proceed if the checksums are identical.

3. Move SLAGball.tgz to the directory where you want to keep it.

4. Run "tar -xzf SLAGball.tgz" to expand the directory tree.  Its root, slagbase, contains 
SLAG.pl, SLAG.sif, and setupSLAGconfigurations.pl.  SLAG.sif is a Singularity container
that can run SLAG.pl.  Subdirectory templates contains 24 templates for configuration files, 
which are the way to communicate file names and settings to SLAG.pl or SLAG.sif.  

5. Subdirectory testcfg contains 24 populated configuration files and 24 shell scripts 
to validate a successful download and function on your system.  These shell scripts 
are meant to be executed from the testout subdirectory.  For example, running 
../testcfg/runCNTLRcanuSLAGtest.sh tests the container-included canu assembler 
with a test database and test query from subdirectory testdata, with output files 
going to testout.  Success is indicated by the presence of filenames with the suffix 
"subset.contigs" in testout or in a family of "*work" subdirectories of testout. 
Assemblers canu, Unicycler, and SPAdes produce "*work" subdirectories that must 
be deleted if you want to repeat the test.  The test scripts and configuration files are 
named by the convention LR = long read (2000 - 6000), SR = short read (2x150), 
MR = medium read (200 - 600), LR1 = long reads with single fragmentation, LR2 = 
long reads with double (overlapping) fragmentation, PI = Racon polishing, CNT = 
invocation of SLAG.sif (the Singularity container), and the assembler name.  Each 
shell script should take a few minutes, except for the PI scripts, which take up to 
one hour.  The test query is phenylalanine ammonia lyase from maize, and the test 
reads are simulated based on the interval from base 320000000 to base 
330000000 in wheat cv. 'Chinese Spring' chromosome 6D.  This interval contains 
one of the wheat copies of the phenylalanine ammonia lyase gene.
6.	If you are an administrator, you can add the path to slagbase to users' default 
$PATH, or include it in the set of loadable modules.
