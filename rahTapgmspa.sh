#!/bin/sh
source /etc/profile
ml bioinfo aTRAM
#atram.py --blast-db=path_to_atram_library/LIBRARY_PREFIX --query=path_to_reference_loci/Locus.fasta --assembler=ASSEMBLER_CHOICE --output-prefix=path_to_output/OUTPUT_PREFIX
echo "atram.py --blast-db=/scratch/brown/ycrane/hemiSRR9125476_htsclepurrenrefforatram2 --query=Zeaphosphoglucomutases.fasta --assembler=spades --output-prefix=/scratch/brown/ycrane/atramout/hemiStanleyphosphoglucomutasesbyatram2 --iterations=21 --protein --cpus=10 --log-file=hemiStanleyphosphoglucomutasesbyatram2log.txt --temp-dir=/scratch/brown/ycrane/atramtemp --timeout=100000 --no-filter --evalue=1e-10 --max-memory=95"
date
atram.py --blast-db=/scratch/brown/ycrane/hemiSRR9125476_htsclepurrenrefforatram2 --query=Zeaphosphoglucomutases.fasta --assembler=spades --output-prefix=/scratch/brown/ycrane/atramout/hemiStanleyphosphoglucomutasesbyatram2 --iterations=21 --protein --cpus=10 --log-file=hemiStanleyphosphoglucomutasesbyatram2log.txt --temp-dir=/scratch/brown/ycrane/atramtemp --timeout=100000 --no-filter --evalue=1e-10 --max-memory=95
date
exit
