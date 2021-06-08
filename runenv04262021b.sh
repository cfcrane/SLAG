#!/bin/bash
source /etc/profile
#ml openmpi bioinfo blast Vmatch SOAPdenovo/240 GenomeThreader
ml bioinfo blast Vmatch SOAPdenovo/240 GenomeThreader
cd /scratch/brown/ycrane
srun --mpi=pmi2 -n 10 /scratch/brown/ycrane/SRAssembler_MPI -q /home/ycrane/Zeahexokinases.fasta -t protein -T /scratch/brown/ycrane/sratempZmhext2 -o /scratch/brown/ycrane/sraoutZmhext2 -r /scratch/brown/ycrane/sraout/processed_reads/ERR328821xcleaned_1library -1 /scratch/brown/ycrane/ERR328821xcleaned_1.fastq -2 /scratch/brown/ycrane/ERR328821xcleaned_2.fastq -Z 400 -A 0 -n 21 -b 1 -x 0 -z 1 -d 500 -i 1000 -m 1000 -M 200000 -p /home/ycrane/testSRAssembler4.conf -s maize
