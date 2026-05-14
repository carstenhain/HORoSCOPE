#!/bin/bash
#SBATCH -J dbg			        						# job name
#SBATCH --mem=256G                                      # memory pool for all cores
#SBATCH -n 16                                           # number of cores
#SBATCH --tmp=10G                                       # select nodes with this min space available (max=200)
#SBATCH -t 5-80:00:00                                   # runtime limit (D-HH:MM:SS)

THREADS=16
SID="/path/to/fastq/NAME"
KMER=61

lighter -t ${THREADS} -r ${SID}_1.fastq.gz -r ${SID}_2.fastq.gz -trim -discard -k 23 3100000000 0.188
ls -1 ${SID}_*.cor.fq.gz > ${SID}.list_reads
bcalm -in ${SID}.list_reads -kmer-size ${KMER} -abundance-min 3
rm ${SID}_*.cor.fq.gz ${SID}.list_reads
gzip ${SID}.unitigs.fa