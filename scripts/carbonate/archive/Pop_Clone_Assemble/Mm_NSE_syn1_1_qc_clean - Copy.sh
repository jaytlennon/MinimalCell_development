#!/bin/bash

#$ -l vmem=50gb walltime=48:00:00

module load samtools; module load cutadapt; module load fastqc; module load bwa; module load java; module load picard
cd /N/dc2/scratch/rzmogerr/scratchy20190816/Mm_NSE/Sample_syn1_1

time cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -o Sample_syn1_1_R1_rmadapter.fastq -p Sample_syn1_1_R2_rmadapter.fastq Sample_syn1_1_R1.fastq Sample_syn1_1_R2.fastq

time cutadapt -q 15,10  -o Sample_syn1_1_R1_filtered.fastq -p Sample_syn1_1_R2_filtered.fastq Sample_syn1_1_R1_rmadapter.fastq Sample_syn1_1_R2_rmadapter.fastq

time cutadapt -u 15 -o Sample_syn1_1_R1_trimmed.fastq -p Sample_syn1_1_R2_trimmed.fastq Sample_syn1_1_R1_filtered.fastq Sample_syn1_1_R2_filtered.fastq

qsub -l walltime=20:00:00,vmem=64gb,nodes=1:ppn=4 Mm_NSE_syn1_1_Assemble.sh

mkdir fastqc

fastqc -o fastqc/ Sample_syn1_1_R1_trimmed.fastq

fastqc -o fastqc/ Sample_syn1_1_R2_trimmed.fastq

exit
