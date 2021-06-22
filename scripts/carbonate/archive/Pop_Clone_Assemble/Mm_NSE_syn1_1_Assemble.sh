#!/bin/bash

#$ -l vmem=50gb walltime=48:00:00 

module load samtools; module load cutadapt; module load fastqc; module load bwa
cd /N/dc2/scratch/rzmogerr/scratchy20190802/Mm_NSE/Sample_syn1_1/

echo "BWA" >&2
time bwa mem /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta Sample_syn1_1_R1_trimmed.fastq Sample_syn1_1_R2_trimmed.fastq > Sample_syn1_1.sam

qsub -l walltime=2:00:00,vmem=20gb Mm_NSE_syn1_1_compress.fastq.sh

echo "SAMTOOLS" >&2
samtools view -bS -T /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta Sample_syn1_1.sam > Sample_syn1_1.bam

samtools sort Sample_syn1_1.bam Sample_syn1_1.sorted

samtools index Sample_syn1_1.sorted.bam


java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=Sample_syn1_1.sorted.bam O=Sample_syn1_1.sorted.fixed.bam SORT_ORDER=coordinate RGID=GeneConv RGLB=bar RGPL=illumina RGSM=Sample_syn1_1 RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=Sample_syn1_1.sorted.fixed.bam O=Sample_syn1_1.sorted.fixed.marked.bam M=Sample_syn1_1.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I Sample_syn1_1.sorted.fixed.marked.bam -T RealignerTargetCreator -o Sample_syn1_1.intervals
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I Sample_syn1_1.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals Sample_syn1_1.intervals --filter_bases_not_stored -o Sample_syn1_1.sorted.fixed.marked.realigned.bam
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I Sample_syn1_1.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o Sample_syn1_1.sorted.fixed.marked.realigned.vcf
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I Sample_syn1_1.sorted.fixed.marked.realigned.bam -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -rf BadCigar --filter_bases_not_stored -knownSites Sample_syn1_1.sorted.fixed.marked.realigned.vcf -o Sample_syn1_1.recal_data.grp
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I Sample_syn1_1.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o Sample_syn1_1.mapped.bam -BQSR Sample_syn1_1.recal_data.grp

exit
