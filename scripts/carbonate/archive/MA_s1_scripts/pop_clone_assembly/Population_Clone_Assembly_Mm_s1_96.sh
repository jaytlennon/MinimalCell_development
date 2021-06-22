#!/bin/bash

#$ -q rcc-30d

# SamMuri_2Y
date
cd /N/dc2/projects/muri2/roy/MA_s1/
#AR=(96)

#for i in "${AR[@]}"
for ((i=1;i<=96;i++));
do
    module load gatk/3.8; module load samtools; module load cutadapt; module load fastqc; module load bwa; module load picard
	cd s1_${i}
	cat *R1_001.fastq > s1_${i}_R1.fastq
	cat *R2_001.fastq > s1_${i}_R2.fastq
	
	# clean and trim
	echo "#!/bin/bash" > MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "#$ -l vmem=50gb walltime=2:00:00" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa; module load java; module load picard" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "cd /N/dc2/projects/muri2/roy/MA_s1/s1_${i}" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "time cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o s1_${i}_R1_rmadapter.fastq -p s1_${i}_R2_rmadapter.fastq s1_${i}_R1.fastq s1_${i}_R2.fastq" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "time cutadapt -u 5 -o s1_${i}_R1_trimmed.fastq s1_${i}_R1_rmadapter.fastq" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "time cutadapt -u 5 -o s1_${i}_R2_trimmed.fastq s1_${i}_R2_rmadapter.fastq" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "qsub -l walltime=2:00:00,vmem=64gb,nodes=1:ppn=4 MA_s1_${i}_Assemble.sh" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "mkdir fastqc" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "fastqc -o fastqc/ s1_${i}_R1_trimmed.fastq" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "fastqc -o fastqc/ s1_${i}_R2_trimmed.fastq" >> MA_s1_${i}_qc_clean.sh
	echo "" >> MA_s1_${i}_qc_clean.sh
	echo "exit" >> MA_s1_${i}_qc_clean.sh
	
#### Assembly ####

	echo "#!/bin/bash" > MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "#$ -l vmem=50gb walltime=24:00:00 " >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "cd /N/dc2/projects/muri2/roy/MA_s1/s1_${i}/" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "echo \"BWA\" >&2" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "time bwa aln /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta s1_${i}_R1_trimmed.fastq > s1_${i}_R1_aln.sai" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "time bwa aln /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta s1_${i}_R2_trimmed.fastq > s1_${i}_R2_aln.sai" >> MA_s1_${i}_Assemble.sh
    echo "" >> MA_s1_${i}_Assemble.sh
	echo "time bwa sampe /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta s1_${i}_R1_aln.sai s1_${i}_R2_aln.sai s1_${i}_R1_trimmed.fastq s1_${i}_R2_trimmed.fastq > s1_${i}.sam" >> MA_s1_${i}_Assemble.sh
   	echo "" >> MA_s1_${i}_Assemble.sh
	echo "qsub -l walltime=1:00:00,vmem=20gb MA_s1_${i}_compress.fastq.sh" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "echo \"SAMTOOLS\" >&2" >> MA_s1_${i}_Assemble.sh
	echo "samtools view -bS -T /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta s1_${i}.sam > s1_${i}.bam" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "samtools sort -o s1_${i}.sorted.bam s1_${i}.bam" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "samtools index s1_${i}.sorted.bam" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	
	        ### GATK ###	
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=s1_${i}.sorted.bam O=s1_${i}.sorted.fixed.bam SORT_ORDER=coordinate RGID=GeneConv RGLB=bar RGPL=illumina RGSM=s1_${i} RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> MA_s1_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=s1_${i}.sorted.fixed.bam O=s1_${i}.sorted.fixed.marked.bam M=s1_${i}.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> MA_s1_${i}_Assemble.sh
	echo "GATK" >&2
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I s1_${i}.sorted.fixed.marked.bam -T RealignerTargetCreator -o s1_${i}.intervals" >> MA_s1_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I s1_${i}.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals s1_${i}.intervals --filter_bases_not_stored -o s1_${i}.sorted.fixed.marked.realigned.bam" >> MA_s1_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I s1_${i}.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o s1_${i}.sorted.fixed.marked.realigned.vcf" >> MA_s1_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I s1_${i}.sorted.fixed.marked.realigned.bam -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -rf BadCigar --filter_bases_not_stored -knownSites s1_${i}.sorted.fixed.marked.realigned.vcf -o s1_${i}.recal_data.grp" >> MA_s1_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I s1_${i}.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o s1_${i}.mapped.bam -BQSR s1_${i}.recal_data.grp" >> MA_s1_${i}_Assemble.sh
	echo "" >> MA_s1_${i}_Assemble.sh
	echo "exit" >> MA_s1_${i}_Assemble.sh	
	
 	# compress.fastq again
	echo "#!/bin/bash" > MA_s1_${i}_compress.fastq.sh
	echo "" >> MA_s1_${i}_compress.fastq.sh
	echo "#$ -l vmem=50gb walltime=4:00:00 " >> MA_s1_${i}_compress.fastq.sh
	echo "module load samtools; module load cutadapt; module load fastqc; module load bwa; module load java; module load picard" >> MA_s1_${i}_compress.fastq.sh
	echo ""
	echo "" >> MA_s1_${i}_compress.fastq.sh
	echo "cd /N/dc2/projects/muri2/roy/MA_s1/s1_${i}/" >> MA_s1_${i}_compress.fastq.sh
	#echo "rm s1_${i}_R*_filtered.fastq" >> MA_s1_${i}_compress.fastq.sh
	echo "rm s1_${i}_R*_rmadapter.fastq" >> MA_s1_${i}_compress.fastq.sh
	echo "" >> MA_s1_${i}_compress.fastq.sh
	echo "exit" >> MA_s1_${i}_compress.fastq.sh
	
	
		
	chmod u+x MA_s1_${i}_qc_clean.sh
	chmod u+x MA_s1_${i}_Assemble.sh
	chmod u+x MA_s1_${i}_compress.fastq.sh

		
	qsub -l walltime=2:00:00,vmem=20gb MA_s1_${i}_qc_clean.sh
	cd ..
	done