#!/bin/bash

#$ -q rcc-30d

# SamMuri_2Y
date
cd /N/dc2/projects/muri2/SamMuri/Year_2
AR=(<Sample Names>)		!!!So, the way u make arrays in UNIX---https://stackoverflow.com/questions/1878882/arrays-in-unix-shell

for i in "${AR[@]}"		!!!And this command is just how you do a for loop for an array in UNIX. Easy peasy.
do
	cd Sample_${i}		!!!So you are gonna want to start by putting all of the R1 R2 pairs, for each sample, into its own carefully named directory. Then youre gonna loop thru the carefully named directories using cd	
	gunzip *
	cat *R1_001.fastq > Sample_${i}_R1.fastq
	cat *R2_001.fastq > Sample_${i}_R2.fastq
        	
	# clean and trim		https://cutadapt.readthedocs.io/en/v1.9.1/guide.html
	echo "#!/bin/bash" > SamMuri_2Y${i}_qc_clean.sh
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "#$ -l vmem=50gb walltime=48:00:00" >> SamMuri_2Y${i}_qc_clean.sh
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "cd /N/dc2/projects/muri2/SamMuri/Year_2/Sample_${i}" >> SamMuri_2Y${i}_qc_clean.sh
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "time cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -o Sample_${i}_R1_rmadapter.fastq -p Sample_${i}_R2_rmadapter.fastq Sample_${i}_R1.fastq Sample_${i}_R2.fastq" >> SamMuri_2Y${i}_qc_clean.sh			!!!TAATANTGCTGCGTTTTCTGAGAATTTAATATAAGTTC 		-a -A is for fwd and rev adapter for paired-end reads.		-o is output.	-p is --paired-output			This hard-coded adapter seqq is just the first 13 bases of the known illumina adapters---see pp 24-25 here: https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-11.pdf		You only ned the first 13? Ja, that is based on reply number 4 on this thread from SEQanswers http://seqanswers.com/forums/archive/index.php/t-17168.html
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "time cutadapt -q 15,10  -o Sample_${i}_R1_filtered.fastq -p Sample_${i}_R2_filtered.fastq Sample_${i}_R1_rmadapter.fastq Sample_${i}_R2_rmadapter.fastq" >> SamMuri_2Y${i}_qc_clean.sh		!!!-q is using cutadapt to also trim low quality ends from reads.
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "time cutadapt -u 15 -o Sample_${i}_R1_trimmed.fastq -p Sample_${i}_R2_trimmed.fastq Sample_${i}_R1_filtered.fastq Sample_${i}_R2_filtered.fastq" >> SamMuri_2Y${i}_qc_clean.sh	!!!-u is removing a fixed number of bases--- 15 in this case
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "qsub -l walltime=20:00:00,vmem=64gb,nodes=1:ppn=4 SamMuri_2Y${i}_Assemble.sh" >> SamMuri_2Y${i}_qc_clean.sh
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "mkdir fastqc" >> SamMuri_2Y${i}_qc_clean.sh
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "fastqc -o fastqc/ Sample_${i}_R1_trimmed.fastq" >> SamMuri_2Y${i}_qc_clean.sh
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "fastqc -o fastqc/ Sample_${i}_R2_trimmed.fastq" >> SamMuri_2Y${i}_qc_clean.sh
	echo "" >> SamMuri_2Y${i}_qc_clean.sh
	echo "exit" >> SamMuri_2Y${i}_qc_clean.sh

#### Assembly ####

	echo "#!/bin/bash" > SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "cd /N/dc2/projects/muri2/SamMuri/Year_2/Sample_${i}/" >> SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "echo \"BWA\" >&2" >> SamMuri_2Y${i}_Assemble.sh
	echo "time bwa mem /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna Sample_${i}_R1_trimmed.fastq Sample_${i}_R2_trimmed.fastq > Sample_${i}.sam" >> SamMuri_2Y${i}_Assemble.sh	!!!so here you need to get your own reference genomes. This means you are gonna chop up into directories by strain. Reference genomes: 8: NCBIs syn1.0	16: NCBIs syn3.0 OR BETTER YET, ask John Glass if they have a better genome. 1-7: 8. 9-15: 16			Or simplify your life and do: 1-8: NCBI syn1.0. 9-16: John Glasss syn3B.
   	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "qsub -l walltime=2:00:00,vmem=20gb SamMuri_2Y${i}_compress.fastq.sh" >> SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "echo \"SAMTOOLS\" >&2" >> SamMuri_2Y${i}_Assemble.sh
	echo "samtools view -bS -T /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna Sample_${i}.sam > Sample_${i}.bam" >> SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "samtools sort Sample_${i}.bam Sample_${i}.sorted" >> SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "samtools index Sample_${i}.sorted.bam" >> SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	
	
        ### GATK ###	
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=Sample_${i}.sorted.bam O=Sample_${i}.sorted.fixed.bam SORT_ORDER=coordinate RGID=GeneConv RGLB=bar RGPL=illumina RGSM=Sample_${i} RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> SamMuri_2Y${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=Sample_${i}.sorted.fixed.bam O=Sample_${i}.sorted.fixed.marked.bam M=Sample_${i}.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> SamMuri_2Y${i}_Assemble.sh
	echo "GATK" >&2
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -I Sample_${i}.sorted.fixed.marked.bam -T RealignerTargetCreator -o Sample_${i}.intervals" >> SamMuri_2Y${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -I Sample_${i}.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals Sample_${i}.intervals --filter_bases_not_stored -o Sample_${i}.sorted.fixed.marked.realigned.bam" >> SamMuri_2Y${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -I Sample_${i}.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o Sample_${i}.sorted.fixed.marked.realigned.vcf" >> SamMuri_2Y${i}_Assemble.sh!!! Megan says you want -glm set to BOTH.		If you are doing clones, you want -ploidy 1---its default avlue is -ploidy 2 which is what you want for popns, it allows for 2 alleles to be called at a site
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I Sample_${i}.sorted.fixed.marked.realigned.bam -R /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -rf BadCigar --filter_bases_not_stored -knownSites Sample_${i}.sorted.fixed.marked.realigned.vcf -o Sample_${i}.recal_data.grp" >> SamMuri_2Y${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -I Sample_${i}.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o Sample_${i}.mapped.bam -BQSR Sample_${i}.recal_data.grp" >> SamMuri_2Y${i}_Assemble.sh
	echo "" >> SamMuri_2Y${i}_Assemble.sh
	echo "exit" >> SamMuri_2Y${i}_Assemble.sh
	
 	# compress.fastq again
	echo "#!/bin/bash" > SamMuri_2Y${i}_compress.fastq.sh
	echo "" >> SamMuri_2Y${i}_compress.fastq.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> SamMuri_2Y${i}_compress.fastq.sh
	echo "" >> SamMuri_2Y${i}_compress.fastq.sh
	echo "cd /N/dc2/projects/muri2/SamMuri/Year_2/Sample_${i}/" >> SamMuri_2Y${i}_compress.fastq.sh
	echo "rm Sample_${i}_R*_filtered.fastq" >> SamMuri_2Y${i}_compress.fastq.sh
	echo "rm Sample_${i}_R*_rmadapter.fastq" >> SamMuri_2Y${i}_compress.fastq.sh
	echo "" >> SamMuri_2Y${i}_compress.fastq.sh
	echo "exit" >> SamMuri_2Y${i}_compress.fastq.sh
	
	
		
	chmod u+x SamMuri_2Y${i}_qc_clean.sh
	chmod u+x SamMuri_2Y${i}_Assemble.sh
	chmod u+x SamMuri_2Y${i}_compress.fastq.sh
		
	qsub -l walltime=2:00:00,vmem=20gb SamMuri_2Y${i}_qc_clean.sh
	cd ..
	
	done
	
	
	ok, here is a back of copy pasted from the attempted GATK run. it was not working.
	
	        ### GATK ###	
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=Sample_${i}.sorted.bam O=Sample_${i}.sorted.fixed.bam SORT_ORDER=coordinate RGID=GeneConv RGLB=bar RGPL=illumina RGSM=Sample_${i} RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=Sample_${i}.sorted.fixed.bam O=Sample_${i}.sorted.fixed.marked.bam M=Sample_${i}.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> Mm_NSE_${i}_Assemble.sh
	echo "GATK" >&2
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I Sample_${i}.sorted.fixed.marked.bam -T RealignerTargetCreator -o Sample_${i}.intervals" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I Sample_${i}.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals Sample_${i}.intervals --filter_bases_not_stored -o Sample_${i}.sorted.fixed.marked.realigned.bam" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I Sample_${i}.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o Sample_${i}.sorted.fixed.marked.realigned.vcf" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I Sample_${i}.sorted.fixed.marked.realigned.bam -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -rf BadCigar --filter_bases_not_stored -knownSites Sample_${i}.sorted.fixed.marked.realigned.vcf -o Sample_${i}.recal_data.grp" >> Mm_NSE_${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /gpfs/home/r/z/rzmogerr/Carbonate/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.fasta -I Sample_${i}.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o Sample_${i}.mapped.bam -BQSR Sample_${i}.recal_data.grp" >> Mm_NSE_${i}_Assemble.sh
	echo "" >> Mm_NSE_${i}_Assemble.sh
	echo "exit" >> Mm_NSE_${i}_Assemble.sh