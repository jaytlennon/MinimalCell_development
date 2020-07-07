#!/bin/bash

cd /N/dc2/projects/muri2/roy/MA_s1/

#AR=(3B_U 3B_V 3B_W 3B_X 3B_Y)

#for i in "${AR[@]}"
for ((i=1;i<=96;i++));
do
	module load breseq
	echo "cd /N/dc2/projects/muri2/roy/MA_s1/" > VPN_s1_BreSeq_${i}
	echo "" >> VPN_s1_BreSeq_${i}
	echo "module load breseq" >> VPN_s1_BreSeq_${i}
	echo "" >> VPN_s1_BreSeq_${i}
	#echo "breseq -j 8 -p -o s1_${i}_breseq -r /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn3A_reference_CP016816.2/Synthetic.bacterium_JCVI-Syn3A.gb Sample_${i}/Sample_${i}_R1_trimmed.fastq" >> s1_BreSeq_${i}
	echo "breseq -j 8 -p -o VPN_s1_${i}_breseq -r /N/dc2/projects/muri2/roy/Mm_NSE/reference/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb s1_${i}/s1_${i}_R1_trimmed.fastq s1_${i}/s1_${i}_R2_trimmed.fastq" >> VPN_s1_BreSeq_${i}
	qsub -l walltime=16:00:00,vmem=100gb,nodes=1:ppn=8 ./VPN_s1_BreSeq_${i}
	done