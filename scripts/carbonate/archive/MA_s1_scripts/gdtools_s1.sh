#!/bin/bash

cd /N/dc2/projects/muri2/roy/MA_s1/
#mkdir MA_s1_outsub/
#AR=(3B_U 3B_V 3B_W 3B_X 3B_Y)

#for i in "${AR[@]}"
for ((i=1;i<=96;i++));
do
	module load breseq
	gdtools subtract -o ../MA_s1_outsub/s1_${i}_sub_annotated.gd s1_${i}_breseq/output/evidence/annotated.gd MA_s1_anc_annotated.gd
	gdtools GD2VCF -o ../MA_s1_GD2VCF/s1_${i}_clean.VCF ../MA_s1_outsub/s1_${i}_sub_annotated.gd -r /N/dc2/projects/muri2/roy_20200302/Mm_NSE/reference/JCVI-syn1.0_reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb
	done