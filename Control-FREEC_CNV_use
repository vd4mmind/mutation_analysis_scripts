#!/bin/sh
#
#$ -N C-FREEC_ova_samples
#$ -cwd
#$ -e err_C-FREEC_ova_samples.txt
#$ -o out_C-FREEC_ova_samples.txt
#$ -S /bin/sh
#$ -M vd4mmind@gmail.com
#$ -m bea
#$ -l h_vmem=25G

/path/freec -conf /path/config_file_FREEC.txt


#cat assess_significance.R | R --slave --args < *_CNVs > < *_ratio.txt >  typical usage of the significance script of CONTROL-FREEC

### High grade tumor
# creating significant genomic regions that have have a CNV
cat /path_to/vdas/client/exome_seq/test_Control_FREEC/scripts_FREEC/assess_significance.R | R --slave --args tumor.realigned.recal.bam_CNVs tumor.realigned.recal.bam_ratio.txt

#cat makeGraph.R | R --slave --args < ploidy > < *_ratio.txt > [< *_BAF.txt >]
#creating the CNV profile graph
cat /path_to/scripts_FREEC/makeGraph.R | R --slave --args 2 tumor.realigned.recal.bam_ratio.txt
