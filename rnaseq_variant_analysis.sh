#!/bin/sh
#
# up is a bash function to simplify vertical file system navigation.
# Copyright (C) 2016 Vivek Das
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
# AUTHOR: Vivek Das.
# RNA-Seq variant calling pieline accoring to GATK Best practices.
# https://www.broadinstitute.org/gatk/guide/article?id=3891
#
# Call with following arguments
# sh rna_seq_variant_pipeline.sh  <output_basename> <fastq folder> <output_folder_loc> [cpus]
#
# Putting the path of STAR aligner
#
# GATK bundle set : one can obtain these from gatk ftp (knonw as gatk bundle)
# ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
#

bn=$1
floc=$2
outloc=$3
nproc=${4:-12}

#Path to reference genome and Index files
star_ref="/path_to/genomes/Homo_sapiens/UCSC_spiked/hg19/Sequence/STARindexes"
ref="/path_to/genomes/Homo_sapiens/UCSC_spiked/hg19/Sequence/genome.fa"

#Path to STAR,java,gatk and picard tools
java_home="/home/vdas/tools/jdk1.8.0_77/bin/java"
STAR="/path_to/softwares/STAR-STAR_2.4.1c/bin/Linux_x86_64_static/STAR"
gatk="/path_to/softwares/GATK-3.5/GenomeAnalysisTK.jar"
picard="/home/vdas/tools/picard-tools-2.2.0/picard.jar"
sambamba="/home/vdas/tools/sambamba/build/sambamba"

#Path to gatk bundle set files
millsIndels="/path_to/GATK_2.8_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
KGIndels="/path_to/GATK_2.8_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf"
dbSNP138="/path_to/GATK_2.8_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf"


#Create an output directory
mkdir -p "$outloc/${bn}_processed"
#fout="$outloc/${bn}_processed"
echo "output directory for fastq $outloc/${bn}"_processed" ..."

#fout=$outloc/$bn"_processed"
fout="$outloc/${bn}_processed"
echo "$fout ..."

echo "performing assembly to create one fastq file for each read mates ..."
zcat $floc/*R1*.fastq.gz > $fout/${bn}_R1.fastq
zcat $floc/*R2*.fastq.gz > $fout/${bn}_R2.fastq

#STAR 2 pass basic mode run
echo -e "["$(date)"]\tAligning.."
$STAR --outFileNamePrefix $fout/${bn} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMattrRGline ID:${bn} CN:RSEQ_GT_LAB LB:PairedEnd PL:Illumina PU:C7MC6ACXX SM:$bn --genomeDir $star_ref --runThreadN 12 --readFilesIn $fout/${bn}_R1.fastq $fout/${bn}_R2.fastq --twopassMode Basic

#samtools sort
echo -e "["$(date)"]\tSorting.."
$sambamba sort -o $fout/${bn}"_sorted.bam" -p -t 12 $fout/${bn}"Aligned.out.bam"

#rm $fout/${bn}"Aligned.out.bam"
rm $fout/${bn}_R1.fastq
rm $fout/${bn}_R2.fastq

#samtools index
echo -e "["$(date)"]\tIndexing.."
$sambamba index -p -t 12 $fout/${bn}"_sorted.bam"

#picard mark duplicates
echo -e "["$(date)"]\tMarking duplicates.."
$java_home -Xmx2g -jar $picard MarkDuplicates I=$fout/${bn}"_sorted.bam" O=$fout/${bn}"_dupMarked.bam" M=$fout/${bn}"_dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 2>$fout/${bn}.MarkDuplicates.log

rm $fout/${bn}"Aligned.out.bam"
rm $fout/${bn}"_sorted.bam"
rm $fout/${bn}"_sorted.bam.bai"

#SplitNCigarReads
echo -e "["$(date)"]\tSpliting reads.."
$java_home -Xmx2g -jar $gatk -T SplitNCigarReads -R $ref -I $fout/${bn}"_dupMarked.bam" -o $fout/${bn}"_split.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2>$fout/${bn}.SplitNCigarReads.log

rm $fout/${bn}"_dupMarked.bam"
rm $fout/${bn}"_dupMarked.bai"

#Create targets for indel realignment
echo -e "["$(date)"]\tCreating targets for indel realignment.."
$java_home -Xmx2g -jar $gatk -T RealignerTargetCreator -R $ref -I $fout/${bn}"_split.bam" -o $fout/${bn}".intervals" -nt 12 -known $millsIndels -known $KGIndels 2>$fout/${bn}.indel.log

#Perform indel realignment
echo -e "["$(date)"]\tPerforming Indel Realignment.."
$java_home -Xmx2g-jar $gatk -T IndelRealigner -R $ref -I $fout/${bn}"_split.bam" -targetIntervals $fout/${bn}".intervals" -known $millsIndels -known $KGIndels -o $fout/${bn}"realigned.bam" 2>$fout/${bn}.indel2.log

rm $fout/${bn}"_split.bam"
rm $fout/${bn}"_split.bai"

#Perform BQSR
echo -e "["$(date)"]\tPerforming BQSR.."
$java_home -Xmx2g -jar $gatk -T BaseRecalibrator -I $fout/${bn}"_realigned.bam" -R $ref -knownSites $KGIndels -knownSites $millsIndels -knownSites $dbSNP138 -o $fout/${bn}"_recal.table" 2>$fout/${bn}.BQSR.log

#Print recalibrated reads
echo -e "["$(date)"]\tPrinting recalibrated reads.."
$java_home -Xmx2g -jar $gatk -T PrintReads -R $ref -I $fout/${bn}"_realigned.bam" -nct 12 -BQSR $fout/${bn}"_recal.table" -o $fout/${bn}"_recal.bam" 2>$fout/${bn}.BQSR2.log

#$java_home -Xmx2g -jar $gatk -T PrintReads -R $genomeFasta/genome.fa -I $dedup_split_bam/processed_indelrealign.bam -nct 50 -BQSR $dedup_split_bam/BQSR.table -o $dedup_split_bam/BQSR_recal.bam

rm $fout/${bn}"_realigned.bam"
rm $fout/${bn}"_realigned.bai"

#Run HaplotypeCaller
echo -e "["$(date)"]\tRunning HaplotypeCaller.."
$java_home -Xmx2g -jar $gatk -T HaplotypeCaller -R $ref -I $fout/${bn}"_recal.bam" -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $fout/${bn}".vcf" 2>$fout/${bn}.HaplotypeCaller.log

#Filter variants for SNPs and then filter for qaulity scores
echo -e "["$(date)"]\tFiltering Point Variants.."

$java_home -Xmx2g -jar $gatk -T SelectVariants -R $ref -V $fout/${bn}".vcf" -selectType SNP -o $fout/${bn}"_raw_snps.vcf"

#### filter over snps for scores

$java_home -Xmx2g -jar $gatk -T VariantFiltration -R $ref -V $fout/${bn}"_raw_snps.vcf" --filterExpression "DP < 5 || SB > -0.1 || QUAL < 50.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > 500" --filterName "DAS_snp_filters" -o $fout/${bn}"_filtered_snps.vcf" 2>$fout/${bn}.VariantFilter_snps.log

echo -e "["$(date)"]\tFiltering INDEL Variants.."
###select indels
$java_home -Xmx2g -jar $gatk -T SelectVariants -R $ref -V $fout/${bn}".vcf" -selectType INDEL -o $fout/${bn}"_raw_indels.vcf"

#### filter over indel files
$java_home -Xmx2g -jar $gatk -T VariantFiltration -R $ref -V $fout/${bn}"_raw_indels.vcf" --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filterName "DAS_indel_filters" -o $fout/${bn}"_filtered_indels.vcf" 2>$fout/${bn}.VariantFilter_indels.log

echo -e "["$(date)"]\tDONE!"
