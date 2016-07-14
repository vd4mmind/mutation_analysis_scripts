#!/bin/sh
#
# MODIFIED: Vivek Das.
# Exome-seq variant calling pieline accoring to GATK Best practices.
# https://www.broadinstitute.org/gatk/guide/article?id=3891
#
# Call with following arguments
# sh exome_auto_gatk_v2..sh  <output_basename> <fastq folder> <output_folder_loc> [cpus]
#
# Putting the bwa mem
#
# GATK bundle set : one can obtain these from gatk ftp (knonw as gatk bundle)
# ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
#

bn=$1
floc=$2
outloc=$3
nproc=${4:-12}

#Path to reference genome and Index files
#star_ref="/scratch/GT/genomes/Homo_sapiens/UCSC_spiked/hg19/Sequence/STARindexes"
bwa_idx="/path/vdas/test_exome/exome/hg19"
ref="/path/vdas/test_exome/exome/hg19.fa"

#Path to STAR,java,gatk and picard tools
java_home="/home/vdas/tools/jdk1.8.0_77/bin/java"
#STAR="/scratch/GT/softwares/STAR-STAR_2.4.1c/bin/Linux_x86_64_static/STAR"
bwamem="/path/softwares/bwa/bwa"
gatk="/path/softwares/GATK-3.5/GenomeAnalysisTK.jar"
picard="/home/vdas/tools/picard-tools-2.2.0/picard.jar"
sambamba="/home/vdas/tools/sambamba/build/sambamba"

#Path to gatk bundle set files
millsIndels="/path/GATK_2.8_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
KGIndels="/path/GATK_2.8_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf"
dbSNP138="/path/GATK_2.8_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf"
#Exon target baits from Agilent SSv4 set (should be upgraded to new version if new agilent kit is used)
ss4exonbaits="/path/vdas/referenceBed/hg19/ss_v4/Exon_SSV4_clean.bed"

#Create an output directory
mkdir -p "$outloc/${bn}_processed"
#fout="$outloc/${bn}_processed"
echo "output directory for fastq $outloc/${bn}"_processed" ..."

#fout=$outloc/$bn"_processed"
fout="$outloc/${bn}_processed"
echo "$fout ..."

echo "performing assembly to create one fastq file for each read mates ..."
zcat $floc/*1.fastq.gz > $fout/${bn}_R1.fastq
zcat $floc/*2.fastq.gz > $fout/${bn}_R2.fastq

#STAR 2 pass basic mode run
#echo -e "["$(date)"]\tAligning.."
#$STAR --outFileNamePrefix $fout/${bn} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMattrRGline ID:${bn} CN:RSEQ_GT_LAB LB:PairedEnd PL:Illumina PU:C7MC6ACXX SM:$bn --genomeDir $star_ref --runThreadN 12 --readFilesIn $fout/${bn}_R1.fastq $fout/${bn}_R2.fastq --twopassMode Basic


# Aligning with bwa mem
echo -e "["$(date)"]\tAligning.."
$bwamem mem -t 8 -M -R "@RG\tID:${bn}\tLB:PairedEnd\tPL:Illumina\tPU:000000000-A442D\tSM:${bn}" $bwa_idx $fout/${bn}_R1.fastq $fout/${bn}_R2.fastq | samtools view -bS -F 0x04 - > $fout/${bn}"_Aligned.out.bam"

#samtools sort
echo -e "["$(date)"]\tSorting.."
$sambamba sort -o $fout/${bn}"_sorted.bam" -p -t 12 $fout/${bn}"_Aligned.out.bam"

#rm $fout/${bn}"_Aligned.out.bam"
rm $fout/${bn}_R1.fastq
rm $fout/${bn}_R2.fastq

#samtools index
echo -e "["$(date)"]\tIndexing.."
$sambamba index -p -t 12 $fout/${bn}"_sorted.bam"

#picard mark duplicates
echo -e "["$(date)"]\tMarking duplicates.."
$java_home -Xmx2g -jar $picard MarkDuplicates I=$fout/${bn}"_sorted.bam" O=$fout/${bn}"_dupMarked.bam" M=$fout/${bn}"_dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 2>$fout/${bn}.MarkDuplicates.log

rm $fout/${bn}"_Aligned.out.bam"
rm $fout/${bn}"_sorted.bam"
rm $fout/${bn}"_sorted.bam.bai"

#SplitNCigarReads
#echo -e "["$(date)"]\tSpliting reads.."
#$java_home -Xmx2g -jar $gatk -T SplitNCigarReads -R $ref -I $fout/${bn}"_dupMarked.bam" -o $fout/${bn}"_split.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2>$fout/${bn}.SplitNCigarReads.log


#Create targets for indel realignment
echo -e "["$(date)"]\tCreating targets for indel realignment.."
$java_home -Xmx2g -jar $gatk -T RealignerTargetCreator -R $ref -I $fout/${bn}"_dupMarked.bam" -o $fout/${bn}".intervals" -nt 12 -known $millsIndels -known $KGIndels -L $ss4exonbaits 2>$fout/${bn}.indel.log

#Perform indel realignment
echo -e "["$(date)"]\tPerforming Indel Realignment.."
$java_home -Xmx2g -jar $gatk -T IndelRealigner -R $ref -I $fout/${bn}"_dupMarked.bam" -targetIntervals $fout/${bn}".intervals" -known $millsIndels -known $KGIndels -o $fout/${bn}"_realigned.bam" 2>$fout/${bn}.indel2.log

#rm $fout/${bn}"_split.bam"
#rm $fout/${bn}"_split.bai"

rm $fout/${bn}"_dupMarked.bam"
rm $fout/${bn}"_dupMarked.bai"

#Perform BQSR
echo -e "["$(date)"]\tPerforming BQSR.."
$java_home -Xmx2g -jar $gatk -T BaseRecalibrator -I $fout/${bn}"_realigned.bam" -R $ref -knownSites $KGIndels -knownSites $millsIndels -knownSites $dbSNP138 -o $fout/${bn}"_recal.table" -L $ss4exonbaits 2>$fout/${bn}.BQSR.log

#Print recalibrated reads
echo -e "["$(date)"]\tPrinting recalibrated reads.."
$java_home -Xmx2g -jar $gatk -T PrintReads -R $ref -I $fout/${bn}"_realigned.bam" -nct 12 -BQSR $fout/${bn}"_recal.table" -o $fout/${bn}"_recal.bam" 2>$fout/${bn}.BQSR2.log

#$java_home -Xmx2g -jar $gatk -T PrintReads -R $genomeFasta/genome.fa -I $dedup_split_bam/processed_indelrealign.bam -nct 50 -BQSR $dedup_split_bam/BQSR.table -o $dedup_split_bam/BQSR_recal.bam

#rm $fout/${bn}"_realigned.bam"
#rm $fout/${bn}"_realigned.bai"

#Run HaplotypeCaller
echo -e "["$(date)"]\tRunning HaplotypeCaller.."
$java_home -Xmx2g -jar $gatk -T HaplotypeCaller -R $ref -I $fout/${bn}"_recal.bam" -L $ss4exonbaits -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $fout/${bn}".vcf" 2>$fout/${bn}.HaplotypeCaller.log

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
