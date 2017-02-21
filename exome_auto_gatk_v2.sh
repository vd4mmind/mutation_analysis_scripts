#!/bin/sh
#
# MODIFIED: Vivek Das.
# Exome-seq variant calling pieline accoring to GATK Best practices.
# https://www.broadinstitute.org/gatk/guide/article?id=3891
#
# Call with following arguments
# sh exome_auto_gatk_v2..sh  <output_basename> <fastq folder> <output_folder_loc> [cpus]
# you can run the script in above mentioned way through another processing script that will log the processing with time at each step
# Putting the bwa mem
#
# GATK bundle set : one can obtain these from gatk ftp (knonw as gatk bundle)
# ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
# The bundles were downloaded from the GATK 2.8 while the GATK 3.5 was used for the analysis
# Please download all the required tools and install them in your HPC and modify the paths accordingly for effective running of the script
# Initializing the variables required for calling the automated script
bn=$1
floc=$2
outloc=$3
nproc=${4:-12}

#Path to reference genome and Index files (This code is based on hg19)
#mention the location of the indexed hg19 file
bwa_idx="/path/vdas/test_exome/exome/hg19"
#mention the location of the hg19 fasta file that will be used as the reference genome
ref="/path/vdas/test_exome/exome/hg19.fa"
#Tools and environment required to run this script
#Path to STAR,java,gatk and picard tools
java_home="/home/vdas/tools/jdk1.8.0_77/bin/java"
#location for bwa mem
bwamem="/path/softwares/bwa/bwa"
#location for the GATK toolkit jar file
gatk="/path/softwares/GATK-3.5/GenomeAnalysisTK.jar"
#location for picard for marking duplicates
picard="/home/vdas/tools/picard-tools-2.2.0/picard.jar"
#location for sambamba which is used for processing the bam files. can be found here http://lomereiter.github.io/sambamba/
sambamba="/home/vdas/tools/sambamba/build/sambamba"

#Path to gatk bundle set files ( these files are used for the different base quality calibration steps that are useful for the variant
#workflow. These bundles were downloaded from GATK 2.8 since they were same as for 3.5
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

# Aligning the paired end read mates to the reference genome hg19 with bwa mem. This is for paired-end data
echo -e "["$(date)"]\tAligning.."
$bwamem mem -t 8 -M -R "@RG\tID:${bn}\tLB:PairedEnd\tPL:Illumina\tPU:000000000-A442D\tSM:${bn}" $bwa_idx $fout/${bn}_R1.fastq $fout/${bn}_R2.fastq | samtools view -bS -F 0x04 - > $fout/${bn}"_Aligned.out.bam"

#sorting the bam files with sambamba
echo -e "["$(date)"]\tSorting.."
$sambamba sort -o $fout/${bn}"_sorted.bam" -p -t 12 $fout/${bn}"_Aligned.out.bam"

#rm $fout/${bn}"_Aligned.out.bam"
rm $fout/${bn}_R1.fastq
rm $fout/${bn}_R2.fastq

# indexing the sorted aligned bam file
echo -e "["$(date)"]\tIndexing.."
$sambamba index -p -t 12 $fout/${bn}"_sorted.bam"

#picard mark duplicates
#Remove the reads due to PCR duplicates or one can mark them with a flag
echo -e "["$(date)"]\tMarking duplicates.."
$java_home -Xmx2g -jar $picard MarkDuplicates I=$fout/${bn}"_sorted.bam" O=$fout/${bn}"_dupMarked.bam" M=$fout/${bn}"_dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 2>$fout/${bn}.MarkDuplicates.log
### Remove the aligned files which are not required downstream for memory management
rm $fout/${bn}"_Aligned.out.bam"
rm $fout/${bn}"_sorted.bam"
rm $fout/${bn}"_sorted.bam.bai"

#SplitNCigarReads
#echo -e "["$(date)"]\tSpliting reads.."
#$java_home -Xmx2g -jar $gatk -T SplitNCigarReads -R $ref -I $fout/${bn}"_dupMarked.bam" -o $fout/${bn}"_split.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2>$fout/${bn}.SplitNCigarReads.log

#Now the GATK 3.5 process handles starts from below
#Step 1 Local realignment around the INDELS (please refer to http://gatkforums.broadinstitute.org/gatk/discussion/38/local-realignment-around-indels#latest)
#Create targets for indel realignment
echo -e "["$(date)"]\tCreating targets for indel realignment.."
$java_home -Xmx2g -jar $gatk -T RealignerTargetCreator -R $ref -I $fout/${bn}"_dupMarked.bam" -o $fout/${bn}".intervals" -nt 12 -known $millsIndels -known $KGIndels -L $ss4exonbaits 2>$fout/${bn}.indel.log
#
#Perform indel realignment
echo -e "["$(date)"]\tPerforming Indel Realignment.."
$java_home -Xmx2g -jar $gatk -T IndelRealigner -R $ref -I $fout/${bn}"_dupMarked.bam" -targetIntervals $fout/${bn}".intervals" -known $millsIndels -known $KGIndels -o $fout/${bn}"_realigned.bam" 2>$fout/${bn}.indel2.log

#rm $fout/${bn}"_split.bam"
#rm $fout/${bn}"_split.bai"

rm $fout/${bn}"_dupMarked.bam"
rm $fout/${bn}"_dupMarked.bai"

#Perform BQSR (base quality recalibration after realigment step. This is a standard test that is performed to Detect systematic errors 
# in base quality scores. 
#Please refer to https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)
echo -e "["$(date)"]\tPerforming BQSR.."
$java_home -Xmx2g -jar $gatk -T BaseRecalibrator -I $fout/${bn}"_realigned.bam" -R $ref -knownSites $KGIndels -knownSites $millsIndels -knownSites $dbSNP138 -o $fout/${bn}"_recal.table" -L $ss4exonbaits 2>$fout/${bn}.BQSR.log

#Print recalibrated reads. This handle writes out sequence read data (for filtering, merging, subsetting etc)
#Please refer to https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php
echo -e "["$(date)"]\tPrinting recalibrated reads.."
$java_home -Xmx2g -jar $gatk -T PrintReads -R $ref -I $fout/${bn}"_realigned.bam" -nct 12 -BQSR $fout/${bn}"_recal.table" -o $fout/${bn}"_recal.bam" 2>$fout/${bn}.BQSR2.log

#Call the variants post realignment and recalibration with Haplotype caller. It calls germline SNPs and indels via local re-assembly of haplotypes
#Run HaplotypeCaller (please refer to https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)
echo -e "["$(date)"]\tRunning HaplotypeCaller.."
$java_home -Xmx2g -jar $gatk -T HaplotypeCaller -R $ref -I $fout/${bn}"_recal.bam" -L $ss4exonbaits -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $fout/${bn}".vcf" 2>$fout/${bn}.HaplotypeCaller.log
# Seperate out the SNPs from the from INDELs followed by filtration
#Filter variants for SNPs and then filter for qaulity scores. Variant filtration methods
echo -e "["$(date)"]\tFiltering Point Variants.."

$java_home -Xmx2g -jar $gatk -T SelectVariants -R $ref -V $fout/${bn}".vcf" -selectType SNP -o $fout/${bn}"_raw_snps.vcf"

#### filter over snps for scores that includes the combination of different Variant filter expression
# Filter variant calls based on INFO and FORMAT annotations
# Please refer to https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php

$java_home -Xmx2g -jar $gatk -T VariantFiltration -R $ref -V $fout/${bn}"_raw_snps.vcf" --filterExpression "DP < 5 || SB > -0.1 || QUAL < 50.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > 500" --filterName "DAS_snp_filters" -o $fout/${bn}"_filtered_snps.vcf" 2>$fout/${bn}.VariantFilter_snps.log

echo -e "["$(date)"]\tFiltering INDEL Variants.."
###select indels with the SelectVariants handle and -selectType process to give raw indels
$java_home -Xmx2g -jar $gatk -T SelectVariants -R $ref -V $fout/${bn}".vcf" -selectType INDEL -o $fout/${bn}"_raw_indels.vcf"

#### filter over raw indel indel files based on ased on INFO and FORMAT annotations
$java_home -Xmx2g -jar $gatk -T VariantFiltration -R $ref -V $fout/${bn}"_raw_indels.vcf" --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filterName "DAS_indel_filters" -o $fout/${bn}"_filtered_indels.vcf" 2>$fout/${bn}.VariantFilter_indels.log

echo -e "["$(date)"]\tDONE!"
