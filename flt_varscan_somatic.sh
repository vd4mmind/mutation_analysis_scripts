#!/bin/bash

var=$1
bam=$2
bn=$3
ref="/path/group/UCSC_"type"/hg19/Sequence/hg19.fa"

echo $var
echo $bam
echo $bn

mkdir -p $bn"_fpFilter"

#cat $var | sed '/^#/d' | sed 1d | awk '{print $1"\t"$2"\t"$2}' > ./$bn_fpFilter/$bn.var

sed 1d $var | awk '{print $1"\t"$2"\t"$2}' > ./$bn"_fpFilter"/$bn.var #for varscan2 somatic

#awk '{print $1"\t"$2"\t"$2}' $var > ./$bn"_fpFilter"/$bn.var #for 5 column format

bam-readcount -q15 -w1 -b15 -l ./$bn"_fpFilter"/$bn.var -f $ref $bam > ./$bn"_fpFilter"/$bn.readCounts

cat $var | sed '/^#/d' | sed 1d | cut -f 1-4 > ./$bn"_fpFilter"/$bn.var #use this for varscan
#cat $var | cut -f 1,2,4,5 > ./$bn"_fpFilter"/$bn.var #for 5 column format

#cat $var | sed '/^#/d'| sed 1d | cut -f 1,2,4,5 > ./$bn_fpFilter/$bn.var #for mutect


################################ This is using fpfilter standalone script
perl /path/tool/varscan/fpfilter.pl --var-file ./$bn"_fpFilter"/$bn.var --readcount-file ./$bn"_fpFilter"/$bn.readCounts --min-depth 14 --min-read-pos 0.1 --min-strandedness 0.01 --max-mmqs-diff 90 --min-var-count 8 --min-var-frac 0.08 --min-var-dist-3 0.1 --max-var-mmqs 100 --output-file ./$bn"_fpFilter"/$bn.fpfilter

grep 'PASS' ./$bn"_fpFilter"/$bn.fpfilter | awk '{OFS="\t" ; print $1,$2,$2,$3,$4,$5,$6,$7}' > ./$bn"_fpFilter"/$bn.pass.annovar
################################


############################## This is using varscan inbuilt fpfilter
#java -jar /path/tool/varscan/VarScan.v2.4.0.jar fpfilter $var ./$bn"_fpFilter"/$bn.readCounts --output-file ./$bn"_fpFilter"/$bn.var --filtered-file ./$bn"_fpFilter"/$bn.failed.var --min-var-count 6 --min-var-freq 0.08 --min-var-readpos 0.1 --min-var-dist3 0.1 --min-strandedness 0.01 --min-strand-reads 5 --min-ref-basequal 30 --min-var-basequal 30 --max-rl-diff 25 --max-var-mmqs 100 --max-mmqs-diff 90 --min-ref-mapqual 30 --min-var-mapqual 30  --keep-failures 1

#awk '{OFS="\t" ; print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$40}' ./$bn"_fpFilter"/$bn.var | sed 's/%//g' | sed 1d | awk '{OFS="\t" ; print $1,$2,$2,$3,$4,$5+$6,$5,$6,$7/100,$8+$9,$8,$9,$10/100,$11}' > ./$bn"_fpFilter"/$bn.annovar

#grep 'PASS$' ./$bn"_fpFilter"/$bn.annovar > ./$bn"_fpFilter"/$bn.pass.annovar
#grep -v 'PASS$' ./$bn"_fpFilter"/$bn.annovar > ./$bn"_fpFilter"/$bn.failed.annovar
###############################

#perl /path/tools/annovar/table_annovar.pl ./$bn"_fpFilter"/$bn.pass.annovar /path/tools/annovar/humandb/ -buildver hg19 -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2014oct_all,snp135,snp135NonFlagged,snp131,clinvar_20150330,cosmic70,ljb2_all -operation g,r,r,f,f,f,f,f,f,f,f -nastring NA --otherinfo --remove --outfile ./$bn"_fpFilter"/$bn

#awk '$6=="exonic" {print $0}' ./$bn"_fpFilter"/$bn.hg19_multianno.txt |  awk '{OFS="\t" ; print $1,$2,$3,$4,$5,$7}' > ./$bn"_fpFilter"/$bn.exonic

#cat ./$bn"_fpFilter"/$bn.hg19_multianno.txt | awk '$6=="splicing" {print $0}' | awk '{OFS="\t" ; print $1,$2,$3,$4,$5,$7}' >> ./$bn"_fpFilter"/$bn.exonic

echo -e "\ndone"
