#!/bin/bash

##########
# This script convert input bed file to promoter, enhancer and all peaks fa file
#  1. rescale to 300bp bed 
#  2. annotatepeak using homer
#  3. split based on annotation file, promoter or enhancer
#  4. make fa files 
#
#  input: arg 1: input bed, arg 2= outdir, arg 3: genome arg 4: bed or fa or both
#  $ /home/zhl022/daima/pipeline/peak_to_enhancer_promoter_split.sh /bed/path /out/dir hg38 both
###########
input_bed=$1
outdir=$2
genome=$3
outtype=$4  #fa or bed or both

[ -z "$genome" ] && genome=hg38|| echo $genome
[ -z "$outtype" ] && outtype=both|| echo $outtype


mkdir -p $outdir

bname=$(basename $input_bed .bed)


## pre make files
rescaled_bed=${outdir}/${bname}_300bp.bed
anno_all_peak_300=${outdir}/${bname}_300bp_anno.txt
anno_enhancer_300=${outdir}/${bname}_300bp_anno_enhancer.txt
anno_promoter_300=${outdir}/${bname}_300bp_anno_promoter.txt
bed_enhancer=${outdir}/${bname}_enhancer.bed
bed_promoter=${outdir}/${bname}_promoter.bed
bed_all_peak=${outdir}/${bname}_all_peaks.bed
fa_enhancer=${outdir}/${bname}_enhancer.fa
fa_promoter=${outdir}/${bname}_promoter.fa
fa_all_peak=${outdir}/${bname}_all_peaks.fa

genome_pre=/home/zes017/homer/data/genomes/
homer_pre=/bioinformatics/homer/bin/

# step 1 rescale
awk '{OFS="\t"; a=int(($2+$3)/2-150); $2=a; $3=a+300;print}' ${input_bed} > $rescaled_bed

#step 2 annotation
${homer_pre}/annotatePeaks.pl $rescaled_bed $genome > $anno_all_peak_300

#step 3 select promoter and enhancer subfiles
awk ' ($8 ~ /intron/) || ($8 ~ /Intergenic/) { print } '  $anno_all_peak_300  > $anno_enhancer_300
awk ' ($8 ~ /promoter-TSS/)   { print } '   $anno_all_peak_300   > $anno_promoter_300

#step 4 extract fa / or bed or both
#homerTools extract $anno_enhancer_300 ${genome_pre}/$genome -fa > $fa_enhancer
#homerTools extract $anno_promoter_300 ${genome_pre}/$genome -fa > $fa_promoter
#homerTools extract $anno_all_peak_300 ${genome_pre}/$genome -fa > $fa_all_peak

if [ $outtype == "fa" ]
then
    ${homer_pre}homerTools extract $anno_enhancer_300 ${genome_pre}/$genome -fa > $fa_enhancer
    ${homer_pre}homerTools extract $anno_promoter_300 ${genome_pre}/$genome -fa > $fa_promoter
    ${homer_pre}homerTools extract $anno_all_peak_300 ${genome_pre}/$genome -fa > $fa_all_peak
elif [ $outtype == "bed" ]
then
    ${homer_pre}pos2bed.pl $anno_enhancer_300 > $bed_enhancer
    ${homer_pre}pos2bed.pl $anno_promoter_300 > $bed_promoter
    ${homer_pre}pos2bed.pl $anno_all_peak_300 > $bed_all_peak
    
else
    ${homer_pre}homerTools extract $anno_enhancer_300 ${genome_pre}/$genome -fa > $fa_enhancer
    ${homer_pre}homerTools extract $anno_promoter_300 ${genome_pre}/$genome -fa > $fa_promoter
    ${homer_pre}homerTools extract $anno_all_peak_300 ${genome_pre}/$genome -fa > $fa_all_peak
    ${homer_pre}pos2bed.pl $anno_enhancer_300 > $bed_enhancer
    ${homer_pre}pos2bed.pl $anno_promoter_300 > $bed_promoter
    ${homer_pre}pos2bed.pl $anno_all_peak_300 > $bed_all_peak
fi
## remove intermediate files # FOR NOW Mar 4, 2022

rm $rescaled_bed 
rm $anno_all_peak_300 
rm $anno_enhancer_300 
rm $anno_promoter_300 