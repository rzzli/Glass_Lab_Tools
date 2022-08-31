#!/bin/bash
########################################
# This script takes a pair of atac tags (must be two atac tag dir for idr)
#     and together with 27ac tag dir (can be any number) to call enhancers.
# Input: atac tag x2; 
#        chip tag (any number)
#        project name (anything you can identify with, e.g. microglia_enhancer1
#        genome as homer accepts (e.g. mm10, hg38)
#        out dir (where to score output files, may not exist)
# Optional input: idr_p value, default 0.05
#                 chip_tag_min, default 30
#                 out peak size, default 300
# Output: bed file
#         fa file
#         npy one hot file (npeak, peaksize(300), 4(ACGT)), sorted based on peakId in bed file.
# usage:./tag_enhancer.sh "atac_tag*" "chip_tag*" my_project_name mm10, out_dir
#
# /home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer.sh "/home/zhl022/daima/projects/make_peaks_Aug30/il4/tags/atac/ctl/*"  "/home/zhl022/daima/projects/make_peaks_Aug30/il4/tags/chip/ctl/*" test1 mm10 ./myo
# or 
# /home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer.sh "/home/zhl022/daima/projects/make_peaks_Aug30/il4/tags/atac/ctl/C57Bl6_F_BMDM_ATAC_notx_rep1_MAH_l20180703 /home/zhl022/daima/projects/make_peaks_Aug30/il4/tags/atac/ctl/C57Bl6_F_BMDM_ATAC_notx_rep2_MAH_l20180703" "/home/zhl022/daima/projects/make_peaks_Aug30/il4/tags/chip/ctl/C57Bl6_F_BMDM_ChIP_H3K27ac_notx_MAH_l20171004 /home/zhl022/daima/projects/make_peaks_Aug30/il4/tags/chip/ctl/C57Bl6_F_BMDM_ChIP_H3K27ac_notx_rep2_MAH_l20171121" test1 mm10 /home/zhl022/daima/projects/make_peaks_Aug30/il4/tags/jpt/myo1
########################################
source /gpfs/data01/glasslab/home/zhl022/anaconda3/etc/profile.d/conda.sh
conda activate rickli
echo $CONDA_DEFAULT_ENV

atac_tags=( $1 )
chip_tags=( $2 )
project_name=$3
genome=$4
out_wk_dir=$5

idr_p=$6
chip_tag_min=$7
out_size=$8

[ -z "$idr_p" ] && idr_p=0.05
[ -z "$chip_tag_min" ] && chip_tag_min=30
[ -z "$out_size" ] && out_size=300

echo $out_wk_dir

arr_chip=("${chip_tags[@]}") #https://stackoverflow.com/questions/19458104/how-do-i-pass-a-wildcard-parameter-to-a-bash-file
chip_size=${#arr_chip[@]}  # number of chips entered, needed for selecting enhancers
mkdir -p $out_wk_dir
####### step 0: set constant variables
homer_pre=/bioinformatics/homer/bin/
idr_script=/home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer/run_idr_homerPeaks_rickli.py
python_path=/home/zhl022/.conda/envs/rickli/bin/python3
genome_pre=/home/zes017/homer/data/genomes/
idr_out=${out_wk_dir}/idr/
enhancer_python=/home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer/enhancer.py
one_hot_python=/home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer/faToOneHot.py

####### 1. atac peaks
for atac in "${atac_tags[@]}";do
    bname=$(basename $atac)
    ${homer_pre}/findPeaks $atac -style factor -o ${out_wk_dir}/${bname}_tag_atac_peaks.txt
    echo $bname
done

####### 2. run_idr
all_peaks=$(ls ${out_wk_dir}/*tag_atac_peaks.txt)
$python_path $idr_script -threshold $idr_p ${all_peaks} $idr_out

####### 3. resize idr and annotate with chip
#idr_result_path=$(ls ${idr_out}/*_idr.tsv)
scale_size_idr=1000
idr_scaled_outfile=${out_wk_dir}/idr_scaled_${scale_size_idr}.tsv
anno_outfile=${out_wk_dir}/${project_name}_anno.txt

awk -v scale_size_idr=$scale_size_idr '{OFS="\t"; NR>1; a=int(($3+$4)/2-scale_size_idr/2); $3=a; $4=a+scale_size_idr;print}' ${idr_out}/*_idr.tsv >${idr_scaled_outfile}
${homer_pre}/annotatePeaks.pl ${idr_scaled_outfile} $genome -d "${chip_tags[@]}" > $anno_outfile

####### 4. resize idr and annotate with chip

$python_path $enhancer_python $anno_outfile $chip_size $chip_tag_min $out_size

####### 5. peaks to bed format, and make fa, make one hot
for peaks in $(ls -d ${out_wk_dir}/*active_peak.txt);do
    bname=$(basename $peaks .txt)
    ${homer_pre}/pos2bed.pl $peaks | awk '{OFS="\t"; NR>1; $2=$2+1; print}' > ${out_wk_dir}/${bname}.bed  # need add 1 to start in bed file
    ${homer_pre}/homerTools extract ${out_wk_dir}/${bname}.bed ${genome_pre}/$genome -fa > ${out_wk_dir}/${bname}.fa
    $python_path $one_hot_python ${out_wk_dir}/${bname}.fa $out_size
done