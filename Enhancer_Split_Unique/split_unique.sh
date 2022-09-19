#!/bin/bash

source /gpfs/data01/glasslab/home/zhl022/anaconda3/etc/profile.d/conda.sh
conda activate rickli

sample1_dir=$1
sample1_name=$2

sample2_dir=$3
sample2_name=$4

output_dir=$5

genome=$6
merge_d=$7
return_full=$8
out_size=$9

[ -z "$merge_d" ] && merge_d=200
[ -z "$return_full" ] && return_full=True
[ -z "$out_size" ] && out_size=300

mkdir -p $output_dir

####### step 0: set constant variables
homer_pre=/bioinformatics/homer/bin/
python_path=/home/zhl022/.conda/envs/rickli/bin/python3
genome_pre=/home/zes017/homer/data/genomes/

merge_split_py=/home/zhl022/daima/to_share/Glass_Lab_Tools/Enhancer_Split_Unique/merge_output_to_unique_bed.py
one_hot_python=/home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer/faToOneHot.py

################################
split=("train" "val" "tv")
region=("enhancer" "promoter" "full")
for spl in ${split[@]};do
    for reg in ${region[@]};do
        ${homer_pre}/mergePeaks -d $merge_d ${sample1_dir}/${sample1_name}_${spl}_${reg}_active_peak.bed ${sample2_dir}/${sample2_name}_${spl}_${reg}_active_peak.bed > ${output_dir}/mergefile_${sample1_name}_${sample2_name}_${spl}_${reg}_active_peak.txt
        $python_path $merge_split_py ${sample1_dir}/${sample1_name}_${spl}_${reg}_active_peak.bed ${sample2_dir}/${sample2_name}_${spl}_${reg}_active_peak.bed ${output_dir}/mergefile_${sample1_name}_${sample2_name}_${spl}_${reg}_active_peak.txt ${output_dir} ${return_full}
    done
done

echo 'done with split unique'

for peaks in $(ls -d ${output_dir}/*bed);do
    bname=$(basename $peaks .bed)

    ${homer_pre}/homerTools extract ${output_dir}/${bname}.bed ${genome_pre}/$genome -fa > ${output_dir}/${bname}.fa
    $python_path $one_hot_python ${output_dir}/${bname}.fa $out_size
done