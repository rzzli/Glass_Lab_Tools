#!/bin/bash
####### step -1: take input variables
source /gpfs/data01/glasslab/home/zhl022/anaconda3/etc/profile.d/conda.sh
conda activate rickli
echo $CONDA_DEFAULT_ENV

atac_tags=( $1 )
project_name=$2
genome=$3
out_wk_dir=$4



chip_tag_min=$5
out_size=$6
idr_p=$7

[ -z "$idr_p" ] && idr_p=0.05 || echo $idr_p
[ -z "$out_size" ] && out_size=300 || echo $out_size
[ -z "$chip_tag_min" ] && chip_tag_min=30 || echo $chip_tag_min

echo $out_wk_dir
mkdir -p $out_wk_dir

####### step 0: set constant variables
homer_pre=/bioinformatics/homer/bin/
idr_script=/home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer/run_idr_homerPeaks_rickli.py
python_path=/home/zhl022/.conda/envs/rickli/bin/python3
genome_pre=/home/zes017/homer/data/genomes/
idr_out=${out_wk_dir}/idr/
#enhancer_python=/home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer/enhancer.py
#one_hot_python=/home/zhl022/daima/to_share/Glass_Lab_Tools/Tag_To_Enhancer/faToOneHot.py
enhancer_python=/home/zhl022/daima/to_share/Glass_Lab_Tools/tf_chip/distal_tf.py
one_hot_python=/home/zhl022/daima/to_share/Glass_Lab_Tools/tf_chip/faToOneHot.py

####### 1. atac peaks
count=1
for atac in "${atac_tags[@]}";do
    bname=${project_name}_${count}
    ${homer_pre}/findPeaks $atac -style factor -o ${out_wk_dir}/${bname}_tag_atac_peaks.txt
    (( count++ ))
    echo $bname
done


####### 2. run_idr
all_peaks=$(ls ${out_wk_dir}/*tag_atac_peaks.txt)
$python_path $idr_script -threshold $idr_p ${all_peaks} $idr_out

####### 3. resize idr and annotate with chip
#idr_result_path=$(ls ${idr_out}/*_idr.tsv)
scale_size_idr=300
idr_scaled_outfile=${out_wk_dir}/idr_scaled_${scale_size_idr}.tsv
anno_outfile=${out_wk_dir}/${project_name}_anno.txt

awk -v scale_size_idr=$scale_size_idr '{OFS="\t"; NR>1; a=int(($3+$4)/2-scale_size_idr/2); $3=a; $4=a+scale_size_idr;print}' ${idr_out}/*_idr.tsv >${idr_scaled_outfile}
${homer_pre}/annotatePeaks.pl ${idr_scaled_outfile} $genome -d "${atac_tags[@]}" > $anno_outfile

####### 4. resize idr and annotate with chip

$python_path $enhancer_python $anno_outfile $out_size $chip_tag_min


####### 5. peaks to bed format, and make fa, make one hot
for peaks in $(ls -d ${out_wk_dir}/*active_peak.txt);do
    bname=$(basename $peaks .txt)
    ${homer_pre}/pos2bed.pl $peaks | awk '{OFS="\t"; NR>1; $2=$2+1; print}' > ${out_wk_dir}/${bname}.bed  # need add 1 to start in bed file
    ${homer_pre}/bed2pos.pl ${out_wk_dir}/${bname}.bed > $peaks
    ${homer_pre}/homerTools extract $peaks ${genome_pre}/$genome -fa > ${out_wk_dir}/${bname}.fa
    $python_path $one_hot_python ${out_wk_dir}/${bname}.fa $out_size
done

