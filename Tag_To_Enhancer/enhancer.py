#!/home/zhl022/.conda/envs/rickli/bin/python3 
import sys
import pandas as pd

anno_path=str(sys.argv[1])
num_chip=int(sys.argv[2])  
chip_cutoff=int(sys.argv[3]) # default 30
scale_size=int(sys.argv[4]) # default 300

full_active_path=anno_path.replace("_anno.txt","_full_active_peak.txt")
enhancer_path=anno_path.replace("_anno.txt","_enhancer_active_peak.txt")
promoter_active_path=anno_path.replace("_anno.txt","_promoter_active_peak.txt")

df=pd.read_csv(anno_path,sep='\t')
df_sub=df[(df.iloc[:,-num_chip:].mean(axis=1))>chip_cutoff]

df_sub_start = df_sub["Start"]
df_sub_end = df_sub["End"]
mid = ((df_sub_start + df_sub_end)//2).astype(int)
df_sub["Start"] = mid - scale_size//2
df_sub["End"] = df_sub["Start"] +300
df_sub_distal = df_sub[(df_sub['Annotation'].str.contains('Intergenic') ) | (df_sub['Annotation'].str.contains('intron') )]
df_sub_promoter = df_sub[(df_sub['Annotation'].str.contains('promoter') )  ]

df_sub.to_csv(full_active_path,header=True,index=False,sep='\t')
df_sub_distal.to_csv(enhancer_path,header=True,index=False,sep='\t')
df_sub_promoter.to_csv(promoter_active_path,header=True,index=False,sep='\t')