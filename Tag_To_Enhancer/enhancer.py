#!/home/zhl022/.conda/envs/rickli/bin/python3 
import sys
import pandas as pd

anno_path=str(sys.argv[1])
num_chip=int(sys.argv[2])  
chip_cutoff=int(sys.argv[3]) # default 30
scale_size=int(sys.argv[4]) # default 300
# set names 
full_active_path=anno_path.replace("_anno.txt","_tv_full_active_peak.txt")  # tv stand for both train and validation
enhancer_path=anno_path.replace("_anno.txt","_tv_enhancer_active_peak.txt")
promoter_active_path=anno_path.replace("_anno.txt","_tv_promoter_active_peak.txt")

full_active_path_train=anno_path.replace("_anno.txt","_train_full_active_peak.txt")  # train
enhancer_path_train=anno_path.replace("_anno.txt","_train_enhancer_active_peak.txt")
promoter_active_path_train=anno_path.replace("_anno.txt","_train_promoter_active_peak.txt")

full_active_path_val=anno_path.replace("_anno.txt","_val_full_active_peak.txt")  # val
enhancer_path_val=anno_path.replace("_anno.txt","_val_enhancer_active_peak.txt")
promoter_active_path_val=anno_path.replace("_anno.txt","_val_promoter_active_peak.txt")

#chip cutoff
df=pd.read_csv(anno_path,sep='\t')
df_sub=df[(df.iloc[:,-num_chip:].mean(axis=1))>chip_cutoff]

heldout=['chr1','chr8','chr18']

# resize
df_sub_start = df_sub["Start"]
df_sub_end = df_sub["End"]
mid = ((df_sub_start + df_sub_end)//2).astype(int)
df_sub["Start"] = mid - scale_size//2
df_sub["End"] = mid + scale_size//2

# select train vs val
df_sub_train = df_sub[~df_sub["Chr"].isin(heldout)]  # train df
df_sub_val = df_sub[df_sub["Chr"].isin(heldout)] # val df

# full distal and promoter
df_sub_distal = df_sub[(df_sub['Annotation'].str.contains('Intergenic') ) | (df_sub['Annotation'].str.contains('intron') )]
df_sub_promoter = df_sub[(df_sub['Annotation'].str.contains('promoter') )  ]

# train distal and promoter
df_sub_distal_train = df_sub_train[(df_sub_train['Annotation'].str.contains('Intergenic') ) | (df_sub_train['Annotation'].str.contains('intron') )]
df_sub_promoter_train = df_sub_train[(df_sub_train['Annotation'].str.contains('promoter') )  ]

# val distal and promoter
df_sub_distal_val = df_sub_val[(df_sub_val['Annotation'].str.contains('Intergenic') ) | (df_sub_val['Annotation'].str.contains('intron') )]
df_sub_promoter_val = df_sub_val[(df_sub_val['Annotation'].str.contains('promoter') )  ]

## save csv
#full
df_sub.to_csv(full_active_path,header=True,index=False,sep='\t')
df_sub_distal.to_csv(enhancer_path,header=True,index=False,sep='\t')
df_sub_promoter.to_csv(promoter_active_path,header=True,index=False,sep='\t')
#train
df_sub_train.to_csv(full_active_path_train,header=True,index=False,sep='\t')
df_sub_distal_train.to_csv(enhancer_path_train,header=True,index=False,sep='\t')
df_sub_promoter_train.to_csv(promoter_active_path_train,header=True,index=False,sep='\t')
#val
df_sub_val.to_csv(full_active_path_val,header=True,index=False,sep='\t')
df_sub_distal_val.to_csv(enhancer_path_val,header=True,index=False,sep='\t')
df_sub_promoter_val.to_csv(promoter_active_path_val,header=True,index=False,sep='\t')