#!/home/zhl022/.conda/envs/rickli/bin/python3 
import sys
import pandas as pd
import os
#########################
# This script takes merged file from homer 
#    and two bed files merged, the output 
#    is the bed file that is unique to each sample 
#    and option to also output the full bed file
# arg 1: bed path for file 1
# arg 2: bed path for file 2
# arg 3: homer merge output
# arg 4: output dir
# arg 5: if output full bed to the output dir
#########################
def flatten(l):
    return [item for sublist in l for item in sublist]


if __name__ == "__main__":

    bedA=str(sys.argv[1])
    bedB=str(sys.argv[2])  
    merge_file=str(sys.argv[3]) # default 30
    outdir=str(sys.argv[4])
    output_full=bool(sys.argv[5])
    
    bnameA=os.path.basename(bedA)
    bnameB=os.path.basename(bedB)
    
    outnameA=bnameA.replace('active','active_unique')
    outnameB=bnameB.replace('active','active_unique')


    in1bed=pd.read_csv(bedA,sep='\t',header=None)
    in2bed=pd.read_csv(bedB,sep='\t',header=None)

    merge_df=pd.read_csv(merge_file,sep='\t' )
    merge_df=merge_df[(merge_df.iloc[:,-2].isnull()) | (merge_df.iloc[:,-1].isnull())]


    # get unique index for peak file A
    uniqA_idx=[idx.split(',') for idx in merge_df.iloc[:,-2].dropna()]
    uniqA_idx=flatten(uniqA_idx)
    
    # get unique index for peak file A
    uniqB_idx=[idx.split(',') for idx in merge_df.iloc[:,-1].dropna()]
    uniqB_idx=flatten(uniqB_idx)
    
    in1bed[in1bed.iloc[:,3].isin(uniqA_idx)].to_csv(outdir+'/'+  outnameA,index=False,header=False,sep='\t')
    in2bed[in2bed.iloc[:,3].isin(uniqB_idx)].to_csv(outdir+'/'+  outnameB,index=False,header=False,sep='\t')
    
    if output_full:
        in1bed.to_csv(outdir+'/'+bnameA,index=False,header=False,sep='\t')
        in2bed.to_csv(outdir+'/'+bnameB,index=False,header=False,sep='\t')