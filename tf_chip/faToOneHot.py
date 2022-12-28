#!/home/zhl022/.conda/envs/rickli/bin/python3 

import os
import numpy as np
import sys
from Bio import SeqIO, Seq
import Bio
import random
random.seed(10)
def seq2onehot(seq):
    onehot_list = []
    nucleotides = ["A", "C" ,"G", "T"]
    for nuc in seq:
        if nuc == "N":
            onehot = [0 for _ in range(len(nucleotides))]
            onehot_list.append(onehot)
        else:
            onehot = [0 for _ in range(len(nucleotides))]
            onehot[nucleotides.index(nuc)] = 1
            onehot_list.append(onehot)
    return onehot_list
def readFasta(fasta_file, skip_duplicate=True, fmt='fasta'):
    ############
    # Read in sequences
    ############Bio.Alphabet.IUPAC module
    #alphabet =  Bio.Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA() # need to use this alphabet for motif score calculation
    alphabet=Bio.Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA()
    id_seq_dict = {} # {sequenceID: fastq sequence}
    duplicate_keys = []
    for seq_record in SeqIO.parse(fasta_file, fmt):  
        seq_record.seq.alphabet = alphabet
        if seq_record.id in id_seq_dict.keys():
            duplicate_keys.append(seq_record.id)
        else:
            id_seq_dict[seq_record.id] = seq_record.seq
    # delete duplicate keys
    if skip_duplicate:
        for dk in duplicate_keys:
            del id_seq_dict[dk]
        if len(duplicate_keys) > 0:
            print('Ignore duplicate keys in %s: %s' % (fasta_file, duplicate_keys))
    return id_seq_dict
def fa_to_one_hot(fa_path):
    fa_dict=readFasta(fa_path)
    fakeys=sorted(list(fa_dict.keys()))
    fa_onehot=np.array(list(map(seq2onehot,[fa_dict[key] for key in fakeys])))
    return fa_onehot


if __name__ == "__main__":
    fa_file=str(sys.argv[1]) # should be only one
    peak_size=int(sys.argv[2]) # 300?
    out_npy=fa_file.replace(".fa",'.npy')
    
    peak_one_hot_array=np.empty(shape=(0,peak_size,4))
    
    temp_one_hot=fa_to_one_hot(fa_file)
    peak_one_hot_array=np.vstack((peak_one_hot_array,temp_one_hot)) 

    with open(out_npy, 'wb') as f:
        np.save(f, peak_one_hot_array)