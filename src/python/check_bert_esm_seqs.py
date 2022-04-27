import pandas as pd
from Bio import SeqIO
import sys
import os


# check that the sequences in the protBERT and esm results are the same as in the fasta sequences, otherwise, print the gene names where the sequences don't match.


#------------------------------ main ------------------------------------------

input_path_bert =  "../../data/protBERT_data_4122/protbert (2).csv"
input_path_esm = "../../data/ESM1b_data_4122/esm (3).csv" #CNN results containing first chain of every protein only
input_path_fasta = "../../data/filtered_4122_seqs.fasta"



# esm predictions:
esm_df = pd.read_csv(input_path_esm)

# grouping df by pdb_id and chain_id: 
grouped = esm_df.groupby([esm_df.pdb_id, esm_df.chain_id])
esm_dict = dict(tuple(grouped))

esm_seqs = {}
for key, value in esm_dict.items():
    pdb_id = key[0] 
    chain_id = key[1]
    df = value
    seq_list = df['wtAA'].tolist()
    seq = ''.join(seq_list)
    name = pdb_id + "_" + chain_id
    esm_seqs[name] = seq


# now we need to parce through the fasta file and compare it with esm_seqs{}:
fasta_seqs = list(SeqIO.parse(input_path_fasta, "fasta"))
for rec in fasta_seqs:
    fasta_sequence = str(rec.seq)
    name = str(rec.name)

    if name in esm_seqs and esm_seqs[name] == fasta_sequence:
        continue
    elif name in esm_seqs and esm_seqs[name] != fasta_sequence:
        print("These esm seqs do not match the fasta seqs :")
        print(name)
    else:
        print("These gene was not found in the esm predictions:")
        print(name, "length:", len(fasta_sequence))


#--------------------------------------------------------------------------------------------

# bert predictions:
bert_df = pd.read_csv(input_path_bert)

# grouping df by pdb_id and chain_id: 
grouped = bert_df.groupby([bert_df.pdb_id, bert_df.chain_id])
bert_dict = dict(tuple(grouped))

bert_seqs = {}
for key, value in bert_dict.items():
    pdb_id = key[0] 
    chain_id = key[1]
    df = value
    seq_list = df['wtAA'].tolist()
    seq = ''.join(seq_list)
    name = pdb_id + "_" + chain_id
    bert_seqs[name] = seq


# now we need to parce through the fasta file and compare it with esm_seqs{}:
fasta_seqs = list(SeqIO.parse(input_path_fasta, "fasta"))
for rec in fasta_seqs:
    fasta_sequence = str(rec.seq)
    name = str(rec.name)

    if name in bert_seqs and bert_seqs[name] == fasta_sequence:
        continue
    elif name in bert_seqs and bert_seqs[name] != fasta_sequence:
        print("These bert seqs do not match the fasta seqs :")
        print(name)
    else:
        print("These genes were not found in the bert predictions:")
        print(name, len(fasta_sequence))

print("Done checking")