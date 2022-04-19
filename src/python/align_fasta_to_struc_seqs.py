import pandas as pd
import gemmi
from typing import Tuple
from typing import List
from typing import Dict
from numpy import isnan
from Bio import SeqIO
import sys
import os


# function that aligns the fasta sequence to the structure sequence:

def fetch_seq_from_nn_pred(nn_pred: pd.DataFrame, chain_id: str) -> Tuple[List[int], str]:

    chain_rows = (
        nn_pred[nn_pred["chain_id"] == chain_id]
        .sort_values(by="pos", ascending=True)
        .drop_duplicates(["pos", "wtAA"])
    )
    seq_nums, seq = list(zip(*chain_rows[["pos", "wtAA"]].values))
    return seq_nums, gemmi.one_letter_code(seq)

def map_fp_seqnums_to_pdb_seqnums(fp_seq: str, nn_pred: pd.DataFrame, chain_id: str) -> Dict[int, int]:
    """
    Returns a dict where:
        the key is the residue position from the fasta sequence
        the value is the residue position from the nn prediction csv (pdb/cif structure)
    """

    pdb_seq_ids, pdb_seq = fetch_seq_from_nn_pred(nn_pred, chain_id)
    alignment = gemmi.align_string_sequences(list(fp_seq), list(pdb_seq), [])

    #printing alignment data:
    #print(alignment.formatted(fp_seq, pdb_seq), end='')  
    #print(alignment.score)

    fp_gaps = alignment.add_gaps(fp_seq, 1)
    pdb_gaps = alignment.add_gaps(pdb_seq, 2)
    pdb_seq_ids_gen = (num for num in pdb_seq_ids)
    fp_pdb_mapping = dict()
    fp_idx = 0

    for idx in range(len(fp_gaps)):
        fp_char = fp_gaps[idx]
        pdb_char = pdb_gaps[idx]
        pdb_seq_id = None
        if pdb_char.isalpha():
            pdb_seq_id = next(pdb_seq_ids_gen)
        if fp_char.isalpha():
            fp_idx += 1
        if (
            pdb_char.isalpha()
            and fp_char.isalpha()
            and pdb_char.upper() == fp_char.upper()
            and not isnan(pdb_seq_id)
        ):
            fp_pdb_mapping[fp_idx] = int(pdb_seq_id)
            
    return fp_pdb_mapping

#---------------------------- main -----------------------------------------------

input_path_1 = "../../output/PSICOV_CNN_output_new/"
input_path_2 = "../../data/PSICOV_seqs.fasta"
output_path = "../../output/aligned_CNN_4122"

nn_df_List = os.listdir(input_path_1)
fasta_seqs = list(SeqIO.parse(input_path_2, "fasta"))

for nn_df_file in nn_df_List:

    gene_pdb = nn_df_file[0:4]
    nn_df = pd.read_csv("../../output/PSICOV_CNN_output_new/" + nn_df_file)

    for rec in fasta_seqs:
        fasta_sequence = str(rec.seq)
        gene_fasta = str(rec.name)

        # check if we are alignming the same protein:
        #if gene_pdb == "1bsg" and gene_pdb == gene_fasta:
        if gene_pdb == gene_fasta:
            #print(gene_pdb)

            #print(len(nn_df))
            #print(len(fasta_sequence))
        
            # the key is the residue position from the fasta sequence
            # the value is the residue position from the nn prediction csv (pdb/cif structure)
            alignment_dict = map_fp_seqnums_to_pdb_seqnums(fasta_sequence, nn_df, "A")
            #print(alignment_dict)

            '''
            old = 0
            for number in list(alignment_dict.values()):
                if number == old:
                    print(number)
                else:
                    old = number
            '''

            values_list = list(alignment_dict.values())
            keys_list = list(alignment_dict.keys())
            #print(len(values_list))
            #print(len(keys_list))

            try:  
                # selecting rows based on condition (keeping only chain A and keeping only rows that align with the sequences)
                nn_df_chainA = nn_df.loc[nn_df['chain_id'] == 'A']
                new_nn_df = nn_df_chainA.loc[nn_df_chainA['pos'].isin(values_list)]
                #print(new_nn_df)
                #print(len(new_nn_df)) 

                # adding a column to pandas dataframe:
                new_nn_df['position'] = keys_list

                new_nn_df.to_csv('../../output/aligned_CNN_PSICOV/' + gene_pdb + '.csv')
            except:
                print("An exception occurred on protein", gene_pdb)


print("Finished aligning all proteins!")


