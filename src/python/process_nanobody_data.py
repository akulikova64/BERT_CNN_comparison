from Bio import SeqIO
import sys
import os
import numpy as np
import math
import csv

# this script makes a CSV of all nanobody data (processed)

def get_predicted_freq_aa(row):
  
  aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

  aa_pred = ""
  current_higest_prob = 0

  for i in range(0, len(row), 1):
    prob = float(row[i])
    if prob > current_higest_prob:
      current_higest_prob = prob
      aa_pred = aa_list[i]

  return current_higest_prob, aa_pred

def reorder_cnn_data(reference, nn_data):
  ref_pos_with_gap = []

  for i, aa in enumerate(reference):
    if aa == "-":
      ref_pos_with_gap.append(i)

  for pos in ref_pos_with_gap:
    nn_data.insert(pos, ["NA", "NA", "NA", "NA", "NA"])

  return nn_data

#-------------------------main-----------------------------------------------------------------
# data paths
cnn_pred = "../../output/cnn_3ogo.csv"
bert_pred = "../../output/bert_prediction_3ogo.csv"
#esm1b_pred = "../../output/esm1b_3ogo.csv"
nanobody_data = "../../data/NanobodyData2.csv"
nanobody_alignments = "../../data/nano_align.fasta"
output_file = "../../output/nanobody_summary_data.csv"

# write new CSV file:
with open(output_file, 'w', newline='') as outputfile:
  writer = csv.writer(outputfile)
  writer.writerow(['nanobody', 'position', 'wt_aa', 'ref_aa', 'change', 'wt_3ogo', 'cnn_aa', 'cnn_freq', 'bert_aa', 'bert_freq', 'mol_weight', 'yeild', 'molar_conc'])

  #-----------------------------------------------------------------------------------------------------
  # saving neural network data to a list "nn_data"
  with open(cnn_pred, newline = '') as cnn_file:
    cnn_reader = csv.reader(cnn_file, delimiter = ',')
    next(cnn_reader)
    with open(bert_pred, newline = '') as bert_file:
      bert_reader = csv.reader(bert_file, delimiter = ',')
      next(bert_reader)

      nn_data = []
    
      # cnn row:
      # "","gene","chain","position","aa_predicted","aa_wt","freq_predicted","freq_wt"
      # 0    1      2         3            4           5          6              7
      # bert/esm row:
      # "" , wt, wtIndex, wtScore, A,  C, D, E, F, G, H,  I,  K,  L,  M,  N,  P,  Q,  R,  S,  T,  V,  W,  Y
      #  0   1     2         3     4   5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21  22  23

      for cnn_row, bert_row in zip(cnn_reader, bert_reader):
        cnn_freq, cnn_aa = cnn_row[6], cnn_row[4]
        bert_freq, bert_aa = get_predicted_freq_aa(bert_row[4:23])
        #esm1b_freq, esm1b_aa =  get_predicted_freq_aa(esm_row[4:23])

        #make sure wt 3ogo amino acid is same
        assert bert_row[1] == cnn_row[5]
        wt_3ogo = cnn_row[5]
        
        # make sure positions are aligned:
        #assert cnn_row[0] == bert_row[2] and cnn_row[0] == esm_row[2]
        assert cnn_row[0] == bert_row[2] 

        nn_data.append([wt_3ogo, cnn_aa, cnn_freq, bert_aa, bert_freq])

  #--------------------------------------------------------------------------------------------
  # saving nanobody alignment as a record object:
  records = list(SeqIO.parse(nanobody_alignments, "fasta")) 
  reference = records[0].seq
  print("reference:", records[0].description)

  # reorder CNN data (add gaps to align to reference sequence) save for later
  nn_data_reordered = reorder_cnn_data(reference, nn_data)

  len_records = len(records)
  for i in range(1, len_records): # parsing though all the nanobodies in the alignment
    nanobody = records[i].description
    data_seq = str(records[i].seq) 

    # getting nanobody experimental data:
    with open(nanobody_data, newline = '') as nanofile:
      nano_reader = csv.reader(nanofile, delimiter = ',')

      # nano row:
      # aa Sequence,Nanobody,Date of purification, MW, Yeild, Molar_conc
      #    0            1         2               3     4       5    
      for nano_row in nano_reader:  # find matching nanobodies in the two datasets
        if nano_row[1].strip() == nanobody.strip():
          mol_weight = nano_row[3]
          yeild = nano_row[4]
          molar_conc = nano_row[5]
          break
        else:
          continue
      
    position = 0
    for i, aa in enumerate(reference):
      wt_aa = data_seq[i]
      ref_aa = aa

      if aa == data_seq[i] and aa != "-":
        change = "no_change"
      elif data_seq[i] == "-" and aa == "-":
        change = "NA"
      elif aa != data_seq[i] and aa != "-" and data_seq[i] != "-":
        change = "mutation"
      elif aa != data_seq[i] and aa == "-" and data_seq[i] != "-":
        # gap in the reference is insertion
        change = "insertion"
      elif aa != data_seq[i] and aa != "-" and data_seq[i] == "-":
        # gap in the alignment seq is a deletion
        change = "deletion"

      for i, nn_list in enumerate(nn_data_reordered):
        if i == position:
          writer.writerow([nanobody, position, wt_aa, ref_aa, change] + nn_list + [mol_weight, yeild, molar_conc])

      position += 1
  
      '''
      with open(esm1b_pred, newline = '') as esm_file:
        esm_reader = csv.reader(esm_file, delimiter = ',')'''