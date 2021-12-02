from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one
from predictor import predict
import pandas as pd
import sys
import os
import numpy as np
import math
import csv

try:
  from keras.models import load_model
  from keras.optimizers import Adam
  from keras.callbacks import ModelCheckpoint
except ImportError:
  from tensorflow.keras.models import load_model
  from tensorflow.keras.optimizers import Adam
  from tensorflow.keras.callbacks import ModelCheckpoint

import time
from datetime import datetime
def timestamp():
  return str(datetime.now().time())

# this script loads a trained model and generates predictions. 

def get_cnn_probs(cnn_pred):

  fileList_cnn = os.listdir(cnn_pred)
  cnn_probs_all = dict() # genes are keys, values are lists of vectors (a list of vectors per genes).
  
  for file in fileList_cnn:
    with open(cnn_pred + file, 'r') as openedFile: # input is CNN output file
      lines = [line.rstrip('\n') for line in openedFile]

    #remove header 
    for i in range(len(lines)):
      lines[i] = lines[i].split(",")
    header = lines[0] 
    del lines[0]
    
    # structure of "line" in lines:
    #pos, aa_id, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS,prGLN,prGLU,prGLY,prHIS,prILE,prLEU,prLYS,prMET,prPHE,prPRO,prSER,prTHR,prTRP,prTYR,prVAL,prHydrophobic,prAromatic,prPolarUncharged,prCationic,prAnionic,prCharged,prSmall,prSulfur,prAcyl,prAlcohol
    # 0     1      2        3       4    5     6      7          8            9         10      11     12    13     14    15    16    17   18    19    20    21     22     23    24   25    26    27    28    29      30             31         32               33        34
    
    aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}

    gene = str(file[0:4]).lower()
    for line in lines:
      cnn_probs = line[10:30]
      position = int(line[0]) + 1
      wt_aa = aaCodes[line[5]]

      if gene not in cnn_probs_all:
        cnn_probs_all[gene] = dict()
        cnn_probs_all[gene][(position, wt_aa)] = cnn_probs
      elif (position, wt_aa) not in cnn_probs_all[gene]:
        cnn_probs_all[gene][(position, wt_aa)] = cnn_probs

  return cnn_probs_all

def get_gene_labels():
  gene_labels = []

  return gene_labels

def get_answer(aa_wt):

  # structure of "line" in lines:
  #pos, aa_id, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS,prGLN,prGLU,prGLY,prHIS,prILE,prLEU,prLYS,prMET,prPHE,prPRO,prSER,prTHR,prTRP,prTYR,prVAL,prHydrophobic,prAromatic,prPolarUncharged,prCationic,prAnionic,prCharged,prSmall,prSulfur,prAcyl,prAlcohol
  # 0     1      2        3       4    5     6      7          8            9         10      11     12    13     14    15    16    17   18    19    20    21     22     23    24   25    26    27    28    29      30             31         32               33        34
  
  aa_list_cnn =   ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
  aa_list = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
  
  answer = aa_wt
  answer_vector = [0] * 20
  answer_vector[aa_list_cnn.index(answer)] = 1
  
  return answer_vector

def get_trans_probs(trans_pred):

  trans_probs_all = {}

  with open(trans_pred, 'r') as openedFile: # input is CNN output file
    lines = [line.rstrip('\n') for line in openedFile]

    #remove header 
    for i in range(len(lines)):
      lines[i] = lines[i].split(",")
    header = lines[0] 
    del lines[0] 

    # trans data:
    # row,aa_wt,position,wt_prob, A, C, D, E, F, G, H,  I,  K,  L,  M,  N,  P,  Q,  R,  S,  T,  V,  W,  Y, gene
    # 0     1     2        3      4  5  6  7  8  9  10  11  12  13  14  15  16 17  18  19  20  21  22  23   24

    for line in lines:
      gene = line[24]
      trans_probs = line[4:24]
      wt_aa = line[1]
      position = line[2]

      if gene not in trans_probs_all:
        trans_probs_all[gene] = dict()
        trans_probs_all[gene][(position, wt_aa)] = trans_probs
      elif (position, wt_aa) not in trans_probs_all[gene]:
        trans_probs_all[gene][(position, wt_aa)] = trans_probs

  return trans_probs_all

def get_combined_data(cnn_probs, trans_probs):

  combined_data = []
  answers = []

  for gene_cnn in cnn_probs:
    for gene_trans in trans_probs:
      if gene_cnn == gene_trans:
        gene = gene_cnn
        cnn_data, trans_data = cnn_probs[gene], trans_probs[gene]

        for cnn_pos_wt, cnn_vector in zip(cnn_data.keys(), cnn_data.values()):
          cnn_pos = int(cnn_pos_wt[0])
          cnn_wt = cnn_pos_wt[1]
          for trans_pos_wt, trans_vector in zip(trans_data.keys(), trans_data.values()):
            trans_pos = int(trans_pos_wt[0])
            trans_wt = trans_pos_wt[1]

            if cnn_pos == trans_pos and cnn_wt == trans_wt:
              combined_vector = cnn_vector + trans_vector
              combined_data.append(combined_vector)

              answer_vector = get_answer(cnn_wt)
              answers.append(answer_vector)
     
  return combined_data, answers

def decode_predictions(predictions):
  """ Finds the winning prediction for each site in the protein/sequence. """
  
  aa_dict = {'H':8, 'E':6, 'D':3,  'R':1, 'K':11, 'S':15, 'T':16, 'N':2, 'Q':5, 'A':0, \
             'V':19, 'L':10, 'I':9, 'M':12, 'F':13, 'Y':18, 'W':17, 'P':14, 'G':7, 'C':4}

  # answer order: 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'
  #                0   1   2  3   4    5   6   7   8   9   10  11  12  13  14  15  16 17  18   19

  winners = []
  for pred in predictions:
      max_prob = 0
      max_prob_i = 0
      for i in range(0, 20):
          if pred[i] > max_prob:
              max_prob = pred[i]
              max_prob_i = i
      for key, value in aa_dict.items():
        if value == max_prob_i:
          winners.append(key)

  return winners

def make_csv(predictions, winners, center_list, model_id, output_path, gene_labels):
  csv_file = output_path + "predictions_" + "_model_" + model_id + ".csv"

  with open(csv_file, 'w', newline='') as file:
      writer = csv.writer(file)
      writer.writerow(['gene','position','wt_aa','pred_aa','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'])
      for i in range(0, len(winners)):
          row = []
          row.append(gene_labels[i])
          row.append(i)
          row.append(center_list[i])
          row.append(winners[i])
          for pred in predictions[i]:
              row.append(pred)
          writer.writerow(row)

#========================================================================================================
# Setting the variables, parameters and data paths/locations:
#========================================================================================================
### data paths/locations
input_path = "../../data/transfer_learning_net/validation/"
output_path = "../data/output/predictions/"
model_path =  "../../data/transfer_learning_net/network_output/model_6_run_1.h5"

### variables
model_id = "34"

#========================================================================================================
# Generating predictions from a trained model
#========================================================================================================

### loading data (combinind cnn and transformer predictions)
cnn_pred = "../../output/PSICOV_CNN_output/"
trans_pred = "../../output/PSICOV_BERT_predictions.csv"
cnn_probs = get_cnn_probs(cnn_pred)
trans_probs  = get_trans_probs(trans_pred)
combined_data, answers = get_combined_data(cnn_probs, trans_probs)
gene_labels = get_gene_labels(cnn_pred)
#get combined data with gene label.
print("Made all the combined vectors and answers:", timestamp(), "\n")

### loading trained model
model = load_model(model_path)
print("Loaded model:", timestamp(), "\n")

### making predictions
predictions = predict(model, combined_data)
print("Predictions made:", timestamp(), "\n")

### decoding predictions
winners = decode_predictions(predictions)
print("Predictions decoded:", timestamp(), "\n")

### saving data in a CSV file
print("Starting to append predictions to CSV", timestamp(), "\n")
make_csv(predictions, winners, output_path, gene_labels)
print("Finished making CSV:", timestamp(), "\n")
print("Completed!:", timestamp(), "\n")