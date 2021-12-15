from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one
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

def predict(model, data):
  """ for making predictions on any dataset """

  predictions = model.predict(
        x = data,
        batch_size = None,
        verbose = 0)

  return predictions

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

def get_gene_labels(cnn_prob, bert_prob, esm_prob, SPLIT):
  gene_labels = []
  exceptions = dict()

  for gene in bert_probs:
    for pos_wt in bert_probs[gene]:
      try: # making sure all vectors have same gene, position and wt_aa:
        bert_vector = bert_probs[gene][pos_wt]
        esm_vector = esm_probs[gene][pos_wt]
        cnn_vector = cnn_probs[gene][pos_wt]
      except:
        if gene not in exceptions:
          exceptions[gene] = 1
        elif gene in exceptions:
          exceptions[gene] += 1
        continue

      # getting gene lablels:
      position = pos_wt[0]
      gene_labels.append([gene, position])
    
  # splitting data:
  training_labels, validation_labels, testing_labels = split_data(gene_labels, SPLIT)

  return training_labels, validation_labels, testing_labels

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

def get_predicted_freq_aa(row, model):
  
  aa_cnn_list = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
  aa_trans_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

  aa_pred = ""
  current_higest_prob = 0

  for i in range(0, len(row), 1):
    prob = float(row[i])
    if prob > current_higest_prob:
      current_higest_prob = prob
      if model == "trans":
        aa_pred = aa_trans_list[i]
      elif model == "cnn":
        aa_pred = aa_cnn_list[i]

  return current_higest_prob, aa_pred

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
      position = int(line[2])

      if gene not in trans_probs_all:
        trans_probs_all[gene] = dict()
        trans_probs_all[gene][(position, wt_aa)] = trans_probs
      elif (position, wt_aa) not in trans_probs_all[gene]:
        trans_probs_all[gene][(position, wt_aa)] = trans_probs

  return trans_probs_all

def get_combined_data(cnn_probs, bert_probs, esm_probs):

  combined_data = []
  answers = []
  exceptions = dict()

  for gene in bert_probs:
    for pos_wt in bert_probs[gene]:
      try: # making sure all vectors have same gene, position and wt_aa:
        bert_vector = bert_probs[gene][pos_wt]
        esm_vector = esm_probs[gene][pos_wt]
        cnn_vector = cnn_probs[gene][pos_wt]
      except:
        if gene not in exceptions:
          exceptions[gene] = 1
        elif gene in exceptions:
          exceptions[gene] += 1
        continue

      # combining synchronied vectors into one:
      combined_vector = cnn_vector + bert_vector + esm_vector
      combined_data.append(combined_vector)

      # creating the answer vector:
      wt_aa = pos_wt[1]
      answer_vector = get_answer(wt_aa)
      answers.append(answer_vector)

  #print(exceptions)
  return combined_data, answers

def get_nn_stats(cnn_probs, bert_probs, esm_probs, SPLIT):
  nn_stats = []
  exceptions = dict()
  
  for gene in bert_probs:
    for pos_wt in bert_probs[gene]:
      try: # making sure all vectors have same gene, position and wt_aa:
        bert_vector = bert_probs[gene][pos_wt]
        esm_vector = esm_probs[gene][pos_wt]
        cnn_vector = cnn_probs[gene][pos_wt]
      except:
        if gene not in exceptions:
          exceptions[gene] = 1
        elif gene in exceptions:
          exceptions[gene] += 1
        continue

      # getting nn data:
      cnn_win_freq, cnn_win_aa = get_predicted_freq_aa(cnn_vector, "cnn")
      bert_win_freq, bert_win_aa = get_predicted_freq_aa(bert_vector, "trans")
      esm_win_freq, esm_win_aa = get_predicted_freq_aa(esm_vector, "trans")
    
      nn_stats.append([float(cnn_win_freq), cnn_win_aa, float(bert_win_freq), bert_win_aa, float(esm_win_freq), esm_win_aa])

  stats_train, stats_val, stats_test = split_data(nn_stats, SPLIT)

  return stats_train, stats_val, stats_test

def decode_predictions(predictions):
  """ Finds the winning prediction for each site in the protein/sequence. """
  
  aa_dict = {'H':8, 'E':6, 'D':3,  'R':1, 'K':11, 'S':15, 'T':16, 'N':2, 'Q':5, 'A':0, \
             'V':19, 'L':10, 'I':9, 'M':12, 'F':13, 'Y':18, 'W':17, 'P':14, 'G':7, 'C':4}

  # answer order: 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'
  #                0   1   2  3   4    5   6   7   8   9   10  11  12  13  14  15  16 17  18   19

  winners = []

  for r, pred in enumerate(predictions):
    max_prob = 0
    max_prob_index = 0

    for i in range(0, 20):
      if pred[i] > max_prob:
        max_prob = pred[i]
        max_prob_index = i
    
    for key, value in aa_dict.items():
      if value == max_prob_index:
        winners.append([key, float(max_prob)])

  return winners

def split_data(data, SPLIT):
  # getting number of elements for each set:
  data_len = len(data)
  train_split = math.floor(data_len * SPLIT[0])
  val_split = math.floor(data_len * SPLIT[1])
  test_split = math.floor(data_len * SPLIT[2])

  # getting indexes:
  start_train, end_train = 0, train_split
  start_val, end_val = end_train, end_train + val_split
  start_test, end_test = end_val, end_val + test_split

  # splitting data:
  training = data[start_train:end_train + 1]
  validation = data[start_val:end_val + 1]
  testing = data[start_test:end_test + 1]

  return training, validation, testing

def make_csv(predictions, winners, center_list, output_path, gene_labels, nn_stats, model_id):

  csv_file = output_path + "predictions_" + "_model_" + model_id + ".csv"

  with open(csv_file, 'w', newline='') as file:
    
      writer = csv.writer(file)
      writer.writerow(['gene','position','wt_aa','pred_aa','pred_freq','cnn_win_freq','cnn_win_aa','bert_win_freq','bert_win_aa','esm_win_freq','esm_win_aa','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'])
     
      for i in range(0, len(winners)):
        row = []
        row.append(gene_labels[i][0]) # gene
        row.append(gene_labels[i][1]) # position
        row.append(center_list[i][0]) # wt aa
        row.append(winners[i][0]) # pred aa
        row.append(winners[i][1]) # pred freq
        for nn_stat in nn_stats[i]:
          row.append(nn_stat) # nn original winners here
        for pred in predictions[i]: # predictions for all 20 aa
            row.append(pred)
        writer.writerow(row)

#========================================================================================================
# Setting the variables, parameters and data paths/locations:
#========================================================================================================
### data paths/locations
input_path = "../../data/transfer_learning_net/validation/" # optional
output_path = "../../data/transfer_learning_net/"
model_path =  "../../data/transfer_learning_net/network_output/model_combo3_new_run_1.h5"
model_id = "combo3_run_1"

#========================================================================================================
# Generating predictions from a trained model
#========================================================================================================

### loading data (combinind cnn and transformer predictions)
cnn_pred = "../../output/PSICOV_CNN_output/"
bert_pred = "../../output/PSICOV_BERT_predictions.csv"
esm_pred = "../../output/PSICOV_ESM1b_predictions.csv"

print()
print("retrieving probabilities from predictions...", timestamp(), "\n")
cnn_probs = get_cnn_probs(cnn_pred)
bert_probs  = get_trans_probs(bert_pred)
esm_probs = get_trans_probs(esm_pred)


# combine probabilities into one vector for each position:
print("combining data...", timestamp(), "\n")
combined_data, answers = get_combined_data(cnn_probs, bert_probs, esm_probs)

np.save(output_path + "combined_data.npy", combined_data, allow_pickle = True)
np.save(output_path + "combined_data_answers.npy", answers, allow_pickle = True)

print("loading combined data and answers...", timestamp(), "\n")
combined_data = np.load(output_path + "combined_data.npy", allow_pickle = True)
answers = np.load(output_path + "combined_data_answers.npy", allow_pickle = True)

#--------- optional -----------------------------------------------------------------------
# spliting data (skip this step if you are not using test data or validation data)
SPLIT = [0.8, 0.2, 0] # needs to match the split that was used for spitting original data

print("starting to split data...", timestamp(), "\n")
training, validation, testing = split_data(combined_data, SPLIT)
train_answers, val_answers, test_answers = split_data(answers, SPLIT)
print("got all split data, getting labels...", timestamp(), "\n")

print("validation:", len(validation))
print("val_answers:", len(val_answers))

#------s----------------------------------------------------------------------------------

# get gene label:
train_labels, val_labels, test_labels = get_gene_labels(cnn_probs, bert_probs, esm_probs, SPLIT) # these are gene labels for each position.
print("Made all the combined vectors and answers.", timestamp(), "\n")
np.save(output_path + "val_labels.npy", val_labels, allow_pickle = True)

print("loading val_labels...", timestamp(), "\n")
val_labels = np.load(output_path + "val_labels.npy", allow_pickle = True)
print("val_labels", len(val_labels))

# get nn stats (freq and aa predicted by each model)
print("Getting freq and aa pridicted by each individual model...", timestamp(), "\n")
nn_stats_train, nn_stats_val, nn_stats_test = get_nn_stats(cnn_probs, bert_probs, esm_probs, SPLIT)
np.save(output_path + "nn_stats_val.npy", nn_stats_val, allow_pickle = True) 

print("loading nn_stats_val...", timestamp(), "\n")
nn_stats_val = np.load(output_path + "nn_stats_val.npy", allow_pickle = True)
print("nn_stats_val:", len(nn_stats_val))


### loading trained model
model = load_model(model_path)
print("Loaded model:", timestamp(), "\n")

### making predictions
predictions = predict(model, validation)
print("Predictions made:", timestamp(), "\n")

### decoding predictions (gets the pred aa and the freq pred- [aa_pred, freq_pred])
winners = decode_predictions(predictions)
print("Predictions decoded:", timestamp(), "\n")

### get centers (decode centers)
center_list = decode_predictions(val_answers)
print("centers/answers decoded:", timestamp(), "\n")

### saving data in a CSV file
print("Starting to append predictions to CSV", timestamp(), "\n")
make_csv(predictions, winners, center_list, output_path, val_labels, nn_stats_val, model_id)
print("Finished making CSV:", timestamp(), "\n")
print("Completed!:", timestamp(), "\n")

