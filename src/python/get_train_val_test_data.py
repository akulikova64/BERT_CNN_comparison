import os
import sys
import csv
import math
import numpy as np

# This program gets the training data for the transfer learining network.
 
def get_cnn_probs(cnn_pred):

  fileList_cnn = os.listdir(cnn_pred)
  cnn_probs_all = {} # genes are keys, values are lists of vectors (a list of vectors per genes).
  

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
    
    gene = str(file[0:4]).lower()
    for line in lines:
      cnn_probs = line[10:30]
      if gene not in cnn_probs_all:
        cnn_probs_all[gene] = cnn_probs

  return cnn_probs_all

def get_cnn_answers(cnn_pred):

  answers = []
  fileList_cnn = os.listdir(cnn_pred)
 
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
    
    aa_list = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

    for line in lines:
      answer = line[5]
      answer_vector = [0] * 20
      answer_vector[aa_list.index(answer)] = 1

    answers.append(answer_vector)
  
  return answers

def get_trans_probs(trans_pred):

  trans_probs_all = {}
  answers = []

  
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
    probs = line[4:24]
    trans_probs = np.array(probs)

    if gene not in trans_probs_all:
        trans_probs_all[gene] = trans_probs
    
  return trans_probs_all, answers

def get_trans_answers(trans_pred):

  answers = []
  fileList_cnn = os.listdir(trans_pred)
 
  for file in fileList_cnn:
    with open(trans_pred + file, 'r') as openedFile: # input is CNN output file
      lines = [line.rstrip('\n') for line in openedFile]

    #remove header 
    for i in range(len(lines)):
      lines[i] = lines[i].split(",")
    header = lines[0] 
    del lines[0]
    
    # trans data:
  # row,aa_wt,position,wt_prob, A, C, D, E, F, G, H,  I,  K,  L,  M,  N,  P,  Q,  R,  S,  T,  V,  W,  Y, gene
  # 0     1     2        3      4  5  6  7  8  9  10  11  12  13  14  15  16 17  18  19  20  21  22  23   24

    aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    for line in lines:
      answer = line[1]
      answer_vector = [0] * 20
      answer_vector[aa_list.index(answer)] = 1

    answers.append(answer_vector)

  return answers

def answers_matching(answers_cnn, answers_trans):

  for a_cnn, a_trans in zip(answers_cnn, answers_trans):
    if a_cnn != a_trans:
      return False
    else:
      return True


def get_combined_data(cnn_probs, trans_probs):

  combined_data = []

  for gene_cnn, gene_trans, cnn_values, trans_values in zip(cnn_probs.keys(), trans_probs.keys(), cnn_probs.values(), trans_probs.values()):
    for pos_cnn, pos_trans in enumerate(zip(cnn_values, trans_values)):
      if gene_cnn == gene_trans:

        cnn_values = [float(item) for item in cnn_values]
        trans_values = [float(item) for item in trans_values]
        combined_vector = cnn_values + trans_values
        combined_data.append(combined_vector)


  return combined_data

def get_combined_answers(answers_cnn, answers_trans):
  combined_answers = []

  return combined_answers

def split_data(probs, SPLIT):
  # getting number of elements for each set:
  data_len = len(probs)
  train_split = math.floor(data_len * SPLIT[1])
  val_split = math.floor(data_len * SPLIT[2])
  test_split = math.floor(data_len * SPLIT[3])

  # getting indexes:
  start_train, end_train = 0, train_split
  start_val, end_val = end_train, end_train + val_split
  start_test, end_test = end_val,  end_val + test_split

  # splitting data:
  training = probs[start_train:end_train + 1]
  validation = probs[start_val:end_val + 1]
  testing = probs[start_test:end_test + 1]

  return training, validation, testing



#---------------main-----------------------------

# data paths
cnn_pred = "../../output/PSICOV_CNN_output/"
trans_pred = "../../output/PSICOV_BERT_predictions.csv"

train_path = "../../data/transfer_learning_net/input/training/"
val_path = "../../data/transfer_learning_net/input/validation/"
test_path = "../../data/transfer_learning_net/input/testing/"

# split is the percentages of training, validation and testing data (from all data)
SPLIT = [0.6, 0.2, 0.2]

# extract the probabilities and wt (answers) from raw files for both models:
cnn_probs, answers_cnn = get_cnn_probs(cnn_pred)
trans_probs, answers_trans  = get_trans_probs(trans_pred)

# check that answers match for both models (will ensure that the probs also match by position and gene)
assert answers_matching(answers_cnn, answers_trans)

# combine probabilities into one vector for each position:
combined_data = get_combined_data(cnn_probs, trans_probs)
combined_answers = get_combined_answers(answers_cnn, answers_trans)

# split data into training, validation and testing
training, validation, testing = split_data(combined_data, SPLIT)
train_answers, val_answers, test_answers = split_data(combined_answers, SPLIT)


# saving data as numpy arrays:
np.save(train_path + "nn_data_train.npy", training, allow_pickle = True)
np.save(train_path + "answers_train.npy", train_answers, allow_pickle = True)
np.save(val_path + "nn_data_val.npy", validation, allow_pickle = True)
np.save(val_path + "answers_val.npy", val_answers, allow_pickle = True)

