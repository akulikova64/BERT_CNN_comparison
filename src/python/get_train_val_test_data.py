import os
from re import T
import sys
import csv
import math
import numpy as np

# This program gets the training data for the transfer learining network.
 
def get_cnn_probs(cnn_pred):

  fileList_cnn = os.listdir(cnn_pred)
  cnn_probs_all = dict() # genes are keys, values are lists of vectors (a list of vectors per genes).
  
  for file in fileList_cnn:
    with open(cnn_pred + file, 'r') as openedFile: # input is CNN output file
      lines = [line.rstrip('\n') for line in openedFile]

    #remove header 
    for i in range(len(lines)):
      lines[i] = lines[i].split(",") #split elements of each row with ","
    header = lines[0] 
    del lines[0]
    
    # structure of "line" in lines (old):
    #pos, aa_id, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS,prGLN,prGLU,prGLY,prHIS,prILE,prLEU,prLYS,prMET,prPHE,prPRO,prSER,prTHR,prTRP,prTYR,prVAL,prHydrophobic,prAromatic,prPolarUncharged,prCationic,prAnionic,prCharged,prSmall,prSulfur,prAcyl,prAlcohol
    # 0     1      2        3       4    5     6      7          8            9         10      11     12    13     14    15    16    17   18    19    20    21     22     23    24   25    26    27    28    29      30             31         32               33        34
    

    # structure of "line" in lines:
    # , Unnamed: 0, model, accuracy, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS, prGLN, prGLU, prGLY, prHIS, prILE, prLEU, prLYS, prMET, prPHE, prPRO, prSER, prTHR, prTRP, prTYR, prVAL, position
    # 0,  1,         2        3        4       5       6     7     8      9         10            11         12    13      14    15     16     17      18     19     20    21    22     23       24     25     26    17      18     29     30     31     32      


    aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}

    gene = str(file[0:4]).lower()
    for line in lines:
      cnn_probs = line[12:31]
      position = int(line[32])
      wt_aa = aaCodes[line[7]]

      if gene not in cnn_probs_all:
        cnn_probs_all[gene] = dict()
        cnn_probs_all[gene][(position, wt_aa)] = cnn_probs
      elif (position, wt_aa) not in cnn_probs_all[gene]:
        cnn_probs_all[gene][(position, wt_aa)] = cnn_probs

  return cnn_probs_all

def get_answer(aa_wt):

  # getting answer for CNN data (answer/wt is exactly the same as trans data):
  # structure of "line" in lines (old):
  #pos, aa_id, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS,prGLN,prGLU,prGLY,prHIS,prILE,prLEU,prLYS,prMET,prPHE,prPRO,prSER,prTHR,prTRP,prTYR,prVAL,prHydrophobic,prAromatic,prPolarUncharged,prCationic,prAnionic,prCharged,prSmall,prSulfur,prAcyl,prAlcohol
  # 0     1      2        3       4    5     6      7          8            9         10      11     12    13     14    15    16    17   18    19    20    21     22     23    24   25    26    27    28    29      30             31         32               33        34
  
  # structure of "line" in lines:
  # , Unnamed: 0, model, accuracy, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS, prGLN, prGLU, prGLY, prHIS, prILE, prLEU, prLYS, prMET, prPHE, prPRO, prSER, prTHR, prTRP, prTYR, prVAL, position
  # 0,  1,         2        3        4       5       6     7     8      9         10            11         12    13      14    15     16     17      18     19     20    21    22     23       24     25     26    17      18     29     30     31     32      


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

    # trans data old (from my code):
    # row,aa_wt,position,wt_prob, A, C, D, E, F, G, H,  I,  K,  L,  M,  N,  P,  Q,  R,  S,  T,  V,  W,  Y, gene
    # 0     1     2        3      4  5  6  7  8  9  10  11  12  13  14  15  16 17  18  19  20  21  22  23   24

    # trans data from Tianlong: (do not use `pos` as the position, use the first col)
    # , model, accuracy, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prA, prR, prN, prD, prC, prQ, prE, prG, prH, prI, prL, prK, prM, prF, prP, prS, prT, prW, prY, prV
    # 0   1      2         3        4       5    6     7      8          9            10        11  12   13   14   15    16  17   18    19   20  21   22   23   24   25   26   27   28   29   30

    for line in lines:
      gene = line[3].lower()
      trans_probs = line[11:30]
      wt_aa = line[6]
      position = int(line[0]) + 1

      if gene not in trans_probs_all:
        trans_probs_all[gene] = dict()
        trans_probs_all[gene][(position, wt_aa)] = trans_probs
      elif (position, wt_aa) not in trans_probs_all[gene]:
        trans_probs_all[gene][(position, wt_aa)] = trans_probs

  return trans_probs_all



def get_combined_data(cnn_probs, bert_probs, esm_probs):

  combined_data = []
  answers = []

  for gene_cnn in cnn_probs:
    for gene_bert in bert_probs:
      for gene_esm in esm_probs:
        if gene_cnn == gene_bert and gene_cnn == gene_esm:
          gene = gene_cnn
          cnn_data, bert_data, esm_data = cnn_probs[gene], bert_probs[gene], esm_probs[gene]

          for cnn_pos_wt, cnn_vector in zip(cnn_data.keys(), cnn_data.values()):
            cnn_pos = int(cnn_pos_wt[0])
            cnn_wt = cnn_pos_wt[1]
            for bert_pos_wt, bert_vector in zip(bert_data.keys(), bert_data.values()):
              bert_pos = int(bert_pos_wt[0])
              bert_wt = bert_pos_wt[1]
              for esm_pos_wt, esm_vector in zip(esm_data.keys(), esm_data.values()):
                esm_pos = int(esm_pos_wt[0])
                esm_wt = esm_pos_wt[1]

                if cnn_pos == bert_pos and  cnn_pos == esm_pos and \
                  cnn_wt == bert_wt and cnn_wt == esm_wt: # checking that the position and wt match.
                  combined_vector = cnn_vector + bert_vector + esm_vector
                  combined_data.append(combined_vector)

                  answer_vector = get_answer(cnn_wt)
                  answers.append(answer_vector)
     
  return combined_data, answers

def get_combined_data_2(cnn_probs, bert_probs, esm_probs):

  combined_data = []
  answers = []
  exceptions = dict()
  genes = [] # genes that are common across all three nets

  '''
  for cnn_key, bert_key, esm_key in zip(cnn_probs['1pko'], bert_probs['1pko'], esm_probs['1pko']):
    print(cnn_key, bert_key, esm_key)
  sys.exit()'''

  for gene in bert_probs:
    for pos_wt in bert_probs[gene]:
      try: # making sure all vectors have same gene, position and wt_aa:
        bert_vector = bert_probs[gene][pos_wt]
        esm_vector = esm_probs[gene][pos_wt]
        cnn_vector = cnn_probs[gene][pos_wt]
        if gene not in genes:
          genes.append(gene)
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

  #print("All three networks did not have the same gene/position or both", exceptions)
  print("number of common genes across all three nets:", len(genes))
  return combined_data, answers

def split_data(probs, SPLIT):
  # getting number of elements for each set:
  data_len = len(probs)
  train_split = int(math.floor(data_len * SPLIT[0]))
  val_split = int(math.floor(data_len * SPLIT[1]))
  test_split = int(math.floor(data_len * SPLIT[2]))

  # getting indexes:
  start_train, end_train = 0, train_split
  start_val, end_val = end_train, end_train + val_split
  start_test, end_test = end_val, end_val + test_split

  # splitting data:
  training = probs[start_train:end_train + 1]
  validation = probs[start_val:end_val + 1]
  testing = probs[start_test:end_test + 1]

  return training, validation, testing


#--------------- main -----------------------------

# data paths
cnn_pred = "../../output/aligned_CNN_4122/"
bert_pred = "../../data/protBERT_data_4122/protbert (2).csv"
esm_pred = "../../data/ESM1b_data_4122/esm (3).csv"

train_path = "../../data/transfer_learning_net/training/"
val_path = "../../data/transfer_learning_net/validation/"
test_path = "../../data/transfer_learning_net/testing/"

# split is the percentages of training, validation and testing data (from all data)
#SPLIT = [0.6, 0.2, 0.2]
SPLIT = [0.95, 0.05, 0.0]

# extract the probabilities and wt (answers) from raw files for both models:
cnn_probs = get_cnn_probs(cnn_pred)
bert_probs  = get_trans_probs(bert_pred)
esm_probs  = get_trans_probs(esm_pred)

# combine probabilities into one vector for each position:
combined_data, answers = get_combined_data_2(cnn_probs, bert_probs, esm_probs)

# split data into training, validation and testing
training, validation, testing = split_data(combined_data, SPLIT)
train_answers, val_answers, test_answers = split_data(answers, SPLIT)

assert len(training) == len(train_answers)
assert len(validation) == len(val_answers)
assert len(testing) == len(test_answers)

# saving data as numpy arrays:
# training
np.save(train_path + "nn_data_train_4122.npy", training)
np.save(train_path + "answers_train_4122.npy", train_answers)

# validation
np.save(val_path + "nn_data_val_4122.npy", validation)
np.save(val_path + "answers_val_4122.npy", val_answers)

# testing
np.save(test_path + "nn_data_test_4122.npy", testing)
np.save(test_path + "answers_test_4122.npy", test_answers)

print()
print("Finished saving data!")














# --------------------------- EXTRA code (ignore) --------------------------------
'''
def get_trans_answers(trans_pred):
  
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

    aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    for line in lines:
      answer = line[1]
      answer_vector = [0] * 20
      answer_vector[aa_list.index(answer)] = 1
      answers.append(answer_vector)

  return answers
'''
'''
def answers_matching(answers_cnn, answers_trans):
  #                 0   1   2   3   4   5   6   7   8   9  10  11   12  13  14  15  16 17  18  19
  aa_list_cnn =   ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
  aa_list_trans = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

  # conversion from cnn aa encoding (key) to trans. aa encoding (value)
  conversion = {0:0, 1:14, 2:11, 3:2, 4:1, 5:13, 6:3, 7:5, 8:6, 9:7, 10:9, 11:8, 12:10, 13:4, 14:12, 15:15, 16:16, 17:18, 18:19, 19:17}

  for a_cnn, a_trans in zip(answers_cnn, answers_trans):
    cnn_index = a_cnn.index(1)
    trans_index = a_trans.index(1)

    if conversion[cnn_index] != trans_index:
      return False
    else:
      return True
'''