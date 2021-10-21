import os
import sys
import csv
import math
import numpy as np

# get ensamble predictions for the cnn and transformer (just average the probabilities)

# trans data:
# row,aa_wt,position,wt_prob, A, C, D, E, F, G, H,  I,  K,  L,  M,  N,  P,  Q,  R,  S,  T,  V,  W,  Y, gene
# 0     1     2        3      4  5  6  7  8  9  10  11  12  13  14  15  16 17  18  19  20  21  22  23   24

# cnn data:
# position,gene,wt_aa,q_A,q_R,q_N, q_D, q_C, q_Q, q_E, q_G, q_H, q_I, q_L, q_K, q_M, q_F, q_P, q_S, q_T, q_W, q_Y, q_V
#     0     1     2    3   4   5    6    7    8    9    10   11   12   13   14  15   16   17   18   19    20   21   22
def get_prob(cnn_row, trans_row, aa):
  cnn_row = [float(item) for item in cnn_row]
  trans_row = [float(item) for item in trans_row]

  trans_aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
  cnn_aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
  trans_prob = 0
  cnn_prob = 0
  highest_prob = 0

  # trans row
  for i in range(20):
    if trans_aa_list[i] == aa:
      trans_prob = float(trans_row[i])

  # cnn row
  for i in range(20):
    if cnn_aa_list[i] == aa:
      cnn_prob = float(cnn_row[i])

  if trans_prob > cnn_prob:
    highest_prob = trans_prob
  else:
    highest_prob = cnn_prob

  return highest_prob

def get_predicted_aa(cnn_row, trans_row):
  cnn_row = [float(item) for item in cnn_row]
  trans_row = [float(item) for item in trans_row]

  trans_aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
  cnn_aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
  
  highest_trans_prob = 0.0
  highest_cnn_prob = 0.0

  trans_aa_pred = ""
  cnn_aa_pred = ""
 
  # getting trans predicted prob and amino acid
  for i in range(20):
    prob = trans_row[i]
  
    if prob > highest_trans_prob:
      highest_trans_prob = prob
      trans_aa_pred = trans_aa_list[i]

  # getting cnn predicted prob and amino acid
  for i in range(20):
    prob = cnn_row[i]
    if prob > highest_cnn_prob:
      highest_cnn_prob = prob
      cnn_aa_pred = cnn_aa_list[i]

  if highest_trans_prob > highest_cnn_prob:
    return trans_aa_pred
  else:
    return cnn_aa_pred

def get_row(cnn_row, trans_row):
  row = []
  cnn_row_probs = cnn_row[3:23]
  trans_row_probs = trans_row[4:24]
  aa_wt = trans_row[1]

  wt_prob = get_prob(cnn_row_probs, trans_row_probs, aa_wt )
  aa_pred = get_predicted_aa(cnn_row_probs, trans_row_probs)
  pred_prob = get_prob(cnn_row_probs, trans_row_probs, aa_pred)
  gene = trans_row[24]

  # new row:
  # 'position','wt_prob','aa_wt', 'pred_prob', 'aa_pred', gene
  row = [trans_row[2], wt_prob, trans_row[1], pred_prob, aa_pred, gene]

  return row



cnn_pred = "../../output/stats_cnn.csv"
trans_pred = "../../output/PSICOV_BERT_predictions.csv"
output_path = "../../output/ens_pred_highest_confidence.csv"

# trans data:
# row,aa_wt,position,wt_prob, A, C, D, E, F, G, H,  I,  K,  L,  M,  N,  P,  Q,  R,  S,  T,  V,  W,  Y, gene
# 0     1     2        3      4  5  6  7  8  9  10  11  12  13  14  15  16 17  18  19  20  21  22  23   24

# cnn data:
# position,gene,wt_aa,q_A,q_R,q_N, q_D, q_C, q_Q, q_E, q_G, q_H, q_I, q_L, q_K, q_M, q_F, q_P, q_S, q_T, q_W, q_Y, q_V
#     0     1     2    3   4   5    6    7    8    9    10   11   12   13   14  15   16   17   18   19    20   21   22

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['position','wt_prob','aa_wt', 'pred_prob', 'aa_pred', 'gene'])

  with open(cnn_pred, newline='') as csvfile:
    cnn_reader = list(csv.reader(csvfile, delimiter=','))
    with open(trans_pred, newline='') as csvfile:
      trans_reader = list(csv.reader(csvfile, delimiter=','))
      for cnn_row in cnn_reader:
        for trans_row in trans_reader:
          if cnn_row[3] != "q_A" and trans_row[1] != "aa_wt":
            if cnn_row[1] == trans_row[24] and cnn_row[0] == trans_row[2]: # check that gene and position is same
              new_row = get_row(cnn_row, trans_row)
              writer.writerow(new_row)

print()
print("Done! \nSaved output to:", output_path)

    



