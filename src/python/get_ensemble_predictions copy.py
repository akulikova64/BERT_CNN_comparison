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
def get_prob(row, aa):
  aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

  for i in range(0, len(aa_list), 1):
    if aa_list[i] == aa:
      prob = float(row[i])

  return prob

def average_probs(cnn_row, trans_row):

  cnn_row = [float(item) for item in cnn_row]
  trans_row = [float(item) for item in trans_row]

  # TRANS row:
  #  A, C, D, E, F, G, H,  I,  K,  L,  M,  N,  P,  Q,  R,  S,  T,  V,  W,  Y
  #  0  1  2  3  4  5  6  7   8   9   10  11  12  13  14  15  16  17  18  19 
  
  # CNN row:
  # A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,  F,  P,  S,  T,  W,  Y,  V
  # 0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  

  A = (cnn_row[0] + trans_row[0]) / 2
  C = (cnn_row[4] + trans_row[1]) / 2
  D = (cnn_row[3] + trans_row[2]) / 2 
  E = (cnn_row[6] + trans_row[3]) / 2 
  F = (cnn_row[13] + trans_row[4]) / 2 
  G = (cnn_row[7] + trans_row[5]) / 2 
  H = (cnn_row[8] + trans_row[6]) / 2
  I = (cnn_row[9] + trans_row[7]) / 2 
  K = (cnn_row[11] + trans_row[8]) / 2 
  L = (cnn_row[10] + trans_row[9]) / 2 
  M = (cnn_row[12] + trans_row[10]) / 2 
  N = (cnn_row[2] + trans_row[11]) / 2 
  P = (cnn_row[14] + trans_row[12]) / 2 
  Q = (cnn_row[5] + trans_row[13]) / 2 
  R = (cnn_row[1] + trans_row[14]) / 2 
  S = (cnn_row[15] + trans_row[15]) / 2 
  T = (cnn_row[16] + trans_row[16]) / 2 
  V = (cnn_row[19] + trans_row[17]) / 2 
  W = (cnn_row[17] + trans_row[18]) / 2 
  Y = (cnn_row[18] + trans_row[19]) / 2

  new_probs = [A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y]

  return new_probs

def get_predicted_aa(row):
  aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

  aa_pred = ""
  current_higest_prob = 0

  for i in range(0, len(row), 1):
    prob = row[i]
    if prob > current_higest_prob:
      current_higest_prob = prob
      aa_pred = aa_list[i]

  return aa_pred

def get_row(cnn_row, trans_row):
  row = []
  
  averaged_row = average_probs(cnn_row[3:23], trans_row[4:24])
  wt_prob = get_prob(averaged_row, trans_row[1])
  aa_pred = get_predicted_aa(averaged_row)
  pred_prob = get_prob(averaged_row, aa_pred)
  gene = trans_row[24]

  # new row:
  # 'position','wt_prob','aa_wt', 'pred_prob', 'aa_pred', A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, gene
  row = [trans_row[2], wt_prob, trans_row[1], pred_prob, aa_pred] + averaged_row + [gene]

  return row

cnn_pred = "../../output/stats_cnn.csv"
trans_pred = "../../output/PSICOV_BERT_predictions.csv"
output_path = "../../output/ensemble_predictions.csv"

# trans data:
# row,aa_wt,position,wt_prob, A, C, D, E, F, G, H,  I,  K,  L,  M,  N,  P,  Q,  R,  S,  T,  V,  W,  Y, gene
# 0     1     2        3      4  5  6  7  8  9  10  11  12  13  14  15  16 17  18  19  20  21  22  23   24

# cnn data:
# position,gene,wt_aa,q_A,q_R,q_N, q_D, q_C, q_Q, q_E, q_G, q_H, q_I, q_L, q_K, q_M, q_F, q_P, q_S, q_T, q_W, q_Y, q_V
#     0     1     2    3   4   5    6    7    8    9    10   11   12   13   14  15   16   17   18   19    20   21   22

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['position','wt_prob','aa_wt', 'pred_prob', 'aa_pred','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','gene'])

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

    



