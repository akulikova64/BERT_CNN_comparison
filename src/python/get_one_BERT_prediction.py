import os
import csv
import sys
from Bio import SeqIO
import pandas as pd
from berteome import berteome

# this script take in one sequence and makes predictions for each position. 

sequence = "MQVQLVESGGALVQPGGSLRLSCAASGFPVNRYSMRWYRQAPGKEREWVAGMSSAGDRSSYEDSVKGRFTISRDDARNTVYLQMNSLKPEDTAVYYCNVNVGFEYWGQGTQVTVSS"
output_path = "../../output/bert_prediction_3ogo.csv"

#with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  #writer = csv.writer(CSV_file)
  #writer.writerow(['row', 'aa_wt','position','wt_prob','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])

all_predictions = berteome.allResiduePredictions(sequence)
all_predictions_DF = berteome.residuePredictionScore(all_predictions, sequence)

all_predictions_DF.to_csv(output_path, header = True) # mode = 'a'