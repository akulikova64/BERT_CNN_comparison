import os
import csv
import sys
from Bio import SeqIO
import pandas as pd
from berteome import berteome

# this script takes in a fasta file of sequences and runds every sequence through the BERT transformer model. 

input_path = "../../data/PSICOV_seqs_reduced.fasta"
output_path = "../../output/PSICOV_BERT_predictions.csv"

'''
with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['row', 'aa_wt','position','wt_prob','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y', 'gene'])
'''

records = list(SeqIO.parse(input_path, "fasta"))
for rec in records:
  sequence = str(rec.seq)
  gene = str(rec.name)

  all_predictions = berteome.allResiduePredictions(sequence)
  all_predictions_DF = berteome.residuePredictionScore(all_predictions, sequence)
  all_predictions_DF['gene'] = gene

  all_predictions_DF.to_csv(output_path, mode = 'a', header = False)

