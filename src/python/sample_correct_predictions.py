import os
import sys
import csv
import math
import numpy as np
import random

# sample the sphere_densities.csv files so that each gene has the same number of mispredictions and 
# correct predictions. The correct predictions need to be sampled at random for each gene. 
# this will only be done for the first file (5A). The other files will be sampled based on the positions
# in the 5A file for consistency

output_path= "../../output/equal_5A.csv"
input_path = "../../output/sphere_densities_5A.csv"
geneList = os.listdir("../../data/pdb/")

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['gene', 'position', 'group', 'num_correct', 'num_mispr', 'prop_mispr'])
  
  for file in geneList:
    gene = file[0:4]
    mispr_count = 0
    correct_positions = []

    with open(input_path, newline = '') as csvfile:
      reader1 = csv.reader(csvfile, delimiter = ',')
      for row in reader1:
        if row[0] == gene:
          if row[2] == "mispr":
            mispr_count += 1
          elif row[2] == "corr":
            correct_positions.append(row[1])
      
      random_rows = random.sample(correct_positions, mispr_count)

    with open(input_path, newline = '') as csvfile:
      reader2 = csv.reader(csvfile, delimiter = ',')
      for row in reader2:
        if row[0] == gene:
          if row[2] == "mispr":
            writer.writerow(row)
          elif row[2] == "corr" and row[1] in random_rows:
            writer.writerow(row)

print()
print("Done! \nOutput file saved to:", output_path)   
            


      
