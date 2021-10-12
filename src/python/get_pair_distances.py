import os
import sys
import csv
import math

# get the distance between all residue pairs for each protein and save each pair to CSV.

# first, I will get two csv's with correctly and incorrectly predicted positions in R.
# next, I will find the distances per gene within both groups (mispredicted and correctly predicted positions)

def get_distance(coord1, coord2):
  coord1 = coord1[1:-1]
  coord2 = coord2[1:-1]

  c1 = coord1.split()
  c2 = coord2.split()

  x1, y1, z1 = float(c1[0]), float(c1[1]), float(c1[2])
  x2, y2, z2 = float(c2[0]), float(c2[1]), float(c2[2])

  distance = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

  return distance


input_path_1 = "../../output/mispredictions.csv"
input_path_2 = "../../output/correct_pred.csv"

output_path_1 = "../../output/mispr_pair_dist.csv"
output_path_2 = "../../output/correct_pair_dist.csv"


geneList = os.listdir("../../data/pdb/")


# getting misprediction distances
with open(output_path_2, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['gene', 'position_1', 'position_2', 'aa_1', 'aa_2', 'distance'])

  with open(input_path_2, newline='') as csvfile:
    reader1 = list(csv.reader(csvfile, delimiter=','))
    reader2 = reader1

    for gene in geneList:
      for row_r1 in reader1:
        if row_r1[1] == 'gene':
          continue
        else:
          for row_r2 in reader2:
            if row_r2[1] == 'gene':
              continue
            elif row_r2[1] != 'gene':
              if row_r1[1] == gene[0:4] and row_r2[1] == gene[0:4]:
              
                distance = get_distance(row_r1[7], row_r2[7])
                writer.writerow([gene[0:4], row_r1[2], row_r2[2], row_r1[4], row_r2[4], distance])

      