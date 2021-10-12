import os
import sys
import csv
import math
import numpy as np

# finds the proportion of mispredictions around each residue within a given radius.

def get_distance(coord1, coord2):
  coord1 = coord1[1:-1]
  coord2 = coord2[1:-1]

  c1 = coord1.split()
  c2 = coord2.split()

  x1, y1, z1 = float(c1[0]), float(c1[1]), float(c1[2])
  x2, y2, z2 = float(c2[0]), float(c2[1]), float(c2[2])

  distance = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

  return distance

def get_mispr_dict():
  # Save
  mispr_dict = {}

  with open("../../output/cnn_wt_max_freq.csv", newline='') as csvfile:
    reader1 = list(csv.reader(csvfile, delimiter=','))
    reader2 = reader1
    for row_r1 in reader1:
      for row_r2 in reader2:
        if row_r1[2] == row_r2[2] & row_r1[0] == row_r2[0]: # check that position and gene is same
          if row_r1[3] != row_r2[3]: # check that amino acids are mismatched (a misprediction)
            gene = row_r1[0]
            mispr_dict[gene].append(row_r1[2])

  np.save('mispr_dict.npy', mispr_dict) 

  return mispr_dict

def get_group(position, gene):
  mispr_positions = get_mispr_dict()

  if position in mispr_positions[gene]:
    group == "mispr"
  else:
    group == "corr"

  return group


sphere_radius = 5

input_path= "../../output/side_chain_centers/"
output_path = "../../output/sphere_densities_" + str(sphere_radius) + "A.csv"

geneList = os.listdir("../../data/pdb/")
coordList = os.listdir("../../output/side_chain_centers/")

# getting misprediction distances
with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['gene', 'position', 'group', 'num_correct', 'num_mispr', 'prop_mispr'])

  for file in coordList:
    gene = file[4:8]
    with open(input_path + file, newline='') as csvfile:
      reader1 = list(csv.reader(csvfile, delimiter=','))
      reader2 = reader1

      for row_r1 in reader1:
        position = row_r1[1]
        group = get_group(position, gene)
        mispr_count = 0
        corr_count = 0

        for row_r2 in reader2:
          if position != row_r2[1]: 

            distance = get_distance(row_r1[6], row_r2[6])
            if distance <= 5:
              contact_group = get_group(row_r2[1], gene)
              if contact_group == "mispr":
                mispr_count += 1
              elif contact_group == "corr":
                corr_count += 1

        writer.writerow([gene, position, group, corr_count, mispr_count, mispr_count/(corr_count + mispr_count)])


                  

      