import json
import sys
import os
import numpy as np
import tables

# this script gets names of proteins for training combinded net.


# paths 
rscb_all_path = "../../data/transfer_learning_net/new_training_data/rcsb_all_names.txt"
cnn_train_path = "../../data/transfer_learning_net/new_training_data/cnn_training_pdbs.txt"
psicov_names_path = "../../data/pdb/"
output_list_path = "../../data/transfer_learning_net/new_training_data/filtered_training_data.txt"
output_unsorted = "../../data/transfer_learning_net/new_training_data/unsorted_data.txt"

with open(output_unsorted, "a", newline='') as unsorted_output:
  with open(output_list_path, 'w', newline='') as output_file:
    with open(rscb_all_path) as all_file:
      with open(cnn_train_path) as cnn_train_path:

        # collecting the cnn training pdb names:
        cnn_train_pdbs = []
        dicts = cnn_train_path.readlines()
        for dict in dicts:
          dict_json = json.loads(dict)
          cnn_train_pdbs.append(dict_json["pdb_code"].upper())
        
        # collecting the psicov pdb names:
        psicov = []
        files = os.listdir(psicov_names_path)
        for file in files:
          psicov.append(file[0:4].upper())

        # parcing through the rcsb_all file and collecting pdbs that are not
        # found in the cnn training data of the PSICOV set.

        # making a list of pdb lists, where each list is contains sequences that are 50% or more similar to each other.
        final_list = []
        similarity_list = []
        for line in all_file:
          split_line = line.split(" ")
          for pdb_code in split_line:
            pdb = pdb_code[0:4]
            if pdb not in psicov and pdb not in similarity_list and pdb not in cnn_train_pdbs:
              similarity_list.append(pdb)
          
          # appending similarity list to final "filtered" list
          unsorted_output.write(str(similarity_list))
          final_list.append(similarity_list)
          similarity_list = []
        
    # final_list structure: final_list[similarity index (i)][index within similarity (j)]
    len_list = len(final_list) # index of row (a specific similarity)
    
    j = 0 # j is index of col (index within similarity group)
    while True:
      num_records = 0
      for i in range(1):
        try:
          output_file.write(str(final_list[i][j]) + "\n")
          num_records += 1
        except: # except out of bound error
          continue

      if num_records == 0:
        break

      j += 1
    


        
    

      

  
  


