import json
import sys
import os
import numpy as np
#import tables

# this script gets names of proteins for training combinded net.


# paths 
rscb_all_path = "../../data/transfer_learning_net/new_training_data/rcsb_all_names.txt"
cnn_train_path = "../../data/transfer_learning_net/new_training_data/cnn_training_pdbs.txt"
psicov_names_path = "../../data/pdb/"
output_list_path = "../../data/transfer_learning_net/new_training_data/filtered_training_data.txt"

with open(output_list_path, 'w') as output_file:
  with open(rscb_all_path) as all_file:
    with open(cnn_train_path) as cnn_train_path:

      # collecting the cnn training pdb names:
      cnn_train_pdbs = []
      dicts = cnn_train_path.readlines()
      for dict in dicts:
        dict_json = json.loads(dict)
        cnn_train_pdbs.append(dict_json["pdb_code"].upper())
      print("Finished collecting old CNN training data names.")
      
      # collecting the psicov pdb names:
      psicov = []
      files = os.listdir(psicov_names_path)
      for file in files:
        psicov.append(file[0:4].upper())
      print("Finished collecting PSICOV pdb names")

      # parcing through the rcsb_all file and collecting pdbs that are not
      # found in the cnn training data of the PSICOV set.

      for line in all_file:
        split_line = line.split(" ")
        pdb = split_line[0][0:4]
        if pdb not in psicov and pdb not in cnn_train_pdbs:
          output_file.write(str(pdb) + "\n")
    

  print("Finished making training pdb list!")


        
    

      

  
  


