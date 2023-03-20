# Sequence and structure based deep learning models represent different aspects of protein biochemistry
## This repository contains the trained model, code, analysis and processed data.

output/ - processed data

data/ - unprocessed data

analysis/figures/ - contains all plots and figures

src/ - R and python scripts for data processing

training_data.txt - dataset used in the training of the combined model (not all proteins in this dataset were used in final training due to technical issues during CNN prediction generation or alignment of structure to sequence.)

model/ - contains the fully trained combined model (model__3306_70_epochs_run_1.h5)
<br />
<br />

### Scripts for training and running the combined model:

1) get_train_val_test_data.py - compiles the training/test and validation data (concatenates esm1-b, protBERT and CNN outputs)

2) transfer_learning_network.py - neural network that trains the combined model

3) get_predictions_transf_learning.py - use this script to upload trained model and generate predicitons

### Scripts for data preparation:

1) get_first_chain_of_CNN_predictions.py - only the first chain of the CNN predictions was used in the case of multimeric proteins. 

2) align_fasta_to_struc_seqs.py - checks that the sequences are aligned to the structures via "position" (important when concatenating the different neural network outputs)

3) cnn_output_to_csv.py - makes sure the CNN output is tidy and matches output from transformer models

---

CC BY-NC 4.0 This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License (see https://creativecommons.org/licenses/by-nc/4.0/legalcode).

Using this work commercially may require a license to intellectual property owned by the Board of Regents of the University of Texas System. Contact licensing@otc.utexas.edu for inquiries.

