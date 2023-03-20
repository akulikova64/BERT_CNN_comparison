## Sequence and structure based deep learning models represent different aspects of protein biochemistry

output/ - processed data

analysis/figures/ - contains all plots and figures

src/ - R and python scripts for data processing

model/ - contains the fully trained combined model


<br />
<br />
### Scripts for training and running the combined model:

1) get_train_val_test_data.py - compiles the training/test and validation data (concatenates esm1-b, protBERT and CNN outputs)

2) transfer_learning_network.py - neural network that trains the combined model

3) get_predictions_transf_learning.py - use this script to upload trained model and generate predicitons. 


<br />
<br />
### Scripts for data preparation:

1) align_fasta_to_struc_seqs.py - checks that the sequences are aligned to the structures via "position" (this is important when concatenating the different neural network outputs)


<br />

---

CC BY-NC 4.0 This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License (see https://creativecommons.org/licenses/by-nc/4.0/legalcode).

Using this work commercially may require a license to intellectual property owned by the Board of Regents of the University of Texas System. Contact licensing@otc.utexas.edu for inquiries.

---

