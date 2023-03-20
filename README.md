## Sequence and structure based deep learning models represent different aspects of protein biochemistry

#### Scripts for training and running the combined model:

get_train_val_test_data.py - compiles the training/test and validation data (concatenates esm1-b, protBERT and CNN outputs)

transfer_learning_network.py - neural network that trains the combined model

get_predictions_transf_learning.py - use this script to upload trained model and generate predicitons. 

