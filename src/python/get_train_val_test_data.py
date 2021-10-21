import os
import sys
import csv
import math
import numpy as np

# This program gets the training data for the transfer learining network.

# data paths
cnn_pred = "../../output/PSICOV_CNN_output/"
trans_pred = "../../output/PSICOV_BERT_predictions.csv"

train_path = "../../data/transfer_learning_net/input/training/"
val_path = "../../data/transfer_learning_net/input/validation/"
test_path = "../../data/transfer_learning_net/input/testing/"

