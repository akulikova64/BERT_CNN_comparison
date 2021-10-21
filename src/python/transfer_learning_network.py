
# neural network that takes in cnn and transformer output and outputs final combined amino acid probabilites.
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from datetime import datetime
import time
import numpy as np
import collections
import math
import csv
import sys
import re
import os

#try:
import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.optimizers import Adam
from keras.callbacks import Callback
from keras.models import load_model
from keras.utils import to_categorical
from keras.callbacks import ModelCheckpoint
from keras.callbacks import CSVLogger

'''
except ImportError:
  import tensorflow
  from tensorflow.keras.models import Sequential
  from tensorflow.keras.layers import Dense
  from tensorflow.keras.callbacks import Callback
  from tensorflow.keras.optimizers import Adam
  from tensorflow.keras.models import load_model
  from tensorflow.keras.utils import to_categorical
  from tensorflow.keras.callbacks import ModelCheckpoint
  from tensorflow.keras.callbacks import CSVLogger
'''

def timestamp():
  return str(datetime.now().time())

def get_current_run(model_id, output_path):
  """ checks output folder for the previous run number and adds one """

  path = output_path + "/"
  fileList = os.listdir(path)

  current_run = 1
  for file in fileList:
    match_1 = re.search(r'model_' + model_id + '_', file)
    if match_1:
      match_2 = re.search(r'run_([0-9]+).h5', file)
      if match_2:
        run = int(match_2.group(1))
        if run >= current_run:
          current_run = run + 1

  return current_run

def compile_model(run, loss, optimizer, metrics, output_path):
  """ loads previous model or compiles a new model """

  last_model = output_path + "/model_" + model_id + "_run_" + str(run-1) + ".h5"
  try:
    model = load_model(last_model)
    print("Loaded last model:", last_model)
  except:
    model = get_model()
    model.compile(loss = loss, optimizer = optimizer, metrics = metrics)
    print("No previous model found. Training as Run 1 (from zero) :")

  print("Model loaded/compiled:", timestamp(), "\n")

  return model

def load_data(training_path, validation_path):
  """ loads training and validation data """

  print("\nStarting to load training data:", timestamp())
  x_train = np.load(training_path + "nn_data_train.npy", allow_pickle = True).tolist()
  y_train = np.load(training_path + "answers_train.npy", allow_pickle = True).tolist()
  print("Finished loading training data:", timestamp())
  x_val = np.load(validation_path + "nn_data_val.npy", allow_pickle = True).tolist()
  y_val = np.load(validation_path + "answers_val.npy", allow_pickle = True).tolist()
  print("Finished loading validation data:", timestamp())

  return x_train, y_train, x_val, y_val

def get_model(GPUS = 1):
  """ model with one dense layer """

  model = Sequential()
  model.add(Dense(1000, activation = 'relu')) # 500 nodes in the last hidden layer
  model.add(Dense(20, activation = 'softmax')) # output layer has 20 possible classes (amino acids 0 - 19)

  return model

def timestamp():
  return str(datetime.now().time())

def get_history(model_id, output_path):
  """ parces the history CSV file """ 

  accuracy, loss, val_accuracy, val_loss = [], [], [], []

  with open(output_path + "/model_" + str(model_id) + "_history_log.csv") as hist_file:
    csv_reader = csv.DictReader(hist_file, delimiter=',')
    for row_values in csv_reader:
      accuracy.append(float(row_values['accuracy']))
      loss.append(float(row_values['loss']))
      val_accuracy.append(float(row_values['val_accuracy']))
      val_loss.append(float(row_values['val_loss']))
  
  return accuracy, loss, val_accuracy, val_loss

def get_plots(run, model_id, BLUR, loss, optimizer, learning_rate, data, output_path, rotations):
  """ creates simple plots of accuracy and loss for training and validation """

  parameter_text = "BLUR = " + str(BLUR) + "\n" + \
                   "loss = " + str(loss) + "\n" \
                   "optimizer = " + str(optimizer)[18:28] + ". \n" \
                   "learning rate = " + str(learning_rate) + "\n" \
                   "training data = " + str(data) + "\n" \
                   "rotations = " + str(rotations) 

  timestr = time.strftime("%m%d-%H%M%S")

  accuracy, loss, val_accuracy, val_loss = get_history(model_id, output_path)

  print("Making plots: ", timestamp(), "\n")

  # plotting accuracy 
  plt.plot(accuracy)
  plt.plot(val_accuracy)
  plt.title('Model ' +  model_id  +  ' Accuracy')
  plt.suptitle(str(datetime.now()), size = 7)
  plt.ylabel('accuracy')
  plt.xlabel('epoch')
  plt.legend(['training', 'validation'], loc = 'upper left')
  plt.annotate(parameter_text, xy = (0.28, 0.84), xycoords = 'axes fraction', size = 7) 
  plt.savefig(output_path + "/Accuracy_model_" + model_id + "_run_" + str(run) + "_" + timestr + ".pdf")
  plt.clf()

  # plotting loss
  plt.plot(loss)
  plt.plot(val_loss)
  plt.title('Model ' +  model_id  + ' Loss')
  plt.suptitle(str(datetime.now()), size = 7)
  plt.ylabel('loss')
  plt.xlabel('epoch')
  plt.legend(['training', 'validaton'], loc = 'upper left')
  plt.annotate(parameter_text, xy = (0.28, 0.84), xycoords = 'axes fraction', size = 7)
  plt.savefig(output_path + "/loss_model_" + model_id + "_run_" + str(run) + "_" + timestr + ".pdf")

  #saving all data in a CSV file
  path = output_path + "/model_" + model_id + "_run_" + str(run) + "_" + timestr + ".csv"
  print("Starting to write CSV file:", timestamp())
  with open(path, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Model_ID", "Epoch", "BLUR", "learning_rate", "rotations", "Acc_train", "Acc_val", "Loss_train", "Loss_val"])
    for i in range(0, len(accuracy)):
      writer.writerow([model_id, i+1, BLUR, learning_rate, rotations, accuracy[i], val_accuracy[i], loss[i], val_loss[i]])
  print("Finished writing CSV file:", timestamp())

def save_model(model, model_id, run, output_path):
  """ saves model with weights """

  current_model = output_path + "/model_" + model_id + "_run_" + str(run) + ".h5"
  model.save(current_model)
  print("Saved current model:", timestamp(), "\n")

# training and saving the model
def train_model(model, run, epochs, x_train, y_train, x_val, y_val, output_path):
  """ calling the model to train """

  print("Starting to train:", timestamp(), "\n")
  
  model.fit(
      x = x_train,
      y = y_train,
      validation_data = (x_val, y_val),
      epochs = epochs, 
      verbose = 1)

  print("Finished training and validation:", timestamp(), "\n")
  
  save_model(model, model_id, run, output_path)
  print("Saved model:", timestamp(), "\n")
  print(model.summary(), "\n")


def get_val_predictions(model, model_id, run, x_val, output_path):
  """ generating validation predicitons """

  print("Starting to predict:", timestamp(), "\n")
  predictions = model.predict_classes(x_val, verbose = 0)
  timestr = time.strftime("%m%d-%H%M%S")
  np.save(output_path + "/predictions_model_" + model_id + "_run_" + str(run) + "_" + timestr + ".npy", predictions)
  print("Finished predicting:", timestamp(), "\n")




#========================================================================================================
# Setting the variables, parameters and data paths/locations:
#========================================================================================================
### data paths/locations
training_path = "../../data/transfer_learning_net/input/training/"
validation_path = "../../data/transfer_learning_net/input/validation/"
output_path = "../../data/transfer_learning_net/output/training_results/"

### variables
EPOCHS = 2 # iterations through the data

model_id = "1"
learning_rate = 0.0001
run = get_current_run(model_id, output_path)

### setting parameters for training
loss ='categorical_crossentropy'
optimizer = Adam(lr = learning_rate)
metrics = ['accuracy']
metrics = ['accuracy']

#========================================================================================================
# Training the network:
#========================================================================================================

### loading training and validation data
x_train, y_train, x_val, y_val = load_data(training_path, validation_path)

### compiling the model
model = compile_model(model_id, run, loss, optimizer, metrics, output_path)

### training and validation
train_model(model, model_id, run, EPOCHS, x_train, y_train, x_val, y_val, output_path)

### generating validation predictions
get_val_predictions(model, model_id, run, x_val, output_path)

### results
get_plots(run, model_id, loss, optimizer, learning_rate, training_path[3:-1], output_path)

print("Training completed!")
