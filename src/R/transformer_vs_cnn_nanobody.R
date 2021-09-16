library(tidyverse)
library(cowplot)
library(data.table)
options(scipen = 999)

# get the 

cnn_data_proc <- read.csv(file = "./output/cnn_wt_max_freq_3ogo.csv", header=TRUE, sep=",")
cnn_data_all <- read.csv(file = "./output/cnn_wt_max_freq_3ogo.csv", header=TRUE, sep=",")

trans_data <- read.csv(file = "./output/nanobody_3ogo.protBERT.csv", header=TRUE, sep=",")
