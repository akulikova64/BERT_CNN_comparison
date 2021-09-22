library(tidyverse)
library(cowplot)
library(data.table)
options(scipen = 999)


data <- read.csv(file = "./output/cnn_wt_max_freq_3ogo.csv", header=TRUE, sep=",")

data2 <- data %>%
  filter(chain == 'G') %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) %>%
  rename(aa_pred_cnn = aa_predicted,
         aa_wt_cnn = aa_wt,
         freq_pred_cnn = freq_predicted,
         freq_wt_cnn = freq_wt)

data_bert <- read.csv(file = "./output/bert_prediction_3ogo.csv", header=TRUE, sep=",")

transf_stats <- data_bert %>%
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, X)) %>%
  rename(freq_pred_transf = pred_prob,
         freq_wt_transf = wtScore,
         aa_wt_transf = wt,
         aa_pred_transf = aa_pred,
         position = wtIndex)

data_cnn_raw <- read.csv(file = "./output/3ogo_final_tot.csv", header=TRUE, sep=",")

data_processed <- data_cnn_raw %>%
  filter(chain_id == 'G') %>%
  select(c(pos, wtAA, wt_prob)) %>%
  rename(position = pos,
         aa_wt = wtAA,
         freq_wt = wt_prob)
  
data_names <- data_cnn_raw %>%
  filter(chain_id == 'G') %>%
  select(c(pos, prALA:prVAL))

names(data_names) <-  c("position", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

data_processed2 <- inner_join(data_processed, data_names) 

cnn_data3 <- data_processed2 %>%
  select(c(position, aa_wt, freq_wt, K, R)) %>%
  mutate(group = "cnn")

data_bert2 <- data_bert %>%
  select(c(wtIndex, wtScore, wt, A:Y)) %>%
  rename(freq_wt = wtScore,
         aa_wt = wt,
         position = wtIndex)

data_bert3 <- data_bert2 %>%
  select(c(position, aa_wt, freq_wt, K, R)) %>%
  mutate(group = "transformer")
  
  
  

for_KR_plots <- rbind(data_bert3, cnn_data3)


  



