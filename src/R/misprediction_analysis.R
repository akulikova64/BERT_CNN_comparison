library(tidyverse)
library(cowplot)
library(data.table)
options(scipen = 999)

# this script will make a plot to compare the cnn and transformer to predict

#looking only at mispredictions, what is the frequency of the predicted amino acid in natural alignments (for both cnn and transformer)

trans_data <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")
cnn_data_all <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

align_data <- read.csv(file = "./output/stats_align_all.csv", header=TRUE, sep=",")

#clean cnn data:

cnn_data <- cnn_data_all %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) %>%
  rename(aa_pred_cnn = aa_predicted,
         aa_wt_cnn = aa_wt,
         freq_pred_cnn = freq_predicted,
         freq_wt_cnn = freq_wt)

#clean trans data:

transf_data <- trans_data %>%
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) %>%
  rename(freq_pred_transf = pred_prob,
         freq_wt_transf = wt_prob,
         aa_wt_transf = aa_wt,
         aa_pred_transf = aa_pred)

joined_data <- inner_join(cnn_data, transf_data)

# clean cnn data:
cnn_data_unjoined <- joined_data %>%
  select(-c(aa_wt_transf, aa_pred_transf, freq_wt_transf, freq_pred_transf)) %>%
  rename(freq_pred = freq_pred_cnn,
         aa_pred = aa_pred_cnn,
         freq_wt = freq_wt_cnn,
         aa_wt = aa_wt_cnn) %>%
  mutate(group = "cnn")

# clean tranformer data:

transf_data_unjoined <- joined_data %>%
  select(-c(aa_wt_cnn, aa_pred_cnn, freq_wt_cnn, freq_pred_cnn)) %>%
  rename(freq_pred = freq_pred_transf,
         aa_pred = aa_pred_transf,
         freq_wt = freq_wt_transf,
         aa_wt = aa_wt_transf) %>%
  mutate(group = "transformer")

# now bind rows:

joined2 <- rbind(cnn_data_unjoined, transf_data_unjoined)


#now we are adding alignment data:

align_data2 <- align_data %>%
  select(c(gene, position, q_A:q_V)) %>%
  nest(align_data = q_A:q_V)

joined3 <- inner_join(joined2, align_data2)

# get mispredictions only:
mispredictions <- joined3 %>%
  filter(aa_pred != aa_wt) %>%
  unnest(align_data) %>%
  mutate(pred_prob = pmax())

# count mispredictions: #these are the totals
count_mispr <- joined3 %>%
  mutate(mispredicted = aa_pred != aa_wt) %>%
  group_by(group) %>%
  count(mispredicted)

get_mispr_group <- function(wt_cnn, pred_cnn, wt_trans, pred_trans) {
  if (wt_cnn != pred_cnn & wt_trans != pred_trans) {
    return('both mispr')
  }
  else if (wt_cnn != pred_cnn & wt_trans == pred_trans) {
    return('cnn mispr')
  }
  else if (wt_cnn == pred_cnn & wt_trans != pred_trans) {
    return('trans mispr')
  }
  else if (wt_cnn == pred_cnn & wt_trans == pred_trans) {
    return('neither mispr')
  }
  NA_character_
 
}
  
joined3_wider <- joined3 %>%
  select(c(gene, position, aa_pred, aa_wt, group)) %>%
  pivot_wider(names_from = group, values_from = c(aa_pred, aa_wt)) %>%
  #wt_cnn, pred_cnn, wt_trans, pred_trans:
  mutate(mispredictions = pmap_chr(list(aa_wt_cnn, aa_pred_cnn, aa_wt_transformer,aa_pred_transformer), get_mispr_group)) 

count_mispr_results <- joined3_wider %>%
  group_by(mispredictions) %>%
  count()

only_both <- joined3_wider %>%
  filter(mispredictions == 'both mispr') %>%
  mutate(match = aa_pred_cnn == aa_pred_transformer)
  #group_by(match) %>%
  #count()






  