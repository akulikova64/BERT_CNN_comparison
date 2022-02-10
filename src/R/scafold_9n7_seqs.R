library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)

# bert predictions for the 9n7 scaffold:

# loading data:
bert_9n7 <- read.csv(file = "./output/bert_9n7_nonresurfaced.csv", header=TRUE, sep=",")
esm1b_9n7 <- read.csv(file = "./output/esm1b9n7_nonresurfaced.csv", header=TRUE, sep=",")

#-------------------------------------------------------------------------------
# getting n-eff=:

get_neff <- function(freq_list) {
  sum <- 0
  for (freq in freq_list)
  {
    freq <- as.numeric(freq)
    sum <- sum + freq * ifelse(freq == 0, 1, log(freq))
  }
  
  entropy <- -sum
  neff <- exp(entropy)
  return(as.numeric(neff))
}

# lets find n-eff for the bert data:
bert_data_n <- bert_9n7 %>%
  select(c(gene, position, A:Y)) %>%
  nest(data = c(A:Y)) %>%
  mutate(n_eff_bert = map(data, get_neff)) %>%
  select(-c(data, gene))

# lets find n-eff for the esm1b data:
esm_data_n <- esm1b_9n7 %>%
  select(c(gene, position, A:Y)) %>%
  nest(data = c(A:Y)) %>%
  mutate(n_eff_esm1b = map(data, get_neff)) %>%
  select(-c(data, gene))

#---------------------------------------------------------------------------------
#leaving only R and K
RK_bert <- bert_9n7 %>%
  select(c(position, K, R))

RK_esm1b <- esm1b_9n7 %>%
  select(c(position, K, R))

#----------------------------------------------------------------------------------
# clean bert data:
bert_data <- bert_9n7 %>% # choose which transformer model to use here
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) 

bert_data2 <- inner_join(bert_data, RK_bert) %>%
  select(c(gene, position, aa_wt, aa_pred, wt_prob, pred_prob, K, R)) %>%
  rename(aa_wt_bert = aa_wt,
         aa_pred_bert = aa_pred,
         wt_prob_bert = wt_prob,
         pred_prob_bert = pred_prob,
         K_bert = K,
         R_bert = R)

bert_data2 <- inner_join(bert_data2, bert_data_n)

#--------------------------------------------------------------------------------
# clean esm1b data:
esm1b_data <- esm1b_9n7 %>% # choose which transformer model to use here
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) 

esm1b_data2 <- inner_join(esm1b_data, RK_esm1b) %>%
  select(c(gene, position, aa_wt, aa_pred, wt_prob, pred_prob, K, R)) %>%
  rename(aa_wt_esm1b = aa_wt,
         aa_pred_esm1b = aa_pred,
         wt_prob_esm1b = wt_prob,
         pred_prob_esm1b = pred_prob,
         K_esm1b = K,
         R_esm1b = R)

esm1b_data2 <- inner_join(esm1b_data2, esm_data_n)

#---------------------------------------------------------------------------------
# combine all data together:
all_data <- inner_join(bert_data2, esm1b_data2) %>%
  rename(aa_wt = aa_wt_bert) %>%
  select(-aa_wt_esm1b) %>%
  select(c(gene, position, aa_wt, aa_pred_bert, aa_pred_esm1b, wt_prob_bert, wt_prob_esm1b, pred_prob_bert, pred_prob_esm1b, K_bert, K_esm1b, R_bert, R_esm1b, n_eff_bert, n_eff_esm1b))
  
all_data$n_eff_bert <- as.numeric(all_data$n_eff_bert)
all_data$n_eff_esm1b <- as.numeric(all_data$n_eff_esm1b)

#--------------------------------------------------------------------------------
# filter/rank

#filtered based on bert and esm1b (confident mispredictions)
filtered_1 <- all_data %>%
  filter(aa_pred_bert != aa_wt,
         pred_prob_bert >= 0.7,
         aa_pred_esm1b != aa_wt,
         pred_prob_esm1b >= 0.7
         )

# filtered based on n-eff of bert:
filtered_2 <- all_data %>%
  filter(aa_pred_bert != aa_wt) %>%
  arrange(desc(n_eff_bert)) %>%
  top_n(5, n_eff_bert)

# filtered based on n-eff of esm1b:
filtered_3 <- all_data %>%
  filter(aa_pred_esm1b != aa_wt) %>%
  arrange(desc(n_eff_esm1b)) %>%
  top_n(5, n_eff_esm1b)

# ranking by K (bert):
filtered_4 <- all_data %>%
  filter(aa_pred_bert != aa_wt) %>% #mispredictions
  arrange(desc(K_bert)) %>%
  top_n(5, K_bert)

# ranking by K (esm):
filtered_5 <- all_data %>%
  filter(aa_pred_esm1b != aa_wt) %>% #mispredictions
  arrange(desc(K_esm1b)) %>%
  top_n(5, K_esm1b)

# ranking by R (bert):
filtered_6 <- all_data %>%
  filter(aa_pred_bert != aa_wt) %>% #mispredictions
  arrange(desc(R_bert)) %>%
  top_n(5, R_bert)

# ranking by R (esm):
filtered_7 <- all_data %>%
  filter(aa_pred_esm1b != aa_wt) %>% #mispredictions
  arrange(desc(R_esm1b)) %>%
  top_n(5, R_esm1b)

# residue mispredicted to R or K by bert or esm1b:
filtered_8 <- all_data %>%
  filter(aa_wt != "K",
         aa_wt != "R",
         ) %>%
  filter(aa_pred_bert == "K" | aa_pred_bert == "R" | aa_pred_esm1b == "K" | aa_pred_esm1b == "R")

# residues with high R or K conf by bert or esm1b:
filtered_9 <- all_data %>%
  filter(aa_wt != "K",
         aa_wt != "R",
  ) %>%
  filter(K_bert > 0.18 | R_bert > 0.18 | K_esm1b > 0.18 | R_esm1b > 0.18)


#------------------------------------------------------------------------------------
# saving seqs to fasta file:

high_KR_conf <- filtered_9 %>%
  select(c(position, aa_pred_bert, aa_pred_esm1b)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K", "K", "R"))

# combine with old df:
high_KR_conf <- left_join(all_data, high_KR_conf)

# get new sequence:
high_KR_conf <- high_KR_conf %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq))

high_KR_conf_seq <- high_KR_conf$new_seq

#--------------------------------------------------------------------------------------
# saving to a fasta file:
write.fasta(sequences, 
            names, 
            file = "./output/9n7_nn_variations.fasta", 
            open = "w", 
            nbchar = 123, 
            as.string = TRUE)
