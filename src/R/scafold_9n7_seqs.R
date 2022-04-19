library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)
library(seqinr)

# bert predictions for the 9n7 scaffold:

# loading data:
bert_9n7 <- read.csv(file = "./output/bert_9n7_nonresurfaced.csv", header=TRUE, sep=",")
esm1b_9n7 <- read.csv(file = "./output/esm1b9n7_nonresurfaced.csv", header=TRUE, sep=",")
cnn_9n7 <- read.csv(file = "./output/cnn_wt_max_freq_9n7.csv", header=TRUE, sep=",")

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

# lets find n-eff for the cnn data:
cnn_data_n <- cnn_9n7 %>%
  select(c(gene, position, A:Y)) %>%
  nest(data = c(A:Y)) %>%
  mutate(n_eff_cnn = map(data, get_neff)) %>%
  select(-c(data, gene))

#---------------------------------------------------------------------------------
#leaving only R and K
RK_bert <- bert_9n7 %>%
  select(c(position, K, R))

RK_esm1b <- esm1b_9n7 %>%
  select(c(position, K, R))

RK_cnn <- cnn_9n7 %>%
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
#--------------------------------------------------------------------------------
# clean cnn data:
cnn_data <- cnn_9n7 %>% # choose which transformer model to use here
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) 

cnn_data2 <- inner_join(cnn_data, RK_cnn) %>%
  select(c(gene, position, aa_wt, aa_pred, wt_prob, pred_prob, K, R)) %>%
  rename(aa_wt_cnn = aa_wt,
         aa_pred_cnn = aa_pred,
         wt_prob_cnn = wt_prob,
         pred_prob_cnn = pred_prob,
         K_cnn = K,
         R_cnn = R)

cnn_data2 <- inner_join(cnn_data2, cnn_data_n)

#---------------------------------------------------------------------------------
# combine all data together:
all_bert_esm <- inner_join(bert_data2, esm1b_data2) 
all_data <- left_join(all_bert_esm, cnn_data2) %>%
  rename(aa_wt = aa_wt_bert) %>%
  select(-c(aa_wt_esm1b, aa_wt_cnn)) %>%
  select(c(gene, position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn, wt_prob_bert, wt_prob_esm1b, wt_prob_cnn, pred_prob_bert, pred_prob_esm1b, pred_prob_cnn, K_bert, K_esm1b, K_cnn, R_bert, R_esm1b, R_cnn, n_eff_bert, n_eff_esm1b, n_eff_cnn))
  

all_data$n_eff_bert <- as.numeric(all_data$n_eff_bert)
all_data$n_eff_esm1b <- as.numeric(all_data$n_eff_esm1b)
all_data$n_eff_cnn <- as.numeric(all_data$n_eff_cnn)

#---------------------------------------------------------------------------------
# filtering out the CDRs:
all_data <- all_data %>%
  mutate(CDR = ifelse(position %in% c(29:35) | position %in% c(55:61) | position %in% c(101:109), "CDR", "non_CDR"))

all_data <- all_data %>%
  filter(CDR != "CDR")


#--------------------------------------------------------------------------------
# filter/rank

#filtered based on bert and esm1b (confident mispredictions)
# filtered_1 <- all_data %>%
#   filter(aa_pred_bert != aa_wt,
#          pred_prob_bert >= 0.7,
#          aa_pred_esm1b != aa_wt,
#          pred_prob_esm1b >= 0.7
#          )
# 
# # cnn confidence:
# filtered_1cnn <- all_data %>%
#   filter(aa_pred_cnn != aa_wt,
#          pred_prob_cnn >= 0.5
#   )

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

# filtered based on n-eff of cnn:
filtered_3cnn <- all_data %>%
  filter(aa_pred_cnn != aa_wt) %>%
  arrange(desc(n_eff_cnn)) %>%
  top_n(5, n_eff_cnn)

# ranking by K (bert):
filtered_4 <- all_data %>%
  filter(aa_pred_bert != aa_wt) %>% # mispredictions
  arrange(desc(K_bert)) %>%
  top_n(5, K_bert)

# ranking by K (esm):
filtered_5 <- all_data %>%
  filter(aa_pred_esm1b != aa_wt) %>% #mispredictions
  arrange(desc(K_esm1b)) %>%
  top_n(5, K_esm1b)

# ranking by K (cnn):
filtered_5cnn <- all_data %>%
  filter(aa_pred_cnn != aa_wt) %>% #mispredictions
  arrange(desc(K_cnn)) %>%
  top_n(5, K_cnn)

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

# ranking by R (cnn):
filtered_7cnn <- all_data %>%
  filter(aa_pred_cnn != aa_wt) %>% #mispredictions
  arrange(desc(R_cnn)) %>%
  top_n(10, R_cnn)

# residue mispredicted to R or K by bert or esm1b or cnn:

aa_pred <- c("C", "D", "R", "K")

predicted_KR <- function(aa_pred){
  if ("K" %in% aa_pred) {
    return("has_KR")
  }
  if ("R" %in% aa_pred){
    return("has_KR")
  }
  else {
    return("no_KR")
  }
}

print(map(aa_pred, predicted_KR))


filtered_8 <- all_data %>%
  filter(aa_wt != "K", aa_wt != "R") %>%
  nest(aa_pred = c(aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(KR_status = map_chr(aa_pred, predicted_KR)) %>%
  unnest(aa_pred) %>%
  filter(KR_status == "has_KR") %>%
  filter(!position %in% c(116, 15))

# residues with high R or K conf by bert or esm1b or cnn:
cut_off <- 0.20
filtered_9 <- all_data %>%
  filter(aa_wt != "K", aa_wt != "R") %>%
  filter(
    K_bert > cut_off |
    R_bert > cut_off | 
    K_esm1b > cut_off | 
    R_esm1b > cut_off | 
    K_cnn > cut_off | 
    R_cnn > cut_off)

# look at high CNN pred conf with low wt conf
filtered_10 <- all_data %>%
  filter(aa_wt != "K", aa_wt != "R") %>%
  filter(aa_pred_cnn == "K" | aa_pred_cnn == "R")


#------------------------------------------------------------------------------------
# saving seqs to fasta file:

# get the final decision residue at each position (in case the models contradict):
final_decision <- filtered_8 %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R")) %>%
  mutate(mutation = paste0(aa_wt, factor(position), new_residue))

# combine with old df:
final_decision <- left_join(all_data, final_decision)

# get new sequence:
final_decision <- final_decision %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq))

missing_A2K <- filtered_8 %>%
  filter(position != 2) %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R")) 
missing_A2K <- left_join(all_data, missing_A2K)
missing_A2K <- missing_A2K %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq))  

missing_Y62K <- filtered_8 %>%
  filter(position != 62) %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R"))
missing_Y62K <- left_join(all_data, missing_Y62K)
missing_Y62K <- missing_Y62K %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq)) 

missing_V4K <- filtered_8 %>%
  filter(position != 4) %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R"))
missing_V4K <- left_join(all_data, missing_V4K)
missing_V4K <- missing_V4K %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq)) 

missing_T96R <- filtered_8 %>%
  filter(position != 96) %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R"))
missing_T96R <- left_join(all_data, missing_T96R)
missing_T96R <- missing_T96R %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq)) 
 
A2K <- filtered_8 %>%
  filter(position == 2) %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R"))
A2K <- left_join(all_data, A2K)
A2K <- A2K %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq)) 

Y62K <- filtered_8 %>%
  filter(position == 62) %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R"))
Y62K <- left_join(all_data, Y62K)
Y62K <- Y62K %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq)) 

V4K <- filtered_8 %>%
  filter(position == 4) %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R"))
V4K <- left_join(all_data, V4K)
V4K <- V4K %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq)) 

T96R <- filtered_8 %>%
  filter(position == 96) %>%
  select(c(position, aa_wt, aa_pred_bert, aa_pred_esm1b, aa_pred_cnn)) %>%
  mutate(new_residue = ifelse(aa_pred_bert == "K" | aa_pred_esm1b == "K" | aa_pred_cnn == "K", "K", "R"))
T96R <- left_join(all_data, T96R)
T96R <- T96R %>%
  mutate(new_seq = ifelse(is.na(new_residue), aa_wt, new_residue)) %>%
  select(c(position, new_seq)) 



all_mutants_seq <- paste(final_decision$new_seq, collapse = "")

missing_A2K <- paste(missing_A2K$new_seq, collapse = "")
  
missing_Y62K <- paste(missing_Y62K$new_seq, collapse = "")
  
missing_V4K <- paste(missing_V4K$new_seq, collapse = "")
  
missing_T96R <- paste(missing_T96R$new_seq, collapse = "")
  
A2K <- paste(A2K$new_seq, collapse = "")

Y62K <- paste(Y62K$new_seq, collapse = "")

V4K <- paste(V4K$new_seq, collapse = "")

T96R <- paste(T96R$new_seq, collapse = "")

sequences <- list(all_mutants_seq, missing_A2K, missing_Y62K, missing_V4K, missing_T96R, A2K, Y62K, V4K, T96R)
names <- list("9n7_A2K_V4K_Y62K_T96R", "9n7_V4K_Y62K_T96R", "9n7_A2K_V4K_T96R", "9n7_A2K_V4K_Y62K", "9n7_A2K_Y62K_T96R", "9n7_A2K", "9n7_V4K", "9n7_Y62K", "9n7_T96R")
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# saving to a fasta file:
write.fasta(sequences, 
            names, 
            file = "./output/9n7_nn_variations.fasta", 
            open = "w", 
            nbchar = 123, 
            as.string = TRUE)
