library(tidyverse)
library(cowplot)
library(data.table)
options(scipen = 999)

get_aa_name <- function(x) {
  
  if (x == 'ALA') {
    return('A')
  }
  else if (x == 'ARG') {
    return('R')
  }
  else if (x == 'ASN') {
    return('N')
  }
  else if (x == 'ASP') {
    return("D")
  }
  else if (x == "CYS") {
    return("C")
  }
  else if (x == 'GLN') {
    return('Q')
  }
  else if (x == 'GLU') {
    return('E')
  }
  else if (x == 'GLY') {
    return('G')
  }
  else if (x == 'HIS') {
    return('H')
  }
  else if (x == 'ILE') {
    return('I')
  }
  else if (x == 'LEU') {
    return('L')
  }
  else if (x == 'LYS') {
    return('K')
  }
  else if (x == 'MET') {
    return('M')
  }
  else if (x == 'PHE') {
    return('F')
  }
  else if (x == 'PRO') {
    return('P')
  }
  else if (x == 'SER') {
    return('S')
  }
  else if (x == 'THR') {
    return('T')
  }
  else if (x == 'TRP') {
    return('W')
  }
  else if (x == 'TYR') {
    return('Y')
  }
  else if (x == 'VAL') {
    return('V')
  }
  NA_character_
}

get_mut <- function(x) {
  
  if (x == 8) {
    return('to K')
  }
  else if (x == 11) {
    return('to K')
  }
  else if (x == 12) {
    return('to R')
  }
  else if (x == 22) {
    return('to K')
  }
  else if (x == 70) {
    return('to K')
  }
  else if (x == 72) {
    return('to K')
  }
  else if (x == 83) {
    return('to R')
  }
  else if (x == 85) {
    return('to R')
  }
  else if (x == 86) {
    return('to K')
  }
  else {
    return('no change')
  }
}


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

cnn_data3 <- cnn_data3 %>%
  mutate(aa_wt = map_chr(aa_wt, get_aa_name))

data_bert2 <- data_bert %>%
  select(c(wtIndex, wtScore, wt, A:Y)) %>%
  rename(freq_wt = wtScore,
         aa_wt = wt,
         position = wtIndex)

data_bert3 <- data_bert2 %>%
  select(c(position, aa_wt, freq_wt, K, R)) %>%
  mutate(group = "transformer")
  
  
for_KR_plots <- rbind(data_bert3, cnn_data3)
for_KR_plots2 <- for_KR_plots %>% 
  #mutate(aa_wt = fct_relevel(aa_wt, "M","Q","V","Q","L","V","E","S","G","G","A","L","V","Q","P","G","G","S","L","R","L","S","C","A","A","S","G","F","P","V","N","R","Y","S","M","R","W","Y","R","Q","A","P","G","K","E","R","E","W","V","A","G","M","S","S","A","G","D","R","S","S","Y","E","D","S","V","K","G","R","F","T","I","S","R","D","D","A","R","N","T","V","Y","L","Q","M","N","S","L","K","P","E","D","T","A","V","Y","Y","C","N","V","N","V","G","F","E","Y","W","G","Q","G","T","Q","V","T","V","S","S"))
  mutate(wt_K = aa_wt == 'K',
         wt_R = aa_wt == 'R')


for_KR_plots3 <- for_KR_plots2 %>%
  mutate(mut = map_chr(position, get_mut))

labels <- c("M","Q","V","Q","L","V","E","S","G","G","A","L","V","Q","P","G","G","S","L","R","L","S","C","A","A","S","G","F","P","V","N","R","Y","S","M","R","W","Y","R","Q","A","P","G","K","E","R","E","W","V","A","G","M","S","S","A","G","D","R","S","S","Y","E","D","S","V","K","G","R","F","T","I","S","R","D","D","A","R","N","T","V","Y","L","Q","M","N","S","L","K","P","E","D","T","A","V","Y","Y","C","N","V","N","V","G","F","E","Y","W","G","Q","G","T","Q","V","T","V","S","S")

plot_K <- for_KR_plots2 %>%
  ggplot(aes(x = position, y = K, color = group)) +
  geom_line() +
  geom_point(aes(fill = wt_K),
             shape = 21, 
             size = 2) +
  scale_color_manual(values = c("#32a852", "#9d46b8")) +
  scale_fill_manual(values = c("black", "red")) +
  scale_x_continuous(
    name = "Position (wt amino acid)",
    labels = labels,
    breaks = seq(from = 1, to = 116, by = 1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Probability of Lysine (K)",
    limits = c(0, 1.0),
    breaks = seq(from = 0.0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  labs(color = "model", fill = "Is K the wt?") +
  theme_cowplot(9) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 9),
    panel.grid.minor = element_blank())

plot_K

# coloring in points by experimental data (ones that were really mutated):
plot_K_mut <- for_KR_plots3 %>%
  ggplot(aes(x = position, y = K, color = group)) +
  geom_rect(aes(ymin=0.0,ymax=1.0,xmin=26.5,xmax=34.5), fill= "#dedede", color = "#FFFFFF") +
  geom_rect(aes(ymin=0.0,ymax=1.0,xmin=51.5,xmax=60.5), fill= "#dedede", color = "#FFFFFF") +
  geom_rect(aes(ymin=0.0,ymax=1.0,xmin=97.5,xmax=104.5), fill= "#dedede", color = "#FFFFFF") +
  geom_line() +
  geom_point(aes(fill = mut),
             shape = 21, 
             size = 2) +
  scale_color_manual(values = c("#32a852", "#9d46b8")) +
  scale_fill_manual(values = c("black", "red", "yellow")) +
  scale_x_continuous(
    name = "Position (wt amino acid)",
    labels = labels,
    breaks = seq(from = 1, to = 116, by = 1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Probability of Lysine (K)",
    limits = c(0, 1.0),
    breaks = seq(from = 0.0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  labs(color = "model", fill = "Mutations \nin lab") +
  theme_cowplot(9) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 9),
    panel.grid.minor = element_blank())

plot_K_mut

#ggsave(filename = "./analysis/figures/3ogo_lysine_muts.png", plot = plot_K_mut, width = 12, height = 4)

#ggsave(filename = "./analysis/figures/3ogo_lysine.png", plot = plot_K, width = 12, height = 4)

# now lets make a plot for Arginine:

plot_R <- for_KR_plots2 %>%
  ggplot(aes(x = position, y = R, color = group)) +
  geom_line() +
  geom_point(aes(fill = wt_R),
             shape = 21, 
             size = 2) +
  scale_color_manual(values = c("#32a852", "#9d46b8")) +
  scale_fill_manual(values = c("black", "red")) +
  scale_x_continuous(
    name = "Position (wt amino acid)",
    labels = labels,
    breaks = seq(from = 1, to = 116, by = 1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Probability of Arginine (R)",
    limits = c(0, 1.0),
    breaks = seq(from = 0.0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  labs(color = "model", fill = "Is R the wt?") +
  theme_cowplot(9) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 9),
    panel.grid.minor = element_blank())

plot_R

plot_R_mut <- for_KR_plots3 %>%
  ggplot(aes(x = position, y = R, color = group)) +
  geom_rect(aes(ymin=0.0,ymax=1.0,xmin=26.5,xmax=34.5), fill= "#dedede", color = "#FFFFFF") +
  geom_rect(aes(ymin=0.0,ymax=1.0,xmin=51.5,xmax=60.5), fill= "#dedede", color = "#FFFFFF") +
  geom_rect(aes(ymin=0.0,ymax=1.0,xmin=97.5,xmax=104.5), fill= "#dedede", color = "#FFFFFF") +
  geom_line() +
  geom_point(aes(fill = mut),
             shape = 21, 
             size = 2) +
  scale_color_manual(values = c("#32a852", "#9d46b8")) +
  scale_fill_manual(values = c("black", "red", "yellow")) +
  scale_x_continuous(
    name = "Position (wt amino acid)",
    labels = labels,
    breaks = seq(from = 1, to = 116, by = 1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Probability of Arginine (R)",
    limits = c(0, 1.0),
    breaks = seq(from = 0.0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  labs(color = "model", fill = "Mutations \nin lab") +
  theme_cowplot(9) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 9),
    panel.grid.minor = element_blank())

plot_R_mut

#ggsave(filename = "./analysis/figures/3ogo_arginine_muts.png", plot = plot_R_mut, width = 12, height = 4)


#ggsave(filename = "./analysis/figures/3ogo_arginine.png", plot = plot_R, width = 12, height = 4)

#=====================================================================================
#now lets look at n-effective for both models  
#=====================================================================================

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

data_bert_n <- data_bert2 %>%
  select(c(position, aa_wt, A:Y)) %>%
  nest(data = c(A:Y)) %>%
  mutate(n_eff = map_chr(data, get_neff)) %>%
  select(-data) %>%
  mutate(group = "transformer")

data_cnn_n <- data_processed2 %>%
  mutate(aa_wt = map_chr(aa_wt, get_aa_name)) %>%
  nest(data = c(A:V)) %>%
  mutate(n_eff = map_chr(data, get_neff)) %>%
  select(-c(data, freq_wt)) %>%
  mutate(group = "cnn")

for_n_plot <- rbind(data_bert_n, data_cnn_n)

labels <- c("M","Q","V","Q","L","V","E","S","G","G","A","L","V","Q","P","G","G","S","L","R","L","S","C","A","A","S","G","F","P","V","N","R","Y","S","M","R","W","Y","R","Q","A","P","G","K","E","R","E","W","V","A","G","M","S","S","A","G","D","R","S","S","Y","E","D","S","V","K","G","R","F","T","I","S","R","D","D","A","R","N","T","V","Y","L","Q","M","N","S","L","K","P","E","D","T","A","V","Y","Y","C","N","V","N","V","G","F","E","Y","W","G","Q","G","T","Q","V","T","V","S","S")

for_n_plot$n_eff <- as.numeric(for_n_plot$n_eff)

for_n_plot2 <- for_n_plot %>%
  mutate(mut = map_chr(position, get_mut))

plot_n <- for_n_plot2 %>%
  ggplot(aes(x = position, y = n_eff, color = group)) +
  geom_rect(aes(ymin=1,ymax=16,xmin=26.5,xmax=34.5), fill= "#dedede", color = "#FFFFFF") +
  geom_rect(aes(ymin=1,ymax=16,xmin=51.5,xmax=60.5), fill= "#dedede", color = "#FFFFFF") +
  geom_rect(aes(ymin=1,ymax=16,xmin=97.5,xmax=104.5), fill= "#dedede", color = "#FFFFFF") +
  geom_line() +
  geom_point(aes(fill = mut),
             shape = 21, 
             size = 2) +
  scale_color_manual(values = c("#32a852", "#9d46b8")) +
  scale_fill_manual(values = c("black", "red", "yellow")) +
  scale_x_continuous(
    name = "Position (wt amino acid)",
    labels = labels,
    breaks = seq(from = 1, to = 116, by = 1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "N-effective",
    limits = c(1, 16),
    breaks = seq(from = 2, to = 16, by = 2),
    expand = c(0, 0)) +
  labs(color = "model") +
  theme_cowplot(9) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 9),
    panel.grid.minor = element_blank())

plot_n

ggsave(filename = "./analysis/figures/3ogo_n_eff_mut.png", plot = plot_n, width = 12, height = 4)


