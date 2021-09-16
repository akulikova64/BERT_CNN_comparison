library(tidyverse)
library(cowplot)
library(data.table)
options(scipen = 999)

# this script compares BERT to CNN for accuracy in predictions

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
  

match_wt <- joined2 %>%
  mutate(match_predict_wt = aa_pred == aa_wt)

#data entries where the predicted amino acid matches the wt
stats_1 <- match_wt %>%
  group_by(gene, group) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)))

means <- stats_1 %>%
  group_by(group) %>%
  summarise(mean = mean(freq_predict_wt))
#cnn: 75.8% Accuracy

stats_1 <- stats_1 %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(freq_predict_wt * (group == "cnn")),
    dx = rnorm(n(), mean = 0, sd = .05),
    dy = rnorm(n(), mean = 0, sd = .05),
    x_value = as.numeric(factor(group)))  

plot <- stats_1 %>%
  ggplot(aes(x = group, y = freq_predict_wt)) +
  geom_path(
    aes(x = as.numeric(factor(group))+dx, y = freq_predict_wt+dy, group = gene, color = color_y),
    size = 0.25) +
  geom_point(
    aes(x = as.numeric(factor(group))+dx, y = freq_predict_wt+dy, group = gene, fill = color_y),
    shape = 21, 
    color = "black",
    size = 2) +
  scale_x_continuous(
    name = "Neural Network",
    limits = c(0.7,2.3),
    labels = c("CNN", "Transformer"),
    breaks = (seq(from = 1.0, to = 2.0, by = 1))
    ) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.3, 1.0),
    breaks = seq(from = 0.3, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_color_gradient(
      aesthetics = c("color", "fill"), 
      high = "#ffd966", 
      low = "#000066") +
  theme_bw(16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 16),
    panel.grid.minor = element_blank())

plot

ggsave(filename = paste0("./analysis/figures/accuracy_cnn_bert.png"), plot = plot, width = 6, height = 6)
 
#==========================================================================================================
# now lets look at n-eff correlations. CNN vs. Transformer.
#==========================================================================================================

get_neff <- function(freq_list) {
  sum = 0
  for (freq in freq_list)
  {
    freq = as.numeric(freq)
    sum = sum + freq * ifelse(freq == 0, 1, log(freq))
  }
    
  entropy = -sum
  neff = exp(entropy)
  return(neff)
}

# lets find n-eff for the transformer data:


trans_data2 <- trans_data %>%
  select(c(gene, position, A:Y)) %>%
  mutate(n_eff = map(trans_data[A:Y], get_neff))





# preparing alignment data:
align_data2 <- align_data %>%
  select(c(position, gene, n_eff)) %>%
  rename(n_eff_align = n_eff)












