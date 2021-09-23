library(tidyverse)
library(cowplot)
library(data.table)
library(ggforce)
library(colorspace)
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

# total_aa <- cnn_data %>%
#   group_by(aa_wt_cnn) %>%
#   count()
# 
# write.csv(total_aa, "./output/aa_counts_PSICOV.csv")

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
#transformer: 71.6% accuracy

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
    limits = c(0.1, 1.1),
    breaks = seq(from = 0.1, to = 1.0, by = 0.1),
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

ggsave(filename = paste0("./analysis/figures/accuracy_cnn_bert_150.png"), plot = plot, width = 6, height = 6)
 
#========================================================================
#lets look at CNN confidence vs. Accuracy for both models:
get_pred_bin <- function(x) {
  
  if (x > 0 & x <= 0.2) {
    return("(0-0.2]")
  }
  else if (x > 0.2 & x <= 0.4) {
    return("(0.2-0.4]")
  }
  else if (x > 0.4 & x <= 0.6) {
    return("(0.4-0.6]")
  }
  else if (x > 0.6 & x <= 0.8) {
    return("(0.6-0.8]")
  }
  else if (x > 0.8 & x <= 1.0) {
    return("(0.8-1.0]")
  }
  NA_character_
}


joined2$freq_pred <- as.numeric(joined2$freq_pred)

joined_conf <- joined2 %>%
  mutate(pred_bin = map_chr(freq_pred, get_pred_bin))

 
  
match_wt <- joined_conf %>%
  mutate(match_predict_wt = aa_pred == aa_wt)

#data entries where the predicted amino acid matches the wt
stats_1 <- match_wt %>%
  group_by(gene, group, pred_bin) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt))) %>%
  na.omit()

stat_data_1 <- stats_1 %>%
  select(-gene) %>%
  group_by(pred_bin, group) %>%
  summarise(estimate = mean(freq_predict_wt),
            std_error = sd(freq_predict_wt)/sqrt(length(freq_predict_wt))) %>%
  na.omit()

# fill = "#8c7b9d", color = "#655775"



plot_conf <- stats_1 %>%
  ggplot(aes(y = freq_predict_wt, x = pred_bin, fill = group)) +
  geom_violin(alpha = 0.7, size = 0.5, bw = 0.03, position = position_dodge(width = 0.8)) + 
  geom_pointrange(data = stat_data_1, aes(x = pred_bin,
                                          y = estimate,
                                          ymin = estimate - 1.96*std_error,
                                          ymax = estimate + 1.96*std_error),
                  color = "black", 
                  alpha = 0.6, 
                  size = 0.3, 
                  position = position_dodge(width = 0.8)) +
  geom_sina(aes(color = group), 
            size = 0.5,
            position = position_dodge(width = 0.8),
            show.legend = F) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) +
  theme(plot.title = element_text(hjust = 0, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "right") +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(from = 0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  labs(fill = "model") +
  scale_x_discrete(
    name = "Confidence (predicted probability)") +
  scale_fill_manual(values = c("#32a852", "#9d46b8")) +
  scale_color_manual(values = c("#154221", "#3e184a"))

plot_conf

ggsave(filename = paste0("./analysis/figures/cnn_bert_conf_acc_150.png"), plot = plot_conf, width = 11, height = 5)

#==========================================================================================================
# now lets look at the RSA vs transformer confidence
#==========================================================================================================

#loading in RSA data:
sa_data <- read.csv(file = "./output/SASA_scores.csv", header=TRUE, sep=",")

sa_data2 <- sa_data %>%
  select(c(gene, position, SASA_rel_total))
  

trans_data <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")

transf_data <- trans_data %>%
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) %>%
  rename(freq_pred_transf = pred_prob,
         freq_wt_transf = wt_prob,
         aa_wt_transf = aa_wt,
         aa_pred_transf = aa_pred) %>%
  select(gene, position, freq_pred_transf)

for_2D_hist <- inner_join(transf_data, sa_data2)

plot_2D_hist <- for_2D_hist %>%
  ggplot(aes(x = SASA_rel_total, y = freq_pred_transf)) +
  #geom_pointrange() +
  geom_hex(bins = 30) +
  scale_x_continuous(
    name = "Relative Solvent Accessibiity (Å^2)",
    #limits = c(0.0, 1.7),
    #breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Transformer confidence",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.1),
    expand = c(0, 0)) + 
  scale_fill_binned_sequential(palette = "Teal", limits = c(0, 50)) +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_2D_hist

ggsave(filename = paste("./analysis/figures/trans_conf_vs_SA_150.png"), plot = plot_2D_hist, width = 8, height = 6)

joined_data <- inner_join(cnn_data, transf_data)
cnn_data_clean <- joined_data %>%
  select(c(gene, position, freq_pred_cnn))

for_2D_hist <- inner_join(cnn_data_clean, sa_data2)

plot_2D_hist <- for_2D_hist %>%
  ggplot(aes(x = SASA_rel_total, y = freq_pred_cnn)) +
  #geom_pointrange() +
  geom_hex(bins = 30) +
  scale_x_continuous(
    name = "Relative Solvent Accessibiity (Å^2)",
    #limits = c(0.0, 1.7),
    #breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "CNN confidence",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.1),
    expand = c(0, 0)) + 
  scale_fill_binned_sequential(palette = "Teal", limits = c(0, 50)) +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_2D_hist

ggsave(filename = paste("./analysis/figures/cnn_conf_vs_SA_150.png"), plot = plot_2D_hist, width = 8, height = 6)



#==========================================================================================================
# now lets look at n-eff correlations. CNN vs. Transformer.
#==========================================================================================================

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

# lets find n-eff for the transformer data:
trans_data2 <- trans_data %>%
  select(c(gene, position, A:Y)) %>%
  nest(data = c(A:Y)) %>%
  mutate(n_eff_tranf = map(data, get_neff)) %>%
  select(-data)


# preparing alignment data:
align_data2 <- align_data %>%
  select(c(position, gene, n_eff)) %>%
  rename(n_eff_align = n_eff)

# preparing cnn data:
cnn_data_new <- read.csv(file = "./output/stats_cnn.csv", header=TRUE, sep=",")

cnn_data2 <- cnn_data_new %>%
  select(c(gene, position, n_eff)) %>%
  rename(n_eff_cnn = n_eff)

joined_for_cor <- inner_join(trans_data2, cnn_data2)
joined_for_cor2 <- inner_join( joined_for_cor, align_data2)

joined_for_cor2$n_eff_align <- as.numeric(joined_for_cor2$n_eff_align)
joined_for_cor2$n_eff_tranf <- as.numeric(joined_for_cor2$n_eff_tranf)
joined_for_cor2$n_eff_cnn <- as.numeric(joined_for_cor2$n_eff_cnn)

cor_transf <- joined_for_cor2 %>%
  group_by(gene) %>%
  summarise(cor = cor(n_eff_tranf, n_eff_align)) %>%
  mutate(group = "Transformer")

cor_cnn <- joined_for_cor2 %>%
  group_by(gene) %>%
  summarise(cor = cor(n_eff_cnn, n_eff_align)) %>%
  mutate(group = "CNN")

data_cor <- rbind(cor_transf, cor_cnn)

data_cor <- data_cor %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (group == "CNN")),
    dx = rnorm(n(), mean = 0, sd = .05),
    dy = rnorm(n(), mean = 0, sd = .05),
    x_value = as.numeric(factor(group)))  


plot_cor <- data_cor %>%
  ggplot(aes(x = group, y = cor)) +
  geom_path(
    aes(x = as.numeric(factor(group))+dx, y = cor+dy, group = gene, color = color_y),
    size = 0.25) +
  geom_point(
    aes(x = as.numeric(factor(group))+dx, y = cor+dy, group = gene, fill = color_y),
    shape = 21, 
    color = "black",
    size = 2) +
  scale_x_continuous(
    #name = "Neural Network",
    name = "",
    limits = c(0.7,2.3),
    labels = c("CNN \n (structure)", "Transformer \n (sequence)"),
    breaks = (seq(from = 1.0, to = 2.0, by = 1))
  ) +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.3, 0.9),
    breaks = seq(from = 0.0, to = 0.9, by = 0.1),
    expand = c(0, 0)) +
  scale_color_gradient(
    aesthetics = c("color", "fill"), 
    high = "#ffd966", 
    low = "#080845") +
  theme_bw(16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 16),
    panel.grid.minor = element_blank())

plot_cor

ggsave(filename = paste0("./analysis/figures/neff_cor_cnn_bert_150.png"), plot = plot_cor, width = 6, height = 6)


#lets get mean correlations:
means_cor <- data_cor %>%
  select(c(gene, cor, group)) %>%
  group_by(group) %>%
  summarise(mean = mean(cor))
#CNN mean cor: 0.325
#transformer mean cor: 0.404




