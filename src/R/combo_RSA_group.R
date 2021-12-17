# plotting RSA vs group

library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)
library(ggridges)

# loading data:
combo_data <- read.csv(file = "./data/transfer_learning_net/predictions__model_combo3_run_1.csv", header=TRUE, sep=",")

#loading in RSA data:
sa_data <- read.csv(file = "./output/SASA_scores.csv", header=TRUE, sep=",")

sa_data2 <- sa_data %>%
  select(c(gene, position, SASA_rel_total))

#categories:
# all three nets unanimous
# only cnn's pick
# only berts pick
# only esm pick
# transformers unanimous
# completely unique

get_group <- function(pred_aa, cnn_pick, bert_pick, esm_pick) {
  
  if (pred_aa == cnn_pick & pred_aa != bert_pick & pred_aa != esm_pick) {
    return("only CNN")
  }
  if (pred_aa == cnn_pick & pred_aa == bert_pick & pred_aa != esm_pick) {
    return("CNN and BERT")
  }
  if (pred_aa != cnn_pick & pred_aa == bert_pick & pred_aa != esm_pick) {
    return("only BERT")
  }
  if (pred_aa != cnn_pick & pred_aa == bert_pick & pred_aa == esm_pick) {
    return("BERT and ESM1b")
  }
  if (pred_aa != cnn_pick & pred_aa != bert_pick & pred_aa == esm_pick) {
    return("only ESM1b")
  }
  if (pred_aa == cnn_pick & pred_aa != bert_pick & pred_aa == esm_pick) {
    return("CNN and ESM1b")
  }
  if (pred_aa == cnn_pick & pred_aa == bert_pick & pred_aa == esm_pick) {
    return("all models unanimous")
  }
  if (pred_aa != cnn_pick & pred_aa != bert_pick & pred_aa != esm_pick) {
    return("unique")
  }
  return("NA")
}

combo_data2 <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa))

# lets look at aa dist of each group:
with_groups <- combo_data2 %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  select(c(gene, position, wt_aa, group))

with_sa <- inner_join(with_groups, sa_data2)

with_sa$group <- as.character(with_sa$group)
with_sa$group <- as.factor(with_sa$group)


fills <- c("#d95a4c", "#e8d146", "#88cfc8", "#e0801f", "#498545", "#8e4e94", "#b3d638", "#0ca0c9")

ridgeline_plot <- with_sa %>%
  ggplot(aes(y = fct_rev(fct_reorder(group, SASA_rel_total)), x = SASA_rel_total, fill = group)) +
  geom_density_ridges(alpha = 0.8) +
  scale_fill_manual(
    values = fills
  ) +
  scale_x_continuous(
    name = "RSA",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0, 0)) + 
  labs(fill = "chosen model") +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))
  

ridgeline_plot

ggsave(filename = "./analysis/figures/combo_groups_RSA.png", plot = ridgeline_plot, width = 10, height = 7)


#----------------------------------------------------------------------------
# only correct predictions:


correct <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa == pred_aa)

# lets look at aa dist of each group:
with_groups <- correct %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  select(c(gene, position, wt_aa, group))

with_sa <- inner_join(with_groups, sa_data2)

with_sa$group <- as.character(with_sa$group)
with_sa$group <- as.factor(with_sa$group)


fills <- c("#d95a4c", "#e8d146", "#88cfc8", "#e0801f", "#498545", "#8e4e94", "#b3d638", "#0ca0c9")

ridgeline_plot_cor <- with_sa %>%
  ggplot(aes(y = fct_rev(fct_reorder(group, SASA_rel_total)), x = SASA_rel_total, fill = group)) +
  geom_density_ridges(alpha = 0.8) +
  scale_fill_manual(
    values = fills
  ) +
  scale_x_continuous(
    name = "RSA",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0, 0)) + 
  labs(fill = "chosen model") +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))


ridgeline_plot_cor

#----------------------------------------------------------------------------
# only mispredictions:


mis <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa != pred_aa)

# lets look at aa dist of each group:
with_groups <- mis %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  select(c(gene, position, wt_aa, group))

with_sa <- inner_join(with_groups, sa_data2)

with_sa$group <- as.character(with_sa$group)
with_sa$group <- as.factor(with_sa$group)


fills <- c("#d95a4c", "#e8d146", "#88cfc8", "#e0801f", "#498545", "#8e4e94", "#b3d638", "#0ca0c9")

ridgeline_plot_mis <- with_sa %>%
  ggplot(aes(y = fct_rev(fct_reorder(group, SASA_rel_total)), x = SASA_rel_total, fill = group)) +
  geom_density_ridges(alpha = 0.8) +
  scale_fill_manual(
    values = fills
  ) +
  scale_x_continuous(
    name = "RSA",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0, 0)) + 
  labs(fill = "chosen model") +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))


ridgeline_plot_mis

#------------------------------------
# histogram of RSA for combo plot:

combo_data2 <- combo_data %>%
  select(c(gene, position, pred_freq))

with_sa <- inner_join(combo_data2, sa_data2)

hist <- with_sa %>%
  ggplot(aes(x = SASA_rel_total, y = pred_freq)) +
  #geom_pointrange() +
  geom_hex(bins = 15) +
  scale_x_continuous(
    name = "RSA (Å^2)",
    limits = c(0.0, 1.0),
    #breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Confidence",
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
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

hist

#-------------------------------------------------------------------------
# get hist for CNN and transformer models:
#-------------------------------------------------------------------------
bert_data <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")
cnn_data <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
esm_data <- read.csv(file = "./output/PSICOV_ESM1b_predictions.csv", header=TRUE, sep=",")
align_data <- read.csv(file = "./output/stats_align_all.csv", header=TRUE, sep=",")
combo_data <- read.csv(file = "./data/transfer_learning_net/predictions__model_combo3_run_1.csv", header=TRUE, sep=",")

# clean cnn data:
cnn_data_clean <- cnn_data %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) %>%
  rename(freq_pred = freq_predicted,
         aa_pred = aa_predicted) %>%
  mutate(group = "CNN")

# clean bert data:

bert_data2 <- bert_data %>% # choose which transformer model to use here
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) %>%
  rename(freq_pred_transf = pred_prob,
         freq_wt_transf = wt_prob,
         aa_wt_transf = aa_wt,
         aa_pred_transf = aa_pred)

bert_data_clean <- bert_data2 %>%
  rename(freq_pred = freq_pred_transf,
         aa_pred = aa_pred_transf,
         freq_wt = freq_wt_transf,
         aa_wt = aa_wt_transf) %>%
  mutate(group = "BERT")

# clean esm1b data:

esm_data2 <- esm_data %>% # choose which transformer model to use here
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) %>%
  rename(freq_pred_transf = pred_prob,
         freq_wt_transf = wt_prob,
         aa_wt_transf = aa_wt,
         aa_pred_transf = aa_pred)

esm_data_clean <- esm_data2 %>%
  rename(freq_pred = freq_pred_transf,
         aa_pred = aa_pred_transf,
         freq_wt = freq_wt_transf,
         aa_wt = aa_wt_transf) %>%
  mutate(group = "ESM1b")

# joining all data:
joined <- rbind(cnn_data_clean, bert_data_clean)
joined2 <- rbind(joined, esm_data_clean)

joined2 <- joined2 %>%
  select(-freq_wt)

combo_data2 <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, pred_freq)) %>%
  rename(aa_wt = wt_aa,
         aa_pred = pred_aa,
         freq_pred = pred_freq) %>%
  mutate(group = "Combined")

joined3 <- rbind(joined2, combo_data2)

gene_combo <- data.frame(unique(combo_data2$gene)) %>%
  rename(gene = unique.combo_data2.gene.)

cnn_hist <- cnn_data_clean %>%
  select(c(gene, position, freq_pred)) %>%
  inner_join(sa_data2) %>%
  inner_join(gene_combo)

bert_hist <- bert_data_clean %>%
  select(c(gene, position, freq_pred)) %>%
  inner_join(sa_data2) %>%
  inner_join(gene_combo)

esm_hist <- esm_data_clean %>%
  select(c(gene, position, freq_pred)) %>%
  inner_join(sa_data2) %>%
  inner_join(gene_combo)

#-----------------------------------------------
#Cnn hist

cnn_hist_plot <- cnn_hist %>%
  ggplot(aes(x = SASA_rel_total, y = freq_pred)) +
  #geom_pointrange() +
  geom_hex(bins = 15) +
  scale_x_continuous(
    name = "RSA (Å^2)",
    limits = c(0.0, 1.0),
    #breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Confidence",
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
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

cnn_hist_plot

bert_hist_plot <- bert_hist %>%
  ggplot(aes(x = SASA_rel_total, y = freq_pred)) +
  #geom_pointrange() +
  geom_hex(bins = 15) +
  scale_x_continuous(
    name = "RSA (Å^2)",
    limits = c(0.0, 1.0),
    #breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Confidence",
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
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

bert_hist_plot

esm_hist_plot <- esm_hist %>%
  ggplot(aes(x = SASA_rel_total, y = freq_pred)) +
  #geom_pointrange() +
  #geom_hex(bins = 30) +
  geom_hex(bins = 15) +
  scale_x_continuous(
    name = "RSA (Å^2)",
    limits = c(0.0, 1.0),
    #breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Confidence",
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
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

esm_hist_plot

# save all data:

abcd <- plot_grid(hist, cnn_hist_plot, bert_hist_plot, esm_hist_plot, ncol = 2, nrow = 2, labels = c('a', 'b', 'c', 'd'))

ggsave(filename = "./analysis/figures/RSA_all_models.png", plot = abcd, width = 11, height = 11)


