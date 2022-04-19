# aa class dist plots


library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)
divergingx_palettes(n = 7, plot = TRUE)

# loading data:
combo_data <- read.csv(file = "./data/transfer_learning_net/predictions__model_combo3_run_1.csv", header=TRUE, sep=",")

correct_predictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa == pred_aa)

#add aa class label: 
calc_class <- function(x) {
  aliphatic = c("M", "L", "I", "V", "A")
  small_polar = c("C", "S", "T", "N", "Q")
  negative = c("D", "E")
  positive = c("R", "K")
  aromatic = c("H", "Y", "W", "F")
  unique = c("P", "G")
  
  if (x %in% aliphatic) {
    return("aliphatic")
  }
  if (x %in% small_polar) {
    return("small polar")
  }
  if (x %in% negative) {
    return("negative")
  }
  if (x %in% positive) {
    return("positive")
  }
  if (x %in% aromatic) {
    return("aromatic")
  }
  if (x %in% unique) {
    return("G or P")
  }
  return("not found")
}


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
  return(NA_character_)
}



# lets look at aa dist of each group:
with_groups <- correct_predictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  select(c(gene, position, wt_aa, group, class))

for_cor_plot <- with_groups %>%
  count(class, group) 
#%>%
  group_by(class) %>%
  mutate(sum = sum(n)) %>%
  mutate(freq = n/sum) 


for_cor_plot$group <- as.character(for_cor_plot$group)
for_cor_plot$group <- as.factor(for_cor_plot$group)
for_cor_plot$class <- as.character(for_cor_plot$class)
for_cor_plot$class <- as.factor(for_cor_plot$class)

for_cor_plot <- for_cor_plot %>%
  mutate(class = fct_relevel(class, "aliphatic", "small polar", "negative", "positive", "aromatic", "G or P")) %>%
  mutate(group = fct_relevel(group, "all models unanimous", "only CNN", "CNN and BERT", "BERT and ESM1b", "CNN and ESM1b", "only BERT", "unique", "only ESM1b"))

fills <- c("#d95a4c", "#e8d146", "#88cfc8", "#e0801f", "#498545", "#8e4e94", "#b3d638", "#0ca0c9")

cor_plot <- for_cor_plot %>%
  ggplot(aes(x = freq, y = fct_rev(group), fill = group)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(class), nrow = 1, ncol = 6) +
  scale_fill_manual(
    values = fills
  ) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.8),
    breaks = seq(0.0, 0.8, by = 0.4),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

cor_plot

#----------------------------------------------------------------------
# misprediction plot
#----------------------------------------------------------------------

mispredictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa != pred_aa)

# lets look at aa dist of each group:
with_groups <- mispredictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  select(c(gene, position, wt_aa, group, class))

for_mis_plot <- with_groups %>%
  group_by(class, group) %>%
  count() %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(sum = sum(n)) %>%
  mutate(freq = n/sum)

for_mis_plot$group <- as.character(for_mis_plot$group)
for_mis_plot$group <- as.factor(for_mis_plot$group)
for_mis_plot$class <- as.character(for_mis_plot$class)
for_mis_plot$class <- as.factor(for_mis_plot$class)

for_mis_plot <- for_mis_plot %>%
  mutate(class = fct_relevel(class, "aliphatic", "small polar", "negative", "positive", "aromatic", "G or P")) %>%
  mutate(group = fct_relevel(group, "all models unanimous", "only CNN", "CNN and BERT", "BERT and ESM1b", "CNN and ESM1b", "only BERT", "unique", "only ESM1b"))

fills <- c("#d95a4c", "#e8d146", "#88cfc8", "#e0801f", "#498545", "#8e4e94", "#b3d638", "#0ca0c9")

mis_plot <- for_mis_plot %>%
  ggplot(aes(x = freq, y = fct_rev(group), fill = group)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(class), nrow = 1, ncol = 6) +
  scale_fill_manual(
    values = fills
  ) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.8),
    breaks = seq(0.0, 0.8, by = 0.4),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

mis_plot

#-----------------------------------------------------------
#odds ratio plot

odds <- for_mis_plot %>%
  rename(freq_mis = freq) %>% 
  select(c(class, group, freq_mis)) %>%
  inner_join(for_cor_plot) %>%
  rename(freq_cor = freq) %>%
  select(c(class, group, freq_mis, freq_cor)) %>%
  mutate(freq_ratio = freq_cor/freq_mis)

odds_plot <- odds %>%
  ggplot(aes(x = freq_ratio, y = fct_rev(group), fill = group)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(class), nrow = 1, ncol = 6) +
  scale_fill_manual(
    values = fills
  ) +
  scale_x_log10(
    name = "Frequency correct over \n frequency mispredicted",
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

odds_plot

abc <- plot_grid(cor_plot, mis_plot, odds_plot, ncol = 1, nrow = 3, labels = c('a', 'b', 'c'))

#ggsave(filename = "./analysis/figures/group_aa_class_dist.png", plot = abc, width = 13, height = 8)

#---------------------------------
# redo of odds plot:

odds_plot <- odds %>%
  ggplot(aes(x = freq_ratio, y = fct_rev(group), fill = group)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(class), nrow = 3, ncol = 2) +
  scale_fill_manual(
    values = fills
  ) +
  scale_x_log10(
    name = "Frequency correct over \n frequency mispredicted",
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

odds_plot

#ggsave(filename = "./analysis/figures/odds_plot_classes_groups.png", plot = odds_plot, width = 12, height = 8)

#-------------------------------------------------------------------------------------
# aa dist with new groups
#-------------------------------------------------------------------------------------


#categories:
# all three nets unanimous
# only cnn's pick
# only berts pick
# only esm pick
# transformers unanimous
# completely unique


get_group2 <- function(pred_aa, cnn_pick, bert_pick, esm_pick) {
  
  # CNN aa matches pred aa (only CNN or also one of the transformers)
  if (pred_aa == cnn_pick & pred_aa != bert_pick & pred_aa != esm_pick |
      pred_aa == cnn_pick & pred_aa == bert_pick & pred_aa != esm_pick |
      pred_aa == cnn_pick & pred_aa != bert_pick & pred_aa == esm_pick) {
    return("prediction matches \nCNN")
  }
  # a transformer (or both) match the prediction (CNN does not match)
  if (pred_aa == bert_pick & pred_aa != esm_pick & pred_aa != cnn_pick |
      pred_aa != bert_pick & pred_aa == esm_pick & pred_aa != cnn_pick|
      pred_aa == bert_pick & pred_aa == esm_pick & pred_aa != cnn_pick){
    return("prediction matches \ntransformer/s")
  }
  # unique prediction (non match)
  if (pred_aa != cnn_pick & pred_aa != bert_pick & pred_aa != esm_pick) {
    return("unique prediction")
  }
  # all models unanimous
  if (pred_aa == cnn_pick & pred_aa == bert_pick & pred_aa == esm_pick) {
    return("prediction matches \nall models")
  }
  return("NA")
}

# lets look at aa dist of each group:
with_groups <- correct_predictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group2)) %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  select(c(gene, position, wt_aa, group, class))

for_cor_plot <- with_groups %>%
  group_by(wt_aa, group) %>%
  count() %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(sum = sum(n)) %>%
  mutate(freq = n/sum)

for_cor_plot$group <- as.character(for_cor_plot$group)
for_cor_plot$group <- as.factor(for_cor_plot$group)
for_cor_plot$wt_aa <- as.character(for_cor_plot$wt_aa)
for_cor_plot$wt_aa <- as.factor(for_cor_plot$wt_aa)

for_cor_plot <- for_cor_plot %>%
  mutate(class = map_chr(wt_aa, calc_class))

for_cor_plot$class <- as.character(for_cor_plot$class)
for_cor_plot$class <- as.factor(for_cor_plot$class)
  

for_cor_plot <- for_cor_plot %>%
  mutate(class = fct_relevel(class, "aliphatic", "small polar", "negative", "positive", "aromatic", "G or P")) %>%
  mutate(group = fct_relevel(group, "prediction matches \nCNN", "prediction matches \ntransformer/s", "prediction matches \nall models", "unique prediction")) %>%
  mutate(wt_aa = fct_relevel(wt_aa, "I","V","A","L","G","E","T","D","F","K","P","R","S","N","Y","M","H","Q","C","W"))


fills <- c("#990008", "#0a2575", "#b35900", "#1a6600", "#5c0679", "#9e9e2e")

cor_plot <- for_cor_plot %>%
  ggplot(aes(y = freq, x = wt_aa, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(group)) +
  scale_fill_manual(
    values = fills
  ) +
  scale_y_continuous(
    name = "Frequency",
    # limits = c(0.0, 0.8),
    # breaks = seq(0.0, 0.8, by = 0.4),
    expand = expansion(add = c(0,0.005))) + 
  scale_x_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

cor_plot

#ggsave(filename = "./analysis/figures/four_groups_aa_dist2.png", plot = cor_plot, width = 11, height = 8.5)

#--------------------------------------------------------------------------------------
#lets look at the distribution of secondary structure for the four groups above.

sec_struc <- read.csv(file = "./output/second_struc.csv", header=TRUE, sep=",")
















