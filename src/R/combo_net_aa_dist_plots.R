# this script looks at amino acid and aa class distribution within choice groups

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
    return("small_polar")
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
    return("unique")
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
  return("NA")
}

# lets look at aa dist of each group:
with_groups <- correct_predictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  select(c(gene, position, wt_aa, group, class))

# let's calc frequencies:
freqs <- with_groups %>%
  group_by(group, wt_aa) %>%
  count() %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(sum = sum(n)) %>%
  ungroup() %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  mutate(freq = n/sum)

freqs$group <- as.character(freqs$group)
freqs$group <- as.factor(freqs$group)
freqs$class <- as.character(freqs$class)
freqs$class <- as.factor(freqs$class)

ordered <- freqs %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  mutate(wt_aa = fct_relevel(wt_aa, "G","L","V","A","E","P","I","D","F","T","K","R","Y","S","N","H","Q","W","M","C"))
  
# distribution plot:

fills <- c("#990008", "#0a2575", "#b35900", "#1a6600", "#5c0679", "#9e9e2e")

plot_correct <- ordered %>%
  ggplot(aes(y = freq, x = wt_aa, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(group), nrow = 4, ncol = 2) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")
    ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0.0, 0.30),
    breaks = seq(0.0, 0.30, by = 0.1),
    expand = c(0, 0)) + 
  scale_x_discrete(
    name = "Wild type amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_correct

#ggsave(filename = "./analysis/figures/cor_groups_aa_dist.png", plot = plot_correct, width = 10, height = 9)


#----------------------------------------------------------------------------------------
# now I need a plot that shows the distribution of groups for each amino acid type

#plot with all sites:


combo_data2 <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) 


# lets look at group dist of each aa:
with_groups <- combo_data2 %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  select(c(gene, position, wt_aa, group, class))

# let's calc frequencies:
freqs <- with_groups %>%
  group_by(group, wt_aa) %>%
  count() %>%
  ungroup() %>%
  group_by(wt_aa) %>%
  mutate(sum = sum(n)) %>%
  ungroup() %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  mutate(freq = n/sum)

freqs$group <- as.character(freqs$group)
freqs$group <- as.factor(freqs$group)
freqs$class <- as.character(freqs$class)
freqs$class <- as.factor(freqs$class)

ordered <- freqs %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  mutate(wt_aa = fct_rev(fct_relevel(wt_aa, "G","P","L","F","V","W","D","A","Y","I","E","T","C","S","H","R","K","N","M","Q"))) %>%
  mutate(group = fct_rev(fct_relevel(group, "all models unanimous", "BERT and ESM1b", "only CNN", "CNN and BERT", "only BERT", "CNN and ESM1b", "only ESM1b", "unique")))


# distribution plot:

fills <- c("#d95a4c", "#e8d146", "#88cfc8", "#e0801f", "#498545", "#8e4e94", "#b3d638", "#0ca0c9")

plot_all <- ordered %>%
  ggplot(aes(y = freq, x = fct_rev(wt_aa), fill = fct_rev(group))) +
  geom_col(alpha = 1.0, position = "stack") +
  #facet_wrap(vars(group), nrow = 4, ncol = 2) +
  scale_fill_manual(
    values = fills
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) + 
  scale_x_discrete(
    name = "Wild type amino acid",
    expand = c(0.0, 0.0)) + 
  theme_cowplot(16) +
  labs(fill = "chosen model") +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_all

#--------------------------------------------
# lets make a within class plot:

#plot with correct sites:


correct_predictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa == pred_aa)

# lets look at group dist of each aa:
with_groups <- correct_predictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  select(c(gene, position, wt_aa, group, class))

# let's calc frequencies:
freqs <- with_groups %>%
  group_by(group, wt_aa) %>%
  count() %>%
  ungroup() %>%
  group_by(wt_aa) %>%
  mutate(sum = sum(n)) %>%
  ungroup() %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  mutate(freq = n/sum)

freqs$group <- as.character(freqs$group)
freqs$group <- as.factor(freqs$group)
freqs$class <- as.character(freqs$class)
freqs$class <- as.factor(freqs$class)

ordered <- freqs %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  mutate(wt_aa = fct_rev(fct_relevel(wt_aa, "G","P","L","F","V","W","D","A","Y","I","E","T","C","S","H","R","K","N","M","Q"))) %>%
  mutate(group = fct_rev(fct_relevel(group, "all models unanimous", "BERT and ESM1b", "only CNN", "CNN and BERT", "only BERT", "CNN and ESM1b", "only ESM1b", "unique")))


# distribution plot:

fills <- c("#d95a4c", "#e8d146", "#88cfc8", "#e0801f", "#498545", "#8e4e94", "#b3d638", "#0ca0c9")

plot_cor_pred <- ordered %>%
  ggplot(aes(y = freq, x = fct_rev(wt_aa), fill = fct_rev(group))) +
  geom_col(alpha = 1.0, position = "stack") +
  #facet_wrap(vars(group), nrow = 4, ncol = 2) +
  scale_fill_manual(
    values = fills
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) + 
  scale_x_discrete(
    name = "Wild type amino acid",
    expand = c(0.0, 0.0)) + 
  theme_cowplot(16) +
  labs(fill = "chosen model") +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_cor_pred

#---------------------------------------
#plot with mispredicted sites:


mispredictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa != pred_aa)

# lets look at group dist of each aa:
with_groups <- mispredictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  select(c(gene, position, wt_aa, group, class))

# let's calc frequencies:
freqs <- with_groups %>%
  group_by(group, wt_aa) %>%
  count() %>%
  ungroup() %>%
  group_by(wt_aa) %>%
  mutate(sum = sum(n)) %>%
  ungroup() %>%
  mutate(class = map_chr(wt_aa, calc_class)) %>%
  mutate(freq = n/sum)

freqs$group <- as.character(freqs$group)
freqs$group <- as.factor(freqs$group)
freqs$class <- as.character(freqs$class)
freqs$class <- as.factor(freqs$class)

ordered <- freqs %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  mutate(wt_aa = fct_rev(fct_relevel(wt_aa, "G","P","L","F","V","W","D","A","Y","I","E","T","C","S","H","R","K","N","M","Q"))) %>%
  mutate(group = fct_rev(fct_relevel(group, "all models unanimous", "BERT and ESM1b", "only CNN", "CNN and BERT", "only BERT", "CNN and ESM1b", "only ESM1b", "unique")))


# distribution plot:

fills <- c("#d95a4c", "#e8d146", "#88cfc8", "#e0801f", "#498545", "#8e4e94", "#b3d638", "#0ca0c9")


plot_mispred <- ordered %>%
  ggplot(aes(y = freq, x = fct_rev(wt_aa), fill = fct_rev(group))) +
  geom_col(alpha = 1.0, position = "stack") +
  #facet_wrap(vars(group), nrow = 4, ncol = 2) +
  #scale_fill_brewer(palette = "Paired") +
  scale_fill_manual(
    values = fills
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) + 
  scale_x_discrete(
    name = "Wild type amino acid",
    expand = c(0.0, 0.0)) + 
  theme_cowplot(16) +
  labs(fill = "chosen model") +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_mispred

#----------------------------------------------------------------------------
#saving plots:

abc <- plot_grid(plot_all + theme(legend.position="none"), 
                    plot_cor_pred + theme(legend.position="none"), 
                    plot_mispred + theme(legend.position="none"), 
                    nrow = 3, 
                    ncol = 1, 
                    align = "hv", 
                    axis = "bl",
                    labels = c('a', 'b', 'c'))

legend <- get_legend(
  # create some space to the left of the legend
  plot_all + theme(legend.box.margin = margin(0, 0, 0, 12)))

with_legend <- plot_grid(abc, legend, rel_widths = c(3, 1))
with_legend

#ggsave(filename = "./analysis/figures/group_dist_aa.png", plot = with_legend, width = 11, height = 10)












