# This script finds the frequencies at with the combined 
# model chooses between the predictions of three models.

library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)


# loading data:
combo_data <- read_csv(file = "./data/transfer_learning_net/predictions__model_combo3_run_1.csv")

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

# lets plot freq of each group:

with_groups <- correct_predictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group))

with_groups$group <- unlist(with_groups$group)

counts <- with_groups %>%
  count(group)

group_freqs_cor <- counts %>%
  mutate(sum = sum(n)) %>%
  mutate(freq = n/sum) %>%
  rename(count = n) %>%
  select(c(group, count, freq)) 
  #mutate(label_text = paste0(count, " (", round(freq*100, 1), "%)"))

group_freqs_cor$group <- as.character(group_freqs_cor$group)

groups_only <- with_groups %>%
  select(group)

cor <- group_freqs_cor
# correct sites bar plot:
cor_bar_plot <- group_freqs_cor %>%
  ggplot(aes(y = fct_reorder(group, freq), x = freq)) +
  geom_col(fill = "#839ea8") +
  scale_y_discrete(
    name = "",
    expand = c(0,0),
    ) +
  scale_x_continuous(
    name = "Proportion of sites",
    limits = c(0.0, 0.56),
    breaks = seq(from = 0.0, to = 0.5, by = 0.1),
    expand = expansion(add = c(0, 0.05))
   ) +
  geom_text(aes(label = count), 
            hjust = - 0.1,
            color = "black",
            size = 3) +
  theme_cowplot(10) +
  theme(
    panel.grid.major.x = element_line(colour="gray90", size=0.1),
    axis.text = element_text(
      color = "black", 
      size = 10)) 
  #theme(plot.margin = unit(c(0,2,0,0), "lines")) #top, right, bottom, left

cor_bar_plot
  
#--------------------------------------------------
# lets make the same plot for mispredictions:

mispredictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa != pred_aa)

# lets plot freq of each group:

with_groups <- mispredictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group))

counts <- with_groups %>%
  select(group) %>%
  count(group)

group_freqs_mis <- counts %>%
  mutate(sum = sum(n)) %>%
  mutate(freq = n/sum) %>%
  rename(count = n) %>%
  select(c(group, count, freq)) 
  #mutate(label_text = paste0(count, " (", round(freq*100, 1), "%)"))

group_freqs_mis$group <- as.character(group_freqs_mis$group)

groups_only <- with_groups %>%
  select(group)

# mispredicted sites bar plot:

order <- c("all models unanimous", "BERT and ESM1b", "only CNN", "CNN and BERT", "only BERT", "CNN and ESM1b", "only ESM1b", "unique")

mis <- group_freqs_mis
mis_bar_plot <- group_freqs_mis %>%
  ggplot(aes(y = fct_rev(fct_relevel(group, order)), x = freq)) +
  geom_col(fill = "#839ea8") +
  scale_y_discrete(
    name = "",
    expand = c(0,0),
  ) +
  scale_x_continuous(
    name = "Proportion of sites",
    limits = c(0.0, 0.56),
    breaks = seq(from = 0.0, to = 0.5, by = 0.1),
    expand = expansion(add = c(0, 0.05))
  ) +
  geom_text(aes(label = count), 
            hjust = - 0.1,
            color = "black",
            size = 3) +
  theme_cowplot(10) +
  theme(
    panel.grid.major.x = element_line(colour="gray90", size=0.1),
    axis.text = element_text(
      color = "black", 
      size = 10)) 
mis_bar_plot


#--------------------------------------------------
# lets make the same plot for all predictions:

all_predictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa))

# lets plot freq of each group:

with_groups <- all_predictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group))

counts <- with_groups %>%
  select(group) %>%
  count(group)

group_freqs <- counts %>%
  mutate(sum = sum(n)) %>%
  mutate(freq = n/sum) %>%
  rename(count = n) %>%
  select(c(group, count, freq)) 
#mutate(label_text = paste0(count, " (", round(freq*100, 1), "%)"))

group_freqs$group <- as.character(group_freqs$group)

groups_only <- with_groups %>%
  select(group)

# all sites bar plot:

order <- c("all models unanimous", "BERT and ESM1b", "only CNN", "only BERT", "CNN and BERT", "CNN and ESM1b", "only ESM1b", "unique")
all <- group_freqs
all_bar_plot <- group_freqs %>%
  ggplot(aes(y = fct_rev(fct_relevel(group, order)), x = freq)) +
  geom_col(fill = "#839ea8") +
  scale_y_discrete(
    name = "",
    expand = c(0,0),
  ) +
  scale_x_continuous(
    name = "Proportion of sites",
    limits = c(0.0, 0.56),
    breaks = seq(from = 0.0, to = 0.5, by = 0.1),
    expand = expansion(add = c(0, 0.05))
  ) +
  geom_text(aes(label = count), 
            hjust = - 0.1,
            color = "black",
            size = 3) +
  theme_cowplot(10) +
  theme(
    panel.grid.major.x = element_line(colour="gray90", size=0.1),
    axis.text = element_text(
      color = "black", 
      size = 10)) +
  theme(plot.margin = unit(c(0,2,0,0), "lines")) #top, right, bottom, left

all_bar_plot

#--------------------------------------------------
# lets make the filled plot (proportions):

all_predictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  mutate(prediction = ifelse(wt_aa == pred_aa, "correct", "misprediction"))

# lets plot freq of each group:

with_groups <- all_predictions %>%
  mutate(group = pmap_chr(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  select(prediction, group)

filled_bar_plot <- with_groups %>%
  ggplot(aes(y = fct_rev(fct_relevel(group, "all models unanimous", "CNN and BERT", "CNN and ESM1b", "BERT and ESM1b", "only CNN", "only BERT", "only ESM1b", "unique")), fill = fct_rev(prediction))) +
  geom_bar(position = "fill") +
  scale_y_discrete(
    name = "",
    expand = c(0,0),
  ) +
  scale_x_continuous(
    name = "Proportion of sites",
    limits = c(0.0, 0.56),
    breaks = seq(from = 0.0, to = 0.5, by = 0.1),
    expand = expansion(add = c(0, 0.05))
  ) +
  geom_text(aes(label = count),
            hjust = - 0.1,
            color = "black",
            size = 3) +
  labs(fill = "prediction") +
  theme_cowplot(10) +
  theme(
    panel.grid.major.x = element_line(colour="gray90", size=0.1),
    axis.text = element_text(
      color = "black", 
      size = 10)) 
  #theme(plot.margin = unit(c(0,2,0,0), "lines")) #top, right, bottom, left

filled_bar_plot



s#--------------------------------------------------------------------------------------
# lets make a faceted plot:

all <- all %>%
  mutate(data = "all data")

mis <- mis %>%
  mutate(data = "mispredictions")

cor <- cor %>%
  mutate(data = "correct \npredictions")

all_mis_cor <- rbind(all, cor, mis)

facet_plot <- all_mis_cor %>%
  ggplot(aes(y = fct_rev(fct_relevel(group, order)), x = freq)) +
  geom_col(fill = "#66102a",
           color = "#3b0414",
           alpha = 0.6,
           size = 0.3
  ) +
  facet_wrap(~data) +
  scale_y_discrete(
    name = "",
    expand = c(0,0),
  ) +
  scale_x_continuous(
    name = "Proportion of sites",
    limits = c(0.0, 0.56),
    breaks = seq(from = 0.0, to = 0.5, by = 0.1),
    expand = expansion(add = c(0, 0.05))
  ) +
  geom_text(aes(label = count), 
            hjust = - 0.1,
            color = "black",
            size = 3) +
  theme_cowplot(10) +
  theme(
    panel.grid.major.x = element_line(colour="gray90", size=0.1),
    axis.text = element_text(
      color = "black", 
      size = 10)) 

facet_plot
#ggsave(filename = "./analysis/figures/choices_barplot_facet.png", plot = facet_plot, width = 9, height = 5)



#Lets get odds ratio for these plots

odds <- group_freqs_mis %>%
  rename(freq_mis = freq, 
         count_mis = count) %>%
  inner_join(group_freqs_cor) %>%
  rename(freq_cor = freq,
         count_cor = count) %>%
  mutate(count_ratio = count_cor/count_mis,
         freq_ratio = freq_cor/freq_mis)
  
odds_bar_plot <- odds %>%
  ggplot(aes(y = fct_reorder(group, count_ratio), x = count_ratio)) +
  geom_col(fill = "#839ea8") +
  scale_y_discrete(
    name = "",
    expand = c(0,0),
  ) +
  scale_x_log10(
    name = "Correct count over \nmispredicted count",
    #limits = c(0, 3.5),
    #breaks = seq(0, 3, by = 0.5),
    expand = c(0.01, 0.1)) + 
  geom_text(aes(label = round(freq_ratio, 2)),
            hjust = -0.1,
            color = "black",
            size = 3) +
  theme_cowplot(10) +
  theme(
    panel.grid.major.x = element_line(colour="gray90", size=0.1),
    axis.text = element_text(
      color = "black", 
      size = 10))

odds_bar_plot 

cor_and_mis_plots <- plot_grid(cor_bar_plot, mis_bar_plot, all_bar_plot, ncol = 3, nrow = 1, labels = c('a', 'b', 'c'))
#with_odds <- plot_grid(cor_and_mis_plots, odds_bar_plot, nrow = 2, ncol = 1, labels = c('', 'c'), scale = 0.9)
#with_odds
cor_and_mis_plots

ggsave(filename = "./analysis/figures/choices_barplot_cor_mis_all.png", plot = cor_and_mis_plots, width = 9, height = 5)


#-------------------------------------------------------------------------------
# Now I want to know what percent in each of the 8 groups is correct.
# Question: which scenario gives the highest proportion of correct predictions? 
# and which groups have the lowest?

combo_data2 <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa))
 
with_groups <- combo_data2 %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  mutate(status = ifelse(pred_aa == wt_aa, "correct", "mispredicted"))

counts <- with_groups %>%
  group_by(group) %>%
  count(status) %>%
  pivot_wider(names_from = status, values_from = n) %>%
  mutate(sum = correct + mispredicted) %>%
  mutate(freq_cor = correct/sum) %>%
  select(c(group, freq_cor)) %>%
  mutate(label_text = paste0(round(freq_cor*100, 1), "%"))


counts$group <- unlist(counts$group)
counts <- as.data.frame(counts) 
counts$group <- as.factor(counts$group)
counts$label_text <- as.character(counts$label_text)

# within group freq bar plot:

#order <- c("all models unanimous", "BERT and ESM1b", "only CNN", "CNN and BERT", "only BERT", "CNN and ESM1b", "only ESM1b", "unique")


within_bar_plot <- counts %>%
  ggplot(aes(y = fct_reorder(group, freq_cor), x = freq_cor)) +
  geom_col(fill = "#8ec8cb", color = "#316785", width = 0.8, size = 0.3) +
  scale_y_discrete(
    name = "",
    expand = c(0,0)
  ) +
  scale_x_continuous(
    name = "Proportion correct",
    expand = c(0.0,0.0),
    limits = c(0.0, 1.05),
    breaks = (seq(from = 0.0, to = 1.0, by = 0.2))
  ) +
  geom_text(aes(label = label_text),
            #vjust = -0.5,
            hjust = -0.1,
            color = "black") +
  theme_cowplot(11) +
  theme(
    panel.background = element_blank(),
    #panel.grid.major.x = element_line(colour="gray80", size=0.1),
    #panel.grid.minor.x = element_line(colour="gray90", size=0.1),
    axis.text = element_text(
      color = "black", 
      size = 11))

within_bar_plot 

#ggsave(filename = "./analysis/figures/choices_barplot_within_group3.png", plot = within_bar_plot, width = 8.5, height = 5)

with_groups$group <- unlist(with_groups$group)
with_groups <- as.data.frame(with_groups) 
with_groups$group <- as.factor(with_groups$group)
with_groups$label_text <- as.character(with_groups$label_text)

# get bar plot with proportion of sites:
within_bar_plot2 <- with_groups %>%
  ggplot(aes(y = group, fill = status)) +
  geom_bar(width = 0.8, 
           size = 0.3, 
           position = "stack"
           ) +
  # scale_y_discrete(
  #   name = "",
  #   expand = c(0,0)
  # ) +
  # scale_x_continuous(
  #   name = "Proportion of sites",
  #   expand = c(0.0,0.0),
  #   limits = c(0.0, 1.05),
  #   breaks = (seq(from = 0.0, to = 1.0, by = 0.2))
  # ) +
  # geom_text(aes(label = label_text),
  #           #vjust = -0.5,
  #           hjust = -0.1,
  #           color = "black") +
  theme_cowplot(11) +
  theme(
    panel.background = element_blank(),
    #panel.grid.major.x = element_line(colour="gray80", size=0.1),
    #panel.grid.minor.x = element_line(colour="gray90", size=0.1),
    axis.text = element_text(
      color = "black", 
      size = 11))

within_bar_plot2



#-----------------------------------------------------------------------------------------------
# now find the percent of positions where the highest conf aa was chosen. What is that percent?
#-----------------------------------------------------------------------------------------------
  

#-----------------------------------------------------------------------
#finding unique and mispredicted positions (MISSED OPPORTUNITIES)
mispredictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa != pred_aa) 

unique <- mispredictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group)) %>%
  filter(group == "unique")

unique2 <- unique %>%
  mutate(missed = ifelse(wt_aa == cnn_win_aa | wt_aa == bert_win_aa | wt_aa == esm_win_aa, "missed_cor", "no_cor"))

freq_unique <- unique2 %>%
  group_by(missed) %>%
  count()

# out of the 55 positions that were unique
# in 38.2% of positions the correct amino acid was one of the predicted by a model. 



# look at all mispredicted positions:
mispredictions <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) %>%
  filter(wt_aa != pred_aa) 

mis <- mispredictions %>%
  mutate(missed = ifelse(wt_aa == cnn_win_aa | wt_aa == bert_win_aa | wt_aa == esm_win_aa, "missed_cor", "no_cor"))

freq_mis <- mis %>%
  group_by(missed) %>%
  count() %>%
  ungroup() %>%
  mutate(sum = sum(n)) %>%
  mutate(freq = n/sum)
# missed the correct aa in 36% of mispredictions. 

#-------------------------------------------------------------------------
# in the positions where the answer was not there, what is the % correct?
combo_data2 <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa)) 

not_there <- combo_data2 %>%
  mutate(missed = ifelse(wt_aa == cnn_win_aa | wt_aa == bert_win_aa | wt_aa == esm_win_aa, "missed_cor", "no_cor")) %>%
  mutate(status = ifelse(wt_aa == pred_aa, "correct", "mispredicted"))

freq_not_there <- not_there %>%
  group_by(missed, status) %>%
  count()

# only 28 out of 548 positions where the answer was not there were correct (5%)  

# out of all mispredicted sites, (814) only 64% of the time, the answer was not predicted by any model. 


# next, talk about the unique category. 

