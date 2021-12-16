# This script finds the frequencies at with the combined 
# model chooses between the predictions of three models.

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

for_heatplot_with_classes <- for_heatplot_final %>%
  mutate(wt_class = map_chr(aa_wt, calc_class)) %>%
  mutate(class = fct_relevel(wt_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))



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
  if (pred_aa == cnn_pick & pred_aa == bert_pick & pred_aa == esm_pick) {
    return("all models unanimous")
  }
  if (pred_aa != cnn_pick & pred_aa != bert_pick & pred_aa != esm_pick) {
    return("unique")
  }
  return("NA")
}

with_groups <- correct_predictions %>%
  mutate(group = pmap(list(pred_aa, cnn_win_aa, bert_win_aa, esm_win_aa), get_group))
 
  
  
  
  
  
  
  
  
  
  
  


