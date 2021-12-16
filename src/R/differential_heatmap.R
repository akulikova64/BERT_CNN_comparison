library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)
divergingx_palettes(n = 7, plot = TRUE)

#differential heat maps

# loading data
bert_data <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")
cnn_data <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
esm_data <- read.csv(file = "./output/PSICOV_ESM1b_predictions.csv", header=TRUE, sep=",")
combo_data <- read.csv(file = "./data/transfer_learning_net/predictions__model_combo3_run_1.csv", header=TRUE, sep=",")

genes <- data.frame(unique(combo_data$gene)) %>%
  rename(gene = unique.combo_data.gene.)


# clean cnn data:
cnn_data_clean <- cnn_data %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) %>%
  rename(freq_pred = freq_predicted,
         aa_pred = aa_predicted) %>%
  mutate(group = "CNN") %>%
  select(-freq_wt)

cnn_data_clean <- inner_join(cnn_data_clean, genes)

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
  mutate(group = "BERT") %>%
  select(-freq_wt)

bert_data_clean <- inner_join(bert_data_clean, genes)

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
  mutate(group = "ESM1b") %>%
  select(-freq_wt)

esm_data_clean <- inner_join(esm_data_clean, genes)

combo_data_clean <- combo_data %>%
  select(c(gene, position, wt_aa, pred_aa, pred_freq)) %>%
  rename(aa_wt = wt_aa,
         aa_pred = pred_aa,
         freq_pred = pred_freq) %>%
  mutate(group = "Combined")

# getting counts of wt-predicted amino acid pairs
cnn <- cnn_data_clean %>%
  select(c(group, aa_wt, aa_pred)) %>%
  group_by(aa_wt, aa_pred) %>%
  summarize(count = n())

bert <- bert_data_clean %>%
  select(c(group, aa_wt, aa_pred)) %>%
  group_by(aa_wt, aa_pred) %>%
  summarize(count = n())

esm <- esm_data_clean %>%
  select(c(group, aa_wt, aa_pred)) %>%
  group_by(aa_wt, aa_pred) %>%
  summarize(count = n())

combo <- combo_data_clean %>%
  select(c(group, aa_wt, aa_pred)) %>%
  group_by(aa_wt, aa_pred) %>%
  summarize(count = n())

# calculating the frequency:

cnn <- cnn %>%
  group_by(aa_wt) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  rename(cnn_freq = freq) %>%
  select(-c(count))

bert <- bert %>%
  group_by(aa_wt) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  rename(bert_freq = freq) %>%
  select(-c(count))

esm <- esm %>%
  group_by(aa_wt) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  rename(esm_freq = freq) %>%
  select(-c(count))

combo <- combo %>%
  group_by(aa_wt) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  rename(combo_freq = freq) %>%
  select(-c(count))


cnn_joined <- full_join(cnn, combo)

cnn_joined <- cnn_joined %>%
  mutate(cnn_freq = ifelse(is.na(cnn_freq), 0, cnn_freq),
         combo_freq = ifelse(is.na(combo_freq), 0, combo_freq)) %>%
  mutate(freq = combo_freq - cnn_freq)

for_heatplot_final <- cnn_joined %>%
  mutate(
    predicted = fct_rev(fct_relevel(aa_pred, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")),
    wt = fct_relevel(aa_wt, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")) 

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

font_size <- 10

plot_cnn_aa <- ggplot() +
  geom_tile(data = for_heatplot_with_classes, aes(x = wt, y = predicted, fill = freq)) + 
  scale_fill_continuous_diverging(
    palette = "Blue-Red 3", 
    #palette = "Tropic",
    rev = TRUE, 
    p2 = 1
  ) +
  scale_x_discrete(
    name = "Wild type residue",
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Predicted residue",
    expand = c(0,0)) +
  labs(fill = "Change in \nfrequency") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = font_size),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_cnn_aa

#----------------------------------------
# now bert heat map:

bert_joined <- full_join(bert, combo)

bert_joined <- bert_joined %>%
  mutate(bert_freq = ifelse(is.na(bert_freq), 0, bert_freq),
         combo_freq = ifelse(is.na(combo_freq), 0, combo_freq)) %>%
  mutate(freq = combo_freq - bert_freq)

for_heatplot_final <- bert_joined %>%
  mutate(
    predicted = fct_rev(fct_relevel(aa_pred, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")),
    wt = fct_relevel(aa_wt, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")) 

for_heatplot_with_classes <- for_heatplot_final %>%
  mutate(wt_class = map_chr(aa_wt, calc_class)) %>%
  mutate(class = fct_relevel(wt_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

font_size <- 10

plot_bert_aa <- ggplot() +
  geom_tile(data = for_heatplot_with_classes, aes(x = wt, y = predicted, fill = freq)) + 
  scale_fill_continuous_diverging(
    palette = "Blue-Red 3", 
    #palette = "Tropic",
    rev = TRUE, 
    p2 = 1
  ) +
  scale_x_discrete(
    name = "Wild type residue",
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Predicted residue",
    expand = c(0,0)) +
  labs(fill = "Change in frequency \n from BERT model") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = font_size),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_bert_aa

#----------------------------------------
# now esm heat map:

esm_joined <- full_join(esm, combo)

esm_joined <- esm_joined %>%
  mutate(esm_freq = ifelse(is.na(esm_freq), 0, esm_freq),
         combo_freq = ifelse(is.na(combo_freq), 0, combo_freq)) %>%
  mutate(freq = combo_freq - esm_freq)

for_heatplot_final <- esm_joined %>%
  mutate(
    predicted = fct_rev(fct_relevel(aa_pred, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")),
    wt = fct_relevel(aa_wt, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")) 

for_heatplot_with_classes <- for_heatplot_final %>%
  mutate(wt_class = map_chr(aa_wt, calc_class)) %>%
  mutate(class = fct_relevel(wt_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

font_size <- 10

plot_esm_aa <- ggplot() +
  geom_tile(data = for_heatplot_with_classes, aes(x = wt, y = predicted, fill = freq)) + 
  scale_fill_continuous_diverging(
    palette = "Blue-Red 3", 
    #palette = "Tropic",
    rev = TRUE, 
    p2 = 1
  ) +
  scale_x_discrete(
    name = "Wild type residue",
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Predicted residue",
    expand = c(0,0)) +
  labs(fill = "Change in frequency \n from ESM1b model") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = font_size),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_esm_aa

#-------------------------------------------------------------------------------
# now I need a heatmap of amino acid classes
#-------------------------------------------------------------------------------

# getting counts of wt-predicted class pairs
cnn <- cnn_data_clean %>%
  select(c(group, aa_wt, aa_pred)) %>%
  mutate(wt_class = map_chr(aa_wt, calc_class)) %>%
  mutate(wt_class = fct_relevel(wt_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  mutate(pred_class = map_chr(aa_pred, calc_class)) %>%
  mutate(pred_class = fct_relevel(pred_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  select(c(group, pred_class, wt_class)) %>%
  group_by(wt_class, pred_class) %>%
  summarize(count = n())

bert <- bert_data_clean %>%
  select(c(group, aa_wt, aa_pred)) %>%
  mutate(wt_class = map_chr(aa_wt, calc_class)) %>%
  mutate(wt_class = fct_relevel(wt_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  mutate(pred_class = map_chr(aa_pred, calc_class)) %>%
  mutate(pred_class = fct_relevel(pred_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  select(c(group, pred_class, wt_class)) %>%
  group_by(wt_class, pred_class) %>%
  summarize(count = n())

esm <- esm_data_clean %>%
  select(c(group, aa_wt, aa_pred)) %>%
  mutate(wt_class = map_chr(aa_wt, calc_class)) %>%
  mutate(wt_class = fct_relevel(wt_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  mutate(pred_class = map_chr(aa_pred, calc_class)) %>%
  mutate(pred_class = fct_relevel(pred_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  select(c(group, pred_class, wt_class)) %>%
  group_by(wt_class, pred_class) %>%
  summarize(count = n())

combo <- combo_data_clean %>%
  select(c(group, aa_wt, aa_pred)) %>%
  mutate(wt_class = map_chr(aa_wt, calc_class)) %>%
  mutate(wt_class = fct_relevel(wt_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  mutate(pred_class = map_chr(aa_pred, calc_class)) %>%
  mutate(pred_class = fct_relevel(pred_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")) %>%
  select(c(group, pred_class, wt_class)) %>%
  group_by(wt_class, pred_class) %>%
  summarize(count = n())

# calculating the frequency:
cnn <- cnn %>%
  group_by(wt_class) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  rename(cnn_freq = freq) %>%
  select(-c(count))

bert <- bert %>%
  group_by(wt_class) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  rename(bert_freq = freq) %>%
  select(-c(count))

esm <- esm %>%
  group_by(wt_class) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  rename(esm_freq = freq) %>%
  select(-c(count))

combo <- combo %>%
  group_by(wt_class) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup() %>%
  rename(combo_freq = freq) %>%
  select(-c(count))




#plots:

cnn_joined <- full_join(cnn, combo)
cnn_joined <- cnn_joined %>%
  mutate(cnn_freq = ifelse(is.na(cnn_freq), 0, cnn_freq),
         combo_freq = ifelse(is.na(combo_freq), 0, combo_freq)) %>%
  mutate(freq = combo_freq - cnn_freq) 

font_size <- 10

plot_cnn_class <- ggplot() +
  geom_tile(data = cnn_joined, aes(x = wt_class, y = fct_rev(pred_class), fill = freq)) + 
  scale_fill_continuous_diverging(
    palette = "Blue-Red 3", 
    #palette = "Tropic",
    rev = TRUE, 
    p2 = 1
  ) +
  scale_x_discrete(
    name = "Wild type class",
    expand = c(0,0),
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_y_discrete(
    name = "Predicted class",
    expand = c(0,0),
    labels = c("unique", "aromatic", "positive", "negative", "small polar", "aliphatic")) +
  labs(fill = "Change in frequency \n from CNN model") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(
      color = "black", 
      size = font_size),
    axis.text.x = element_text(
      angle = 45, hjust = 1),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_cnn_class

# bert plot:
bert_joined <- full_join(bert, combo)
bert_joined <- bert_joined %>%
  mutate(bert_freq = ifelse(is.na(bert_freq), 0, bert_freq),
         combo_freq = ifelse(is.na(combo_freq), 0, combo_freq)) %>%
  mutate(freq = combo_freq - bert_freq) 

font_size <- 10

plot_bert_class <- ggplot() +
  geom_tile(data = bert_joined, aes(x = wt_class, y = fct_rev(pred_class), fill = freq)) + 
  scale_fill_continuous_diverging(
    palette = "Blue-Red 3", 
    #palette = "Tropic",
    rev = TRUE, 
    p2 = 1
  ) +
  scale_x_discrete(
    name = "Wild type class",
    expand = c(0,0),
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_y_discrete(
    name = "Predicted class",
    expand = c(0,0),
    labels = c("unique", "aromatic", "positive", "negative", "small polar", "aliphatic")) +
  labs(fill = "Change in frequency \n from BERT model") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(
      color = "black", 
      size = font_size),
    axis.text.x = element_text(
      angle = 45, hjust = 1),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_bert_class

# esm plot:
esm_joined <- full_join(esm, combo)
esm_joined <- esm_joined %>%
  mutate(esm_freq = ifelse(is.na(esm_freq), 0, esm_freq),
         combo_freq = ifelse(is.na(combo_freq), 0, combo_freq)) %>%
  mutate(freq = combo_freq - esm_freq) 

font_size <- 10

plot_esm_class <- ggplot() +
  geom_tile(data = esm_joined, aes(x = wt_class, y = fct_rev(pred_class), fill = freq)) + 
  scale_fill_continuous_diverging(
    palette = "Blue-Red 3", 
    #palette = "Tropic",
    rev = TRUE, 
    p2 = 1
  ) +
  scale_x_discrete(
    name = "Wild type class",
    expand = c(0,0),
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_y_discrete(
    name = "Predicted class",
    expand = c(0,0),
    labels = c("unique", "aromatic", "positive", "negative", "small polar", "aliphatic")) +
  labs(fill = "Change in frequency \n from ESM1b model") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(
      color = "black", 
      size = font_size),
    axis.text.x = element_text(
      angle = 45, hjust = 1),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_esm_class

abcdef <- plot_grid(plot_cnn_aa + theme(legend.position="none"), 
                    plot_cnn_class + theme(legend.position="none"), 
                    plot_bert_aa + theme(legend.position="none"), 
                    plot_bert_class + theme(legend.position="none"), 
                    plot_esm_aa + theme(legend.position="none"), 
                    plot_esm_class + theme(legend.position="none"), 
                    nrow = 3, 
                    ncol = 2, 
                    align = "hv", 
                    axis = "bl",
                    labels = c('a', 'b', 'c', 'd', 'e', 'f'))

legend <- get_legend(
  # create some space to the left of the legend
  plot_cnn_aa + theme(legend.box.margin = margin(0, 0, 0, 12)))

with_legend <- plot_grid(abcdef, legend, rel_widths = c(3, 1))
with_legend

ggsave(filename = "./analysis/figures/divergent_heatmaps.png", plot = with_legend, width = 9, height = 10)










