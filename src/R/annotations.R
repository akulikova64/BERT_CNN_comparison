library(tidyverse)
library(cowplot)
library(data.table)
library(ggforce)
library(colorspace)
library(gridExtra)
library(grid)
options(scipen = 999)



annot <- read.csv(file = "./data/test.csv", header=TRUE, sep=",")
trans_data <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")
cnn_data_all <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

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
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)),
            mean_conf = mean(freq_pred))

means <- stats_1 %>%
  group_by(group) %>%
  summarise(mean = mean(freq_predict_wt))
#cnn: 75.8% Accuracy
#transformer: 71.6% accuracy

with_annot1 <- left_join(stats_1, annot)

for_diag1 <- with_annot1 %>%
  select(c(gene, group, freq_predict_wt, structural_class)) %>%
  pivot_wider(names_from = group, values_from = freq_predict_wt)

for_diag <- stats_1 %>%
  select(c(gene, group, freq_predict_wt)) %>%
  pivot_wider(names_from = group, values_from = freq_predict_wt)

#fills = c("#bf8040", "#ac5396", "#70adc2", "#748f3d", "#cc5933", "#7070c2")
fills = c("#f51da6", "#ffbb00", "#0cba00", "#0048e3")  
colors = c("#b00270", "bf8c00", "087500", "002d8f")

plot <- for_diag1 %>%
  na.omit() %>%
  ggplot(aes(x = transformer, y = cnn, fill = fct_rev(structural_class), color = fct_rev(structural_class))) +
  geom_point(
    shape = 21, 
    color = "black",
    size = 3) +
  geom_abline(slope = 1, 
              intercept = 0,
              color = "maroon"
              ) +
  scale_x_continuous(
    name = "Transformer Accuracy",
    limits = c(0.0, 1.01),
    breaks = (seq(from = 0.0, to = 1.0, by = 0.1)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "CNN Accuracy",
    limits = c(0.0, 1.01),
    breaks = seq(from = 0.0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_fill_manual(values = fills) +
  scale_color_manual(values = colors) +
  labs(fill = "structural class") +
  theme_bw(16) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 16),
    panel.grid.minor = element_blank())

plot

#ggsave(filename = paste0("./analysis/figures/accuracy_cnn_bert_diag_struc.png"), plot = plot, width = 9, height = 6)


#----------------------------------------------------------
#lets look at structural class vs accuracy from annotations

for_str_class <- for_diag %>%
  pivot_longer(cols = c(cnn, transformer), names_to = "group", values_to = "accuracy")

with_annot <- left_join(for_str_class, annot)

fills <- c("#a3ba4e", "#9880b0")
colors <- c("#37420e", "#30154d")

struc_plot <- with_annot %>%
  na.omit() %>%
  ggplot(aes(x = structural_class, y = accuracy, fill = group, color = group)) +
  geom_violin(position = "dodge") +
  geom_sina(position = "dodge") +
  theme_cowplot(16) +
  scale_color_manual(values = colors, ) +
  scale_fill_manual(values = fills, labels = c("cnn", "transformer")) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.2, 1.01),
    breaks = seq(from = 0.2, to = 1.0, by = 0.1),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    name = "Structural Class",
    expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "right")

struc_plot

#ggsave(filename = paste0("./analysis/figures/structural_class_acc.png"), plot = struc_plot, width = 12, height = 6)

#--------------------------------------------------------------------------------------------
# now lets look at the proteins where one model is significantly better/worse than the other


cnn_poor <- with_annot %>%
  pivot_wider(names_from = group, values_from = accuracy) %>%
  filter(cnn < 0.7 & transformer > 0.7) %>%
  select(-c(structural_class, protein_type))

trans_poor <- with_annot %>%
  pivot_wider(names_from = group, values_from = accuracy) %>%
  filter(transformer < 0.4 & cnn > 0.7) %>%
  select(-c(structural_class, protein_type))


table1 <- table(cnn_poor)

#ggsave(filename = paste0("./analysis/figures/cnn_poor_table.png"), plot = table1, width = 12, height = 6)
