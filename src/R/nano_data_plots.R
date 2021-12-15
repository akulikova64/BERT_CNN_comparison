library(tidyverse)
library(cowplot)
library(data.table)
library(ggforce)
library(colorspace)
library(ggpubr)
options(scipen = 999)

# this script makes nanobody data plots

nano_data <- read.csv(file = "./output/nanobody_summary_data.csv", header=TRUE, sep=",")


# get all positions that were mutated from 3ogo in experiment
exp_mut_pos <- nano_data %>%
  filter(change == "mutation") %>%
  mutate(cnn_not_match = wt_3ogo != cnn_aa,
         bert_not_match = wt_3ogo != bert_aa) %>%
  group_by(nanobody) %>%
  summarise(freq_of_mutation_cnn = sum(cnn_not_match, na.rm = TRUE)/sum(!is.na(cnn_not_match)),
            freq_of_mutation_bert = sum(bert_not_match, na.rm = TRUE)/sum(!is.na(bert_not_match))) %>%
  pivot_longer(cols = c(freq_of_mutation_cnn, freq_of_mutation_bert), names_to = "model", values_to = "freq_of_mutation") %>%
  mutate(model_name = ifelse(model == "freq_of_mutation_cnn", "CNN", "BERT"))


bar_plot_1 <- exp_mut_pos %>%
  ggplot(aes(x = model_name, y = freq_of_mutation, fill = fct_rev(model_name), color = fct_rev(model_name))) +
  geom_violin(alpha = 0.5) +
  geom_sina(alpha = 0.5) +
  scale_fill_manual(values = c("#608046", "#785782", "#ab6f07")) +
  scale_color_manual(values = c("#1c4728", "#452d4d", "#5e3f08")) +
  scale_y_continuous(
    name = "Proportion at which the experimentaly \nmutated positions match the \nCNN/BERT suggested mutations",
    limits = c(0.1, 0.8),
    breaks = seq(0.1, 0.8, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "") +
  labs(fill = "model") +
  guides(color = FALSE) +
  theme_cowplot(12)


bar_plot_1

#-------------------------------------------------------------------------------
# Next, we are going to see if the exact mutations match, not just the positions

# get all positions that were mutated from 3ogo in experiment
exp_match_nn <- nano_data %>%
  filter(change == "mutation") %>%
  mutate(cnn_match_exp = wt_aa == cnn_aa,
         bert_match_exp = wt_aa == bert_aa) %>%
  group_by(nanobody) %>%
  summarise(freq_of_cnn_match = sum(cnn_match_exp, na.rm = TRUE)/sum(!is.na(cnn_match_exp)),
            freq_of_bert_match = sum(bert_match_exp, na.rm = TRUE)/sum(!is.na(bert_match_exp))) %>%
  pivot_longer(cols = c(freq_of_cnn_match, freq_of_bert_match), names_to = "model", values_to = "freq_of_match") %>%
  mutate(model_name = ifelse(model == "freq_of_cnn_match", "CNN", "BERT"))


bar_plot_2 <- exp_match_nn %>%
  ggplot(aes(x = model_name, y = freq_of_match, fill = fct_rev(model_name), color = fct_rev(model_name))) +
  geom_violin(alpha = 0.5) +
  geom_sina(alpha = 0.5) +
  scale_fill_manual(values = c("#608046", "#785782")) +
  scale_color_manual(values = c("#1c4728", "#452d4d")) +
  scale_y_continuous(
    name = "Proportion at which experimental \nmutations match the \nCNN/BERT suggested mutations exactly",
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "") +
  labs(fill = "model") +
  guides(color = FALSE) +
  theme_cowplot(12)

bar_plot_2


#-------------------------------------------------------------------------------
# Next, lets correlate yeild with the percent of experimentally mutated positions
# that match the positions of nn suggested mutations.

# we will have 30 points in each plot (one for each nanobody)

exp_mut_pos2 <- nano_data %>%
  filter(change == "mutation") %>%
  mutate(cnn_not_match = wt_3ogo != cnn_aa,
         bert_not_match = wt_3ogo != bert_aa) %>%
  group_by(nanobody) %>%
  summarise(freq_of_mutation_cnn = sum(cnn_not_match, na.rm = TRUE)/sum(!is.na(cnn_not_match)),
            freq_of_mutation_bert = sum(bert_not_match, na.rm = TRUE)/sum(!is.na(bert_not_match)),
            yeild = mean(yeild)) %>%
  pivot_longer(cols = c(freq_of_mutation_cnn, freq_of_mutation_bert), names_to = "model", values_to = "freq_of_mutation") %>%
  mutate(model_name = ifelse(model == "freq_of_mutation_cnn", "CNN", "BERT"))

cor_plot_cnn_1 <- exp_mut_pos2 %>%
  filter(model_name == "CNN") %>%
  ggplot(aes(x = freq_of_mutation, y = yeild)) +
  geom_point(color = "#1c4728") +
  scale_y_continuous(
    name = "Yeild",
    limits = c(0, 7),
    #breaks = seq(0, 0.5, by = 0.1),
    expand = c(0, 0)) +
  scale_x_continuous(
    name = "Freq at which the experimentally \nmutated positions match \nCNN suggestions",
    limits = c(0.18, 0.28),
    #breaks = seq(0, 0.5, by = 0.1),
    expand = c(0, 0)) +
  theme_cowplot(12) +
  geom_smooth(method=lm, se=FALSE, color = "maroon") +
  stat_regline_equation(label.y = 6.3, label.x = 0.185, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5.9, label.x = 0.185, aes(label = ..rr.label..))

cor_plot_cnn_1

cor_plot_bert_1 <- exp_mut_pos2 %>%
  filter(model_name == "BERT") %>%
  ggplot(aes(x = freq_of_mutation, y = yeild)) +
  geom_point(color = "#452d4d") +
  scale_y_continuous(
    name = "Yeild",
    limits = c(0, 7),
    #breaks = seq(0, 0.5, by = 0.1),
    expand = c(0, 0)) +
  scale_x_continuous(
    name = "Freq at which the experimentally \nmutated positions match \nBERT suggestions",
    limits = c(0.5, 0.8),
    #breaks = seq(0, 0.5, by = 0.1),
    expand = c(0, 0)) +
  theme_cowplot(12) +
  geom_smooth(method=lm, se=FALSE, color = "maroon") +
  stat_regline_equation(label.y = 6.3, label.x = 0.51, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5.9, label.x = 0.51, aes(label = ..rr.label..))

cor_plot_bert_1

ggsave(filename = "./analysis/figures/nano_data_bar1.png", plot = bar_plot_1 , width = 6, height = 5)
ggsave(filename = "./analysis/figures/nano_data_bar2.png", plot = bar_plot_2 , width = 6, height = 5)
ggsave(filename = "./analysis/figures/nano_data_yeild_cor_cnn.png", plot = cor_plot_cnn_1 , width = 6, height = 5)
ggsave(filename = "./analysis/figures/nano_data_yeild_cor_bert.png", plot = cor_plot_bert_1 , width = 6, height = 5)








