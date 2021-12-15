# This will be the figure with accuracy for each model

library(tidyverse)
library(cowplot)
library(data.table)
library(ggforce)
library(colorspace)
library(ggpubr)
options(scipen = 999)

# this script compares BERT to CNN for accuracy in predictions
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

match_wt <- joined3 %>%
  mutate(match_predict_wt = aa_pred == aa_wt)

#data entries where the predicted amino acid matches the wt
stats_1 <- match_wt %>%
  group_by(gene, group) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)),
            mean_conf = mean(freq_pred))

means <- stats_1 %>%
  group_by(group) %>%
  summarise(mean = mean(freq_predict_wt))
# CNN: 75.8% accuracy
# BERT transformer: 71.6% accuracy
# ESM1b transformer: 60.58% accuracy

joined_stats <- stats_1 %>%
  select(-c(mean_conf)) %>%
  group_by(group) %>%
  summarise(estimate = mean(freq_predict_wt),
            std_error = sd(freq_predict_wt)/sqrt(length(freq_predict_wt)))


# making figure:

accuracy_all <- stats_1 %>%
  ggplot(aes(y = freq_predict_wt, 
             x = fct_relevel(group, "ESM1b", "BERT", "CNN", "Combined"), 
             fill = fct_relevel(group, "ESM1b", "BERT", "CNN", "Combined"), 
             color = fct_relevel(group, "ESM1b", "BERT", "CNN", "Combined"))) +
  geom_violin(alpha = 0.6, 
              size = 0.4, 
              bw = 0.02, 
              position = position_dodge(width = 1),
              width = 1,
              #scale = "width"
  ) + 
  geom_sina(alpha = 0.5,
            size = 1) +
  geom_pointrange(data = joined_stats, aes(x = group,
                                           y = estimate,
                                           ymin = estimate - 1.96*std_error,
                                           ymax = estimate + 1.96*std_error),
                  color = "black", 
                  alpha = 0.7, 
                  size = 0.4,
                  position = position_dodge(width = 1)) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(14) + 
  theme(plot.title = element_text(hjust = 0, size = 14), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        #panel.grid.minor.y = element_line(color = "grey92", size=0.5),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none") +
  labs(fill = "") +
  guides(color = FALSE) +
  scale_fill_manual(values = c("#ab6f07", "#785782", "#608046", "#66102a")) +
  #scale_fill_manual(values = c("#46b063", "#a169b3", "#db9214")) +
  scale_color_manual(values = c("#5e3f08", "#452d4d", "#1c4728", "#3b0414")) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "")

accuracy_all

#----------------------------------------------------------------------------
# next plots are correlation plots:

#lets plot BERT accuracy vs ESM1b accuracy. I suspect they are not correlated.

wider_stats <- stats_1 %>%
  select(c(group, freq_predict_wt)) %>%
  pivot_wider(names_from = group, values_from = freq_predict_wt)

#cor <- cor(wider_stats$BERT ~ wider_stats$ESM1b)

cor_esm_cnn <- wider_stats %>%
  ggplot(aes(x = ESM1b, y = CNN)) +
  geom_point(
    shape = 21, 
    color = "black",
    fill = "#38a3c7",
    alpha = 0.5,
    size = 2) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "maroon"
  ) +
  #geom_smooth(method=lm, se=FALSE, color = "maroon") +
  scale_y_continuous(
    name = "CNN",
    limits = c(0.0, 1.01),
    breaks = (seq(from = 0.0, to = 1.0, by = 0.2)),
    expand = c(0, 0)
  ) +
  #stat_regline_equation(label.y = 0.95, label.x = 0.02, aes(label = ..eq.label..)) +
  #stat_regline_equation(label.y = 0.89, label.x = 0.02, aes(label = ..rr.label..)) +
  scale_x_continuous(
    name = "ESM1b",
    limits = c(0.0, 1.01),
    breaks = seq(from = 0.0, to = 1.0, by = 0.2),
    expand = c(0, 0)) +
  theme_bw(14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    panel.grid.minor = element_blank())

cor_esm_cnn

cor_bert_cnn <- wider_stats %>%
  ggplot(aes(x = BERT, y = CNN)) +
  geom_point(
    shape = 21, 
    color = "black",
    fill = "#38a3c7",
    alpha = 0.5,
    size = 2) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "maroon"
  ) +
  #geom_smooth(method=lm, se=FALSE, color = "maroon") +
  scale_y_continuous(
    name = "CNN",
    limits = c(0.0, 1.01),
    breaks = (seq(from = 0.0, to = 1.0, by = 0.2)),
    expand = c(0, 0)
  ) +
  #stat_regline_equation(label.y = 0.95, label.x = 0.02, aes(label = ..eq.label..)) +
  #stat_regline_equation(label.y = 0.89, label.x = 0.02, aes(label = ..rr.label..)) +
  scale_x_continuous(
    name = "BERT",
    limits = c(0.0, 1.01),
    breaks = seq(from = 0.0, to = 1.0, by = 0.2),
    expand = c(0, 0)) +
  theme_bw(14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    panel.grid.minor = element_blank())

cor_bert_cnn

cor_bert_esm <- wider_stats %>%
  ggplot(aes(x = BERT, y = ESM1b)) +
  geom_point(
    shape = 21, 
    color = "black",
    fill = "#38a3c7",
    alpha = 0.5,
    size = 2) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "maroon"
  ) +
  #geom_smooth(method=lm, se=FALSE, color = "maroon") +
  scale_y_continuous(
    name = "ESM1b",
    limits = c(0.0, 1.01),
    breaks = (seq(from = 0.0, to = 1.0, by = 0.2)),
    expand = c(0, 0)
  ) +
  #stat_regline_equation(label.y = 0.95, label.x = 0.02, aes(label = ..eq.label..)) +
  #stat_regline_equation(label.y = 0.89, label.x = 0.02, aes(label = ..rr.label..)) +
  scale_x_continuous(
    name = "BERT",
    limits = c(0.0, 1.01),
    breaks = seq(from = 0.0, to = 1.0, by = 0.2),
    expand = c(0, 0)) +
  theme_bw(14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    panel.grid.minor = element_blank())

cor_bert_esm

a <- plot_grid(accuracy_all, nrow = 1, align="h", labels = c('a'))
bcd <- plot_grid(cor_bert_cnn, cor_esm_cnn, cor_bert_esm, nrow = 1, align="h", labels = c('b', 'c', 'd'))
abcd <- plot_grid(a, bcd, nrow = 2, rel_heights = c(3, 2), align = "v", labels = c('', ''))
abcd
ggsave(filename = "./analysis/figures/accuracy_fig4.png", plot = abcd, width = 9, height = 7)


#------------------------------------------------------------------------------------------------
# Let's look at the same plot, but fillter out all PSICOV proteins not found across all 4 models
#------------------------------------------------------------------------------------------------

genes <- data.frame(unique(combo_data$gene)) %>%
  rename(gene = unique.combo_data.gene.)

filtered <- inner_join(joined3, genes)
  
match_wt <- filtered %>%
  mutate(match_predict_wt = aa_pred == aa_wt)

#data entries where the predicted amino acid matches the wt
stats_1 <- match_wt %>%
  group_by(gene, group) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)),
            mean_conf = mean(freq_pred))

means <- stats_1 %>%
  group_by(group) %>%
  summarise(mean = mean(freq_predict_wt))
# CNN: 75.8% accuracy
# BERT transformer: 71.6% accuracy
# ESM1b transformer: 60.58% accuracy

joined_stats <- stats_1 %>%
  select(-c(mean_conf)) %>%
  group_by(group) %>%
  summarise(estimate = mean(freq_predict_wt),
            std_error = sd(freq_predict_wt)/sqrt(length(freq_predict_wt)))


# making figure:

accuracy_all <- stats_1 %>%
  ggplot(aes(y = freq_predict_wt, 
             x = fct_relevel(group, "ESM1b", "BERT", "CNN", "Combined"), 
             fill = fct_relevel(group, "ESM1b", "BERT", "CNN", "Combined"), 
             color = fct_relevel(group, "ESM1b", "BERT", "CNN", "Combined"))) +
  geom_violin(alpha = 0.6, 
              size = 0.4, 
              #bw = 0.02,
              bw = 0.03,
              position = position_dodge(width = 1),
              width = 1,
              #scale = "width"
  ) + 
  geom_sina(alpha = 0.3,
            size = 2) +
  geom_pointrange(data = joined_stats, aes(x = group,
                                           y = estimate,
                                           ymin = estimate - 1.96*std_error,
                                           ymax = estimate + 1.96*std_error),
                  color = "black", 
                  alpha = 0.7, 
                  size = 0.4,
                  position = position_dodge(width = 1)) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(14) + 
  theme(plot.title = element_text(hjust = 0, size = 14), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        #panel.grid.minor.y = element_line(color = "grey92", size=0.5),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none") +
  labs(fill = "") +
  guides(color = FALSE) +
  scale_fill_manual(values = c("#ab6f07", "#785782", "#608046", "#66102a")) +
  #scale_fill_manual(values = c("#46b063", "#a169b3", "#db9214")) +
  scale_color_manual(values = c("#5e3f08", "#452d4d", "#1c4728", "#3b0414")) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "")

accuracy_all

#----------------------------------------------------------------------------
# next plots are correlation plots:

#lets plot BERT accuracy vs ESM1b accuracy. I suspect they are not correlated.

wider_stats <- stats_1 %>%
  select(c(group, freq_predict_wt)) %>%
  pivot_wider(names_from = group, values_from = freq_predict_wt)

#cor <- cor(wider_stats$BERT ~ wider_stats$ESM1b)

cor_esm_cnn <- wider_stats %>%
  ggplot(aes(x = ESM1b, y = CNN)) +
  geom_point(
    shape = 21, 
    color = "black",
    fill = "#38a3c7",
    alpha = 0.5,
    size = 2) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "maroon"
  ) +
  #geom_smooth(method=lm, se=FALSE, color = "maroon") +
  scale_y_continuous(
    name = "CNN",
    limits = c(0.0, 1.01),
    breaks = (seq(from = 0.0, to = 1.0, by = 0.2)),
    expand = c(0, 0)
  ) +
  #stat_regline_equation(label.y = 0.95, label.x = 0.02, aes(label = ..eq.label..)) +
  #stat_regline_equation(label.y = 0.89, label.x = 0.02, aes(label = ..rr.label..)) +
  scale_x_continuous(
    name = "ESM1b",
    limits = c(0.0, 1.01),
    breaks = seq(from = 0.0, to = 1.0, by = 0.2),
    expand = c(0, 0)) +
  theme_bw(14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    panel.grid.minor = element_blank())

cor_esm_cnn

cor_bert_cnn <- wider_stats %>%
  ggplot(aes(x = BERT, y = CNN)) +
  geom_point(
    shape = 21, 
    color = "black",
    fill = "#38a3c7",
    alpha = 0.5,
    size = 2) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "maroon"
  ) +
  #geom_smooth(method=lm, se=FALSE, color = "maroon") +
  scale_y_continuous(
    name = "CNN",
    limits = c(0.0, 1.01),
    breaks = (seq(from = 0.0, to = 1.0, by = 0.2)),
    expand = c(0, 0)
  ) +
  #stat_regline_equation(label.y = 0.95, label.x = 0.02, aes(label = ..eq.label..)) +
  #stat_regline_equation(label.y = 0.89, label.x = 0.02, aes(label = ..rr.label..)) +
  scale_x_continuous(
    name = "BERT",
    limits = c(0.0, 1.01),
    breaks = seq(from = 0.0, to = 1.0, by = 0.2),
    expand = c(0, 0)) +
  theme_bw(14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    panel.grid.minor = element_blank())

cor_bert_cnn

cor_bert_esm <- wider_stats %>%
  ggplot(aes(x = BERT, y = ESM1b)) +
  geom_point(
    shape = 21, 
    color = "black",
    fill = "#38a3c7",
    alpha = 0.5,
    size = 2) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "maroon"
  ) +
  #geom_smooth(method=lm, se=FALSE, color = "maroon") +
  scale_y_continuous(
    name = "ESM1b",
    limits = c(0.0, 1.01),
    breaks = (seq(from = 0.0, to = 1.0, by = 0.2)),
    expand = c(0, 0)
  ) +
  #stat_regline_equation(label.y = 0.95, label.x = 0.02, aes(label = ..eq.label..)) +
  #stat_regline_equation(label.y = 0.89, label.x = 0.02, aes(label = ..rr.label..)) +
  scale_x_continuous(
    name = "BERT",
    limits = c(0.0, 1.01),
    breaks = seq(from = 0.0, to = 1.0, by = 0.2),
    expand = c(0, 0)) +
  theme_bw(14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    panel.grid.minor = element_blank())

cor_bert_esm

a <- plot_grid(accuracy_all, nrow = 1, align="h", labels = c('a'))
bcd <- plot_grid(cor_bert_cnn, cor_esm_cnn, cor_bert_esm, nrow = 1, align="h", labels = c('b', 'c', 'd'))
abcd <- plot_grid(a, bcd, nrow = 2, rel_heights = c(3, 2), align = "v", labels = c('', ''))
abcd
ggsave(filename = "./analysis/figures/accuracy_fig_filtered.png", plot = abcd, width = 9, height = 7)


  
  
