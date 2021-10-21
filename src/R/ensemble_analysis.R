library(tidyverse)
library(cowplot)
library(data.table)
library(ggforce)
library(colorspace)
library(gridExtra)
library(grid)
options(scipen = 999)

# ensemble net analysis:
ensemble <- read.csv(file = "./output/ensemble_predictions.csv", header=TRUE, sep=",")
trans_data_all <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")
cnn_data_all <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")


# clean ensemble data:
ens <- ensemble %>%
  select(c(gene, position, aa_wt, aa_pred, wt_prob, pred_prob)) %>%
  mutate(group = "ensemble")

# clean cnn data:
cnn_data <- cnn_data_all %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) %>%
  mutate(group = "cnn") %>%
  rename(wt_prob = freq_wt,
         pred_prob = freq_predicted,
         aa_pred = aa_predicted)

# clean tranformer data:
trans_data <- trans_data_all %>%
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) %>%
  mutate(group = "transformer") %>%
  select(c(gene, position, aa_pred, aa_wt, pred_prob, wt_prob, group))


# check that genes are present across all three groups (eliminate genes that are not)
t <- trans_data %>%
  group_by(gene) %>%
  count() %>%
  select(gene)
e <- ens %>%
  group_by(gene) %>%
  count() %>%
  select(gene)
c <- cnn_data %>%
  group_by(gene) %>%
  count() %>%
  select(gene)

inner <- inner_join(t, e)
inner2 <- inner_join(c, inner)

# now bind rows
joined <- rbind(cnn_data, trans_data, ens)

joined2 <- inner_join(joined, inner2)

match_wt <- joined2 %>%
  mutate(match_predict_wt = aa_pred == aa_wt)

#data entries where the predicted amino acid matches the wt
stats_1 <- match_wt %>%
  group_by(gene, group) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)),
            mean_conf = mean(pred_prob))

means <- stats_1 %>%
  group_by(group) %>%
  summarise(mean = mean(freq_predict_wt))

#cnn: 75.2% Accuracy
#ensemble: 74.9% Accuracy
#transformer: 67.7% Accuracy

stats_1 <- stats_1 %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(freq_predict_wt * (group == "ensemble")),
    dx = rnorm(n(), mean = 0, sd = .05),
    dy = rnorm(n(), mean = 0, sd = .05),
    x_value = as.numeric(factor(group)))  

plot <- stats_1 %>%
  ggplot(aes(x = fct_relevel(group, "cnn", "ensemble", "transformer"), y = freq_predict_wt)) +
  geom_path(
    aes(x = as.numeric(factor(group))+dx, y = freq_predict_wt+dy, group = gene, color = color_y),
    size = 0.25) +
  geom_point(
    aes(x = as.numeric(factor(group))+dx, y = freq_predict_wt+dy, group = gene, fill = color_y),
    shape = 21, 
    color = "black",
    size = 2) +
  scale_x_continuous(
    name = "",
    #limits = c(0.7,2.3),
    labels = c("CNN", "Ensemble \n (average of prob.)", "Transformer"),
    breaks = (seq(from = 1.0, to = 3.0, by = 1))
  ) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.1, 1.0),
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

ggsave(filename = paste0("./analysis/figures/accuracy_with_ensemble.png"), plot = plot, width = 7.5, height = 6)

#========================================================================================
# Instead of averaging the probabilities, I took the highest confidence prediction
# for each position between the two models.
#========================================================================================

ensemble_conf <- read.csv(file = "./output/ens_pred_highest_confidence.csv", header=TRUE, sep=",")

# clean ensemble data:
ens <- ensemble_conf %>%
  select(c(gene, position, aa_wt, aa_pred, wt_prob, pred_prob)) %>%
  mutate(group = "ensemble")

# get "cnn_data" and "trans_data" from above, the proceed with the next step

# check that genes are present across all three groups (eliminate genes that are not)
t <- trans_data %>%
  group_by(gene) %>%
  count() %>%
  select(gene)
e <- ens %>%
  group_by(gene) %>%
  count() %>%
  select(gene)
c <- cnn_data %>%
  group_by(gene) %>%
  count() %>%
  select(gene)

inner <- inner_join(t, e)
inner2 <- inner_join(c, inner)

# now bind rows
joined <- rbind(cnn_data, trans_data, ens)

joined2 <- inner_join(joined, inner2)

match_wt <- joined2 %>%
  mutate(match_predict_wt = aa_pred == aa_wt)

#data entries where the predicted amino acid matches the wt
stats_1 <- match_wt %>%
  group_by(gene, group) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)),
            mean_conf = mean(pred_prob))

means <- stats_1 %>%
  group_by(group) %>%
  summarise(mean = mean(freq_predict_wt))

#cnn: 75.2% Accuracy
#ensemble: 74.0% Accuracy
#transformer: 67.7% Accuracy

stats_1 <- stats_1 %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(freq_predict_wt * (group == "ensemble")),
    dx = rnorm(n(), mean = 0, sd = .05),
    dy = rnorm(n(), mean = 0, sd = .05),
    x_value = as.numeric(factor(group)))  

plot2 <- stats_1 %>%
  ggplot(aes(x = fct_relevel(group, "cnn", "ensemble", "transformer"), y = freq_predict_wt)) +
  geom_path(
    aes(x = as.numeric(factor(group))+dx, y = freq_predict_wt+dy, group = gene, color = color_y),
    size = 0.25) +
  geom_point(
    aes(x = as.numeric(factor(group))+dx, y = freq_predict_wt+dy, group = gene, fill = color_y),
    shape = 21, 
    color = "black",
    size = 2) +
  scale_x_continuous(
    name = "",
    #limits = c(0.7,2.3),
    labels = c("CNN", "Ensemble \n (Highest conf. chosen)", "Transformer"),
    breaks = (seq(from = 1.0, to = 3.0, by = 1))
  ) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.1, 1.0),
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

plot2

ggsave(filename = paste0("./analysis/figures/accurac_ens_confidence.png"), plot = plot2, width = 7.5, height = 6)




