
# looking if the mispredictions are inspired by what is most abundant in nature.

#=====================================================================================
# Frequency of the predicted in MSA as a function of CNN confidence bins (mispredictions only)
#=====================================================================================

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
trans_bert_data <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")
trans_data_esm <- read.csv(file = "./output/PSICOV_ESM1b_predictions.csv", header=TRUE, sep=",")

natural_data <- read.csv(file = "./data/natural_max_freq_all.csv", header=TRUE, sep=",")
natural_var <- read.csv(file = "./output/stats_align_all.csv", header=TRUE, sep=",")

transf_esm <- trans_data_esm %>% # choose which transformer model to use here
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) %>%
  rename(freq_pred = pred_prob,
         freq_wt = wt_prob,
         aa_wt = aa_wt,
         aa_pred = aa_pred)

transf_bert <- trans_bert_data %>% # choose which transformer model to use here
  mutate(pred_prob = pmax(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)) %>%
  pivot_longer(cols = c(A:Y), names_to = "aa", values_to = "freq") %>%
  mutate(aa_pred = ifelse(pred_prob == freq, aa, NA)) %>%
  na.omit() %>%
  select(-c(aa, freq, row)) %>%
  rename(freq_pred = pred_prob,
         freq_wt = wt_prob,
         aa_wt = aa_wt,
         aa_pred = aa_pred)

#clean cnn data:

cnn_data <- cnn_data_all %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) %>%
  rename(aa_pred = aa_predicted,
         aa_wt = aa_wt,
         freq_pred = freq_predicted,
         freq_wt = freq_wt)
  #mutate(group = "CNN")

#clean natural data:

nat_data_clean <- natural_data %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa, freq))
  # rename(aa_pred = aa_natural_max,
  #        aa_wt = aa_natural_wt,
  #        freq_pred = freq_natural_max,
  #        freq_wt = freq_natural_wt) %>%
  #mutate(group = "natural")


# clean esm data:
esm_data_unjoined <- joined_data %>%
  select(-c(aa_wt_bert, aa_pred_bert, freq_wt_bert, freq_pred_bert)) %>%
  rename(freq_pred = freq_pred_esm,
         aa_pred = aa_pred_esm,
         freq_wt = freq_wt_esm,
         aa_wt = aa_wt_esm)
  #mutate(group = "ESM1b")

# clean bert data:
bert_data_unjoined <- joined_data %>%
  select(-c(aa_wt_esm, aa_pred_esm, freq_wt_esm, freq_pred_esm)) %>%
  rename(freq_pred = freq_pred_bert,
         aa_pred = aa_pred_bert,
         freq_wt = freq_wt_bert,
         aa_wt = aa_wt_bert) 
  #mutate(group = "BERT")

joined_data_cnn <- inner_join(cnn_data, nat_data_clean)
joined_data_esm <- inner_join(transf_esm, nat_data_clean)
joined_data_bert <- inner_join(transf_bert, nat_data_clean)

#----------first, let's get the cnn data-----------------------
joined_data_trimmed <- joined_data_cnn %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))


joined_data2 <- inner_join(joined_data_trimmed, natural_var)

joined_data_clean <- joined_data2 %>%
  select(c(gene, position, freq_pred, aa_pred, aa_wt, q_A:q_V))

mismatches <- joined_data_clean %>%
  filter(aa_pred != aa_wt) %>%
  select(-aa_wt)

# renaming the natural frequencies
mismatches2 <- mismatches %>%
  pivot_longer(cols = c(q_A:q_V), names_to = "nat_aa", values_to = "nat_freq") %>%
  mutate(nat_aa = substr(nat_aa, 3, 3)) %>%
  filter(aa_pred == nat_aa)


get_pred_bin <- function(x) {
  
  if (x > 0 & x <= 0.2) {
    return("(0-0.2]")
  }
  else if (x > 0.2 & x <= 0.4) {
    return("(0.2-0.4]")
  }
  else if (x > 0.4 & x <= 0.6) {
    return("(0.4-0.6]")
  }
  else if (x > 0.6 & x <= 0.8) {
    return("(0.6-0.8]")
  }
  else if (x > 0.8 & x <= 1.0) {
    return("(0.8-1.0]")
  }
  else {
    return(NA)
  }
}


for_plot_cnn <- mismatches2 %>%
  mutate(pred_bin = map_chr(freq_pred, get_pred_bin)) %>%
  select(c(gene, position, nat_freq, pred_bin)) %>%
  na.omit() %>%
  mutate(group = "CNN")



data_summary <- function(x) {
  m <- mean(x)
  ymin <- y-(sd(x)/sqrt(length(x)))
  ymax <- y+(sd(x)/sqrt(length(x)))
  return(c(y=m,ymin=ymin,ymax=ymax))
}

stat_data_cnn <- for_plot_cnn %>%
  select(-c(position, gene)) %>%
  group_by(pred_bin) %>%
  summarise(estimate = mean(nat_freq),
            std_error = sd(nat_freq)/sqrt(length(nat_freq))) %>%
  mutate(group = "CNN")

#---------now lets get esm1b data-------------------------------

joined_data_trimmed <- joined_data_esm %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))


joined_data2 <- inner_join(joined_data_trimmed, natural_var)

joined_data_clean <- joined_data2 %>%
  select(c(gene, position, freq_pred, aa_pred, aa_wt, q_A:q_V))

mismatches <- joined_data_clean %>%
  filter(aa_pred != aa_wt) %>%
  select(-aa_wt)

# renaming the natural frequencies
mismatches2 <- mismatches %>%
  pivot_longer(cols = c(q_A:q_V), names_to = "nat_aa", values_to = "nat_freq") %>%
  mutate(nat_aa = substr(nat_aa, 3, 3)) %>%
  filter(aa_pred == nat_aa)

for_plot_esm <- mismatches2 %>%
  mutate(pred_bin = map_chr(freq_pred, get_pred_bin)) %>%
  select(c(gene, position, nat_freq, pred_bin)) %>%
  na.omit() %>%
  mutate(group = "ESM1b")

stat_data_esm <- for_plot_esm %>%
  select(-c(position, gene)) %>%
  group_by(pred_bin) %>%
  summarise(estimate = mean(nat_freq),
            std_error = sd(nat_freq)/sqrt(length(nat_freq))) %>%
  mutate(group = "ESM1b")

#---------now lets get BERT data-------------------------------

joined_data_trimmed <- joined_data_bert %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))


joined_data2 <- inner_join(joined_data_trimmed, natural_var)

joined_data_clean <- joined_data2 %>%
  select(c(gene, position, freq_pred, aa_pred, aa_wt, q_A:q_V))

mismatches <- joined_data_clean %>%
  filter(aa_pred != aa_wt) %>%
  select(-aa_wt)

# renaming the natural frequencies
mismatches2 <- mismatches %>%
  pivot_longer(cols = c(q_A:q_V), names_to = "nat_aa", values_to = "nat_freq") %>%
  mutate(nat_aa = substr(nat_aa, 3, 3)) %>%
  filter(aa_pred == nat_aa)

for_plot_bert <- mismatches2 %>%
  mutate(pred_bin = map_chr(freq_pred, get_pred_bin)) %>%
  select(c(gene, position, nat_freq, pred_bin)) %>%
  na.omit() %>%
  mutate(group = "BERT")

stat_data_bert <- for_plot_bert %>%
  select(-c(position, gene)) %>%
  group_by(pred_bin) %>%
  summarise(estimate = mean(nat_freq),
            std_error = sd(nat_freq)/sqrt(length(nat_freq))) %>%
  mutate(group = "BERT")

#------combinding all three datasets----------
joined_final <- rbind(for_plot_cnn, for_plot_bert, for_plot_esm)

joined_stats <- rbind(stat_data_cnn, stat_data_bert, stat_data_esm)

#fill = "#99a88a", color = "#4d5841"
# ploting the data:
plot_nat_conf <- joined_final %>%
  ggplot(aes(y = nat_freq, x = pred_bin, fill = fct_relevel(group, "CNN", "BERT", "ESM1b"))) +
  geom_violin(alpha = 0.6, 
              size = 0.4, 
              bw = 0.02, 
              position = position_dodge(width = 0.6),
              scale = "width") + 
  geom_pointrange(data = joined_stats, aes(x = pred_bin,
                                        y = estimate,
                                        ymin = estimate - 1.96*std_error,
                                        ymax = estimate + 1.96*std_error),
                  color = "black", 
                  alpha = 0.7, 
                  size = 0.4,
                  position = position_dodge(width = 0.6)) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(14) + 
  theme(plot.title = element_text(hjust = 0, size = 14), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        panel.grid.minor.y = element_line(color = "grey92", size=0.5),
        legend.position = "right") +
  labs(fill = "model") +
  scale_fill_manual(values = c("#32a852", "#9d46b8", "orange")) +
  #scale_color_manual(values = c("#154221", "#3e184a", "orange")) +
  scale_y_continuous(
    name = "Amino acid frequency",
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "CNN confidence")

plot_nat_conf

#using standard_error
#ggsave(filename = "./analysis/figures/std_error.png", plot = plot_nat_conf, width = 8, height = 5)


ggsave(filename = "./analysis/figures/nat_freq_vs_pred_bert_esm.png", plot = plot_nat_conf, width = 8.5, height = 4.5)
