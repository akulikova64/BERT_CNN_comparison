# secondary structure analysis cnn vs transformer


trans_data <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")
cnn_data_all <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
sec_struc <- read.csv(file = "./output/second_struc.csv", header=TRUE, sep=",")

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

joined2 <- rbind(cnn_data_unjoined, transf_data_unjoined)

match_wt <- joined2 %>%
  mutate(match_predict_wt = aa_pred == aa_wt) %>%
  select(c(gene, position, freq_pred, match_predict_wt, group))

trans_only <- match_wt %>%
  filter(group == "transformer") %>%
  select(-group)

#join with secondary structure:
with_struc <- trans_only %>%
  inner_join(sec_struc) %>%
  mutate(second_struc = ifelse(second_struc == "Coil", "Random coil", second_struc),
         second_struc = ifelse(second_struc == "Strand", "Beta strand", second_struc),
         second_struc = ifelse(second_struc == "AlphaHelix", "Alpha helix", second_struc),
         second_struc = ifelse(second_struc == "310Helix" | second_struc == "PiHelix", "Noncanon. helix", second_struc))

#lets get accuracy per protein:
#data entries where the predicted amino acid matches the wt
stats_1 <- with_struc %>%
  group_by(gene) %>%
  mutate(accuracy = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)))

# adding an accuracy bin:
get_acc_bin <- function(x) {
  
  if (x > 0.2 & x <= 0.4) {
    return("(0.2-0.4]")
  }
  else if (x > 0.4 & x <= 0.61) {
    return("(0.4-0.6]")
  }
  else if (x > 0.61 & x <= 0.8) {
    return("(0.6-0.8]")
  } 
  else if (x > 0.8 & x <= 1.0) {
    return("(0.8-1.0]")
  }
  NA_character_
}

get_bins <- function(acc_bin) {
  
  if (acc_bin == "(0.2-0.4]") {
    return("low accuracy \n proteins")
  }
  if (acc_bin == "(0.8-1.0]") {
    return("high accuracy \n proteins")
  }
  else {
    return(NA)
  }
  NA_character_
}

stats_new <- stats_1 %>%
  mutate(acc_bin = map_chr(accuracy, get_acc_bin)) 

binned_data <- stats_new %>%
  mutate(group = map_chr(acc_bin, get_bins)) %>%
  na.omit()

# finds the count of each structure per bin:
struc_counts <- binned_data %>%
  select(c(gene, position, second_struc, group)) %>%
  group_by(group) %>%
  count(second_struc) %>%
  mutate(struc_count = n) %>%
  select(-n) %>%
  ungroup()


# now I need to add up all aa within each  bin to get bin totals and append this to the aa_counts
for_barplot <- struc_counts %>%
  group_by(group) %>%
  mutate(bin_count = sum(struc_count)) %>%
  ungroup()


#---------------
  
for_barplot_2 <- for_barplot %>%
  mutate(freq = struc_count/bin_count) %>%
  select(-c(struc_count, bin_count)) %>% 
  mutate(second_struc = fct_rev(fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge")))

struc_fills = c("#bf8040", "#ac5396", "#70adc2", "#748f3d", "#cc5933", "#7070c2")

struc_mispr_plot <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.36),
    breaks = seq(0.0, 0.35, by = 0.05),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  facet_wrap(vars(group)) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

struc_mispr_plot

ggsave(filename = "./analysis/figures/struc_high_low_transf.png", plot = struc_mispr_plot, width = 9, height = 5)

#----------------------------------------
# getting the odds ratio:

for_odds <- for_barplot_2 %>%
  pivot_wider(names_from = group, values_from = freq) %>%
  mutate(odds = get('low accuracy \n proteins')/get('high accuracy \n proteins'))

plot_odds <- for_odds %>%
  ggplot(aes(x = odds, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = round(odds, 2), hjust = -0.05)) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_log10(
    name = "Odds ratio (log scale)"
  ) +
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_odds

ggsave(filename = "./analysis/figures/struc_odds_trans.png", plot = plot_odds, width = 9, height = 5)

#---------------------------------------------------------------
# Now lets compare cnn to transformer
# lets use only high confidence positions vs low confidence positions. (no more accuracy)

high_conf <- match_wt %>% # 23,340 positions
  filter(freq_pred >= 0.8)

low_conf <- match_wt %>% #6,858 positions
  filter(freq_pred < 0.4)
  

high <- high_conf %>%
  inner_join(sec_struc) %>%
  mutate(second_struc = ifelse(second_struc == "Coil", "Random coil", second_struc),
         second_struc = ifelse(second_struc == "Strand", "Beta strand", second_struc),
         second_struc = ifelse(second_struc == "AlphaHelix", "Alpha helix", second_struc),
         second_struc = ifelse(second_struc == "310Helix" | second_struc == "PiHelix", "Noncanon. helix", second_struc))


low <- low_conf %>%
  inner_join(sec_struc) %>%
  mutate(second_struc = ifelse(second_struc == "Coil", "Random coil", second_struc),
         second_struc = ifelse(second_struc == "Strand", "Beta strand", second_struc),
         second_struc = ifelse(second_struc == "AlphaHelix", "Alpha helix", second_struc),
         second_struc = ifelse(second_struc == "310Helix" | second_struc == "PiHelix", "Noncanon. helix", second_struc))


# finds the count of each structure per bin:
struc_counts <- high %>%
  select(c(gene, position, second_struc, group)) %>%
  group_by(group) %>%
  count(second_struc) %>%
  mutate(struc_count = n) %>%
  select(-n) %>%
  ungroup()

# now I need to add up all aa within each  bin to get bin totals and append this to the aa_counts
for_barplot <- struc_counts %>%
  group_by(group) %>%
  mutate(bin_count = sum(struc_count)) %>%
  ungroup()


#---------------

for_barplot_3 <- for_barplot %>%
  mutate(freq = struc_count/bin_count) %>%
  select(-c(struc_count, bin_count)) %>% 
  mutate(second_struc = fct_rev(fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge")))

struc_plot <- for_barplot_3 %>%
  ggplot(aes(x = freq, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.36),
    breaks = seq(0.0, 0.35, by = 0.05),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  facet_wrap(vars(fct_rev(group))) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

struc_plot

ggsave(filename = "./analysis/figures/struc_high_transf_cnn.png", plot = struc_plot, width = 9, height = 5)

for_odds <- for_barplot_3 %>%
  pivot_wider(names_from = group, values_from = freq) %>%
  mutate(odds = transformer/cnn)

plot_odds <- for_odds %>%
  ggplot(aes(x = odds, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = round(odds, 2), hjust = -0.05)) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_log10(
    name = "Odds ratio (log scale)"
  ) +
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_odds

ggsave(filename = "./analysis/figures/struc_odds_trans_cnn_high.png", plot = plot_odds, width = 9, height = 5)

#--------------------------------------------------
#now, low confidence predictions

# finds the count of each structure per bin:
struc_counts <- low %>%
  select(c(gene, position, second_struc, group)) %>%
  group_by(group) %>%
  count(second_struc) %>%
  mutate(struc_count = n) %>%
  select(-n) %>%
  ungroup()

# now I need to add up all aa within each  bin to get bin totals and append this to the aa_counts
for_barplot <- struc_counts %>%
  group_by(group) %>%
  mutate(bin_count = sum(struc_count)) %>%
  ungroup()


#---------------

for_barplot_3 <- for_barplot %>%
  mutate(freq = struc_count/bin_count) %>%
  select(-c(struc_count, bin_count)) %>% 
  mutate(second_struc = fct_rev(fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge")))

struc_plot <- for_barplot_3 %>%
  ggplot(aes(x = freq, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.36),
    breaks = seq(0.0, 0.35, by = 0.05),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  facet_wrap(vars(fct_rev(group))) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

struc_plot

ggsave(filename = "./analysis/figures/struc_low_transf_cnn.png", plot = struc_plot, width = 9, height = 5)

for_odds <- for_barplot_3 %>%
  pivot_wider(names_from = group, values_from = freq) %>%
  mutate(odds = transformer/cnn)

plot_odds <- for_odds %>%
  ggplot(aes(x = odds, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = round(odds, 2), hjust = -0.05)) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_log10(
    name = "Odds ratio (log scale)"
  ) +
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_odds

ggsave(filename = "./analysis/figures/struc_odds_trans_cnn_low.png", plot = plot_odds, width = 9, height = 5)

#----------------------------------------------------------------------------------------------
# Accuracy by secondary structure

with_struc <- match_wt %>%
  inner_join(sec_struc) %>%
  mutate(second_struc = ifelse(second_struc == "Coil", "Random coil", second_struc),
         second_struc = ifelse(second_struc == "Strand", "Beta strand", second_struc),
         second_struc = ifelse(second_struc == "AlphaHelix", "Alpha helix", second_struc),
         second_struc = ifelse(second_struc == "310Helix" | second_struc == "PiHelix", "Noncanon. helix", second_struc))

# get accuracy per secondary structure:
stats_1 <- with_struc %>%
  group_by(second_struc, group, gene) %>%
  summarise(accuracy = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)))

for_violin <- stats_1 %>%
  pivot_wider(names_from = group, values_from = accuracy)
  

violin_plot <- stats_1 %>%
  ggplot(aes(x = sec_struc, y = accuracy)) +
  geom_violin()

violin_plot












