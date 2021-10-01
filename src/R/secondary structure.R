# secondary structure analysis cnn vs transformer


trans_data <- read.csv(file = "./output/PSICOV_BERT_predictions.csv", header=TRUE, sep=",")
cnn_data_all <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
sec_struc <- read.csv(file = "./output/second_struc.csv", header=TRUE, sep=",")

cnn_data <- cnn_data_all %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) %>%
  rename(aa_pred_cnn = aa_predicted,
         aa_wt_cnn = aa_wt,
         freq_pred_cnn = freq_predicted,
         freq_wt_cnn = freq_wt)

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

#ggsave(filename = "./analysis/figures/struc_high_low_transf.png", plot = struc_mispr_plot, width = 9, height = 5)

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

#ggsave(filename = "./analysis/figures/struc_odds_trans.png", plot = plot_odds, width = 9, height = 5)

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

#ggsave(filename = "./analysis/figures/struc_high_transf_cnn.png", plot = struc_plot, width = 9, height = 5)

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

#ggsave(filename = "./analysis/figures/struc_odds_trans_cnn_high.png", plot = plot_odds, width = 9, height = 5)

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

#ggsave(filename = "./analysis/figures/struc_low_transf_cnn.png", plot = struc_plot, width = 9, height = 5)

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

#ggsave(filename = "./analysis/figures/struc_odds_trans_cnn_low.png", plot = plot_odds, width = 9, height = 5)

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
  pivot_wider(names_from = group, values_from = accuracy) %>%
  mutate(change = cnn - transformer) 
  
for_violin <- for_violin %>%  
  mutate(second_struc = fct_rev(fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge")))

fills = c("#bf8040", "#ac5396", "#70adc2", "#748f3d", "#cc5933", "#7070c2")
  
#fills = c("#ac5396","#bf8040", "#70adc2","#748f3d", "#cc5933", "#7070c2")
#fills = c("#ff00ff","#ff0000", "#00ff00","#0000ff", "#ffff00", "#000000")
colors = c("#470e39", "#522c05", "#183f4d","#29380b", "#6b1f06", "#0c0c5c")


violin_plot <- for_violin %>%
  ggplot(aes(x = fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge"), 
             y = change, 
             fill = second_struc, 
             color = second_struc
             )) +
  geom_hline(yintercept = 0, color = "#666666", linetype="dashed") +
  theme_cowplot(16) +
  geom_violin() +
  geom_sina() + 
  theme(plot.title = element_text(hjust = 0, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position="none",
        axis.text.x = element_text(
          size = 16,
          color = "black",
          angle = 45, 
          hjust = 1,
          vjust = 1)
        ) + 
  #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  scale_y_continuous(
    name = "Change in accuracy \n from transformer to CNN",
    limits = c(-1.0, 1.0),
    breaks = seq(from = -1.0, to = 1.0, by = 0.2),
    expand = c(0, 0)) +
  labs(fill = "") +
  scale_x_discrete(
    name = "",
    position = "bottom") +
  scale_fill_manual(values = fills) +
  scale_color_manual(values = colors)


violin_plot

ggsave(filename = "./analysis/figures/violin_trans_cnn_struc.png", plot = violin_plot, width = 8, height = 5)

for_bar <- for_violin %>%
  filter(change != 0) %>%
  mutate(change_bin = ifelse(change < 0, "transformer better", "cnn better")) %>%
  group_by(change_bin, second_struc) %>%
  count() %>%
  mutate(count = ifelse(change_bin == "transformer better", -n, n)) %>%
  mutate(text_label = abs(count)) 

bar_chart <- for_bar %>%
  ggplot(aes(y = fct_rev(fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge")), x = count, fill = change_bin)) +
  geom_col(alpha = 0.7) +
  geom_vline(xintercept = 0, color = "#666666", linetype="dashed") +
  scale_fill_manual(values = c("#216610", "#d16d0f")) +
  theme_cowplot(16) +
  scale_x_continuous(
    name = "Protein Count",
    limits = c(-70, 90),
    breaks = seq(-70, 80, by = 10),
    labels = c("70", "60", "50", "40", "30", "20", "10", "0", "10", "20", "30", "40", "50", "60", "70", "80"),
    expand = c(0, 0.03)) +
  labs(fill = "") +
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.03)) + 
  theme(plot.title = element_text(hjust = 0, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position="top") +
  geom_text(aes(label = text_label), 
            vjust = 0.5, 
            hjust = -0.5,
            #position = position_dodge(width = 0.8),
            color = "black") 

bar_chart

ggsave(filename = "./analysis/figures/bar_trans_cnn_struc.png", plot = bar_chart, width = 9, height = 5)

#----------------------------------------------------------------------------
# another violin plot, only with average confidence per protein this time.


stats_2 <- with_struc %>%
  group_by(second_struc, group, gene) %>%
  summarise(mean_conf = mean(freq_pred))

for_violin <- stats_2 %>%
  pivot_wider(names_from = group, values_from = mean_conf) %>%
  mutate(change = cnn - transformer) 

for_violin <- for_violin %>%  
  mutate(second_struc = fct_rev(fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge")))

violin_plot2 <- for_violin %>%
  ggplot(aes(x = fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge"), 
             y = change, 
             fill = second_struc, 
             color = second_struc
  )) +
  geom_hline(yintercept = 0, color = "#666666", linetype="dashed") +
  theme_cowplot(16) +
  geom_violin() +
  geom_sina() + 
  theme(plot.title = element_text(hjust = 0, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position="none",
        axis.text.x = element_text(
          size = 16,
          color = "black",
          angle = 45, 
          hjust = 1,
          vjust = 1)
  ) + 
  #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  scale_y_continuous(
    name = "Change in mean confidence \n from transformer to CNN",
    limits = c(-1.0, 1.0),
    breaks = seq(from = -1.0, to = 1.0, by = 0.2),
    expand = c(0, 0)) +
  labs(fill = "") +
  scale_x_discrete(
    name = "",
    position = "bottom") +
  scale_fill_manual(values = fills) +
  scale_color_manual(values = colors)


violin_plot2

ggsave(filename = "./analysis/figures/violin_trans_cnn_struc_conf.png", plot = violin_plot2, width = 8, height = 5)

# --------------------------------------------------------------------------
# mispredictions: both cnn and transformer

mispr <- with_struc %>%
  filter(match_predict_wt == FALSE)

# finds the count of each structure per bin:
struc_counts <- mispr %>%
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

ggsave(filename = "./analysis/figures/trans_cnn_struc_mispr.png", plot = struc_mispr_plot, width = 8, height = 5)

for_odds <- for_barplot_2 %>%
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
    expand = c(0.03, 0.04)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_odds

ggsave(filename = "./analysis/figures/struc_odds_trans_cnn_mispr.png", plot = plot_odds, width = 8.5, height = 7)

#-----------------------------------------------------------------------
# lets see if we can get amino acid distributions for mispredictions

joined2 <- rbind(cnn_data_unjoined, transf_data_unjoined)

match_wt <- joined2 %>%
  mutate(match_predict_wt = aa_pred == aa_wt) %>%
  select(c(gene, position, aa_wt, match_predict_wt, group))

mispr <- match_wt %>%
  filter(match_predict_wt == FALSE)

aa_counts <- mispr %>%
  select(c(gene, position, aa_wt, group)) %>%
  group_by(group) %>%
  count(aa_wt) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

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
    return("small polar")
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

with_class <- aa_counts %>%
  mutate(class = map_chr(aa_wt, calc_class))

for_barplot <- with_class %>%
  group_by(group) %>%
  mutate(group_count = sum(aa_count)) %>%
  ungroup()

for_barplot_2 <- for_barplot %>%
  mutate(freq = aa_count/group_count) %>%
  select(-c(aa_count, group_count)) %>%
  mutate(aa_wt = fct_rev(fct_relevel(aa_wt, "K", "E", "Q", "S", "N", "R", "D", "A", "T", "V", "I", "F", "L", "H", "M", "Y", "C", "W", "G", "P"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small polar", "negative", "positive", "aromatic", "unique"))

order <- for_barplot_2 %>%
  filter(group == "cnn")


fills <- c("#990008", "#0a2575", "#b35900", "#1a6600", "#5c0679", "#9e9e2e")

plot_a <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = aa_wt, fill = class)) +
  geom_col(alpha = 0.7) +
  facet_wrap(vars(group)) +
  scale_fill_manual(
    values = fills) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.12),
    breaks = seq(0.0, 0.12, by = 0.02),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Wild type amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(14) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "right")

plot_a

ggsave(filename = "./analysis/figures/cnn_trans_mispr_aa.png", plot = plot_a, width = 11, height = 8)  

#----------- odds ratio


for_odds <- for_barplot_2 %>%
  pivot_wider(names_from = group, values_from = freq) %>%
  mutate(odds = transformer/cnn)

plot_odds <- for_odds %>%
  ggplot(aes(x = odds, y = fct_reorder(aa_wt, odds), fill = class)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = round(odds, 2), hjust = -0.3)) +
  scale_fill_manual(
    values = fills) +
  scale_x_log10(
    name = "Odds ratio (log scale)"
  ) +
  geom_hline(yintercept = 12.5, color = "#666666", linetype="dashed") +
  scale_y_discrete(
    name = "",
    expand = c(0.03, 0.04)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "right")

plot_odds

ggsave(filename = "./analysis/figures/aa_odds_trans_cnn_mispr.png", plot = plot_odds, width = 9, height = 7)
