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

trans_only <- match_wt

#join with secondary structure:
with_struc <- joined_data %>%
  inner_join(sec_struc) %>%
  mutate(second_struc = ifelse(second_struc == "Coil", "Random coil", second_struc),
         second_struc = ifelse(second_struc == "Strand", "Beta strand", second_struc),
         second_struc = ifelse(second_struc == "AlphaHelix", "Alpha helix", second_struc),
         second_struc = ifelse(second_struc == "310Helix" | second_struc == "PiHelix", "Noncanon. helix", second_struc))


# finds the count of each structure per bin:
struc_counts <- with_struc %>%
  select(c(gene, position, second_struc, div_group)) %>%
  group_by(div_group) %>%
  count(second_struc) %>%
  mutate(struc_count = n) %>%
  select(-n) %>%
  ungroup()


# now I need to add up all aa within each  bin to get bin totals and append this to the aa_counts
for_barplot <- struc_counts %>%
  group_by(div_group) %>%
  mutate(bin_count = sum(struc_count)) %>%
  ungroup()

# quick barplot with counts:
counts_table <- tibble(bin = c("n-eff > 1.5 \n KL-Div < 0.5", "KL-Div > 5", "n-eff = 1 \n KL-Div ~ 0"), 
                       count = c(1471, 1393, 136))
counts_plot <- counts_table %>%
  ggplot(aes(x = fct_relevel(bin, "n-eff > 1.5 \n KL-Div < 0.5", "KL-Div > 5", "n-eff = 1 \n KL-Div ~ 0"), y = count)) +
  geom_col(fill = "#988981", color = "#70635c", alpha = 0.8) +
  geom_text(aes(label = count, vjust = -0.25)) +
  scale_x_discrete(
    name = "") + 
  scale_y_continuous(
    name = "Count",
    limits = c(0, 1550),
    breaks = seq(0, 1400, by = 200),
    expand = c(0, 0)) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5)
  )

counts_plot

ggsave(filename = "./analysis/figures/kl_bin_counts.png", plot = counts_plot, width = 5.5, height = 4)
