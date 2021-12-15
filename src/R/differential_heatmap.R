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
  labs(fill = "Change in frequency \n from CNN model") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = font_size),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_cnn_aa




















