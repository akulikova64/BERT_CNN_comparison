library(tidyverse)
library(cowplot)
library(data.table)
library(ggforce)
library(colorspace)
library(gridExtra)
library(grid)
options(scipen = 999)


# do the mispredictions cluster in proteins for both the cnn and trans?

cnn_5A <- read.csv(file = "./output/sphere_densities_5A.csv", header=TRUE, sep=",")
cnn_6A <- read.csv(file = "./output/sphere_densities_6A.csv", header=TRUE, sep=",")
cnn_7A <- read.csv(file = "./output/sphere_densities_7A.csv", header=TRUE, sep=",")
cnn_8A <- read.csv(file = "./output/sphere_densities_8A.csv", header=TRUE, sep=",")
cnn_9A <- read.csv(file = "./output/sphere_densities_9A.csv", header=TRUE, sep=",")
cnn_10A <- read.csv(file = "./output/sphere_densities_10A.csv", header=TRUE, sep=",")

equal_5A <- read.csv(file = "./output/equal_5A.csv", header=TRUE, sep=",")

counts <- equal_5A %>%
  group_by(group) %>%
  count()

positions <- equal_5A %>%
  select(c(gene, position, group))

equal_5A <- equal_5A %>%
  mutate(radius = "5A")
equal_6A <- inner_join(cnn_6A, positions) %>%
  mutate(radius = "6A")
equal_7A <- inner_join(cnn_7A, positions) %>%
  mutate(radius = "7A")
equal_8A <- inner_join(cnn_8A, positions) %>%
  mutate(radius = "8A")
equal_9A <- inner_join(cnn_9A, positions) %>%
  mutate(radius = "9A")
equal_10A <- inner_join(cnn_9A, positions) %>%
  mutate(radius = "10A")

joined_equal = rbind(equal_5A, equal_6A, equal_7A, equal_8A, equal_9A, equal_10A)


# lets just look at the means:
prop_means <- joined_equal %>%
  group_by(group, radius) %>%
  summarise(mean_prop = mean(prop_mispr))

mispr_count_means <- joined_equal %>%
  group_by(group, radius) %>%
  summarise(mean_mispr_count = mean(num_mispr))

# perhaps we need to filter out the positions where 
# there is no other amino acid but the focal in the sphere.

prop_means <- joined_equal %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(group, radius) %>%
  summarise(mean_prop = mean(prop_mispr))

mispr_count_means <- joined_equal %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(group, radius) %>%
  summarise(mean_mispr_count = mean(num_mispr))

#grouping by gene:

cnn_grouped <- joined_equal %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(gene, group, radius) %>%
  summarise(mean_prop = mean(prop_mispr))

stat_data_1 <- cnn_grouped %>%
  select(-gene) %>%
  group_by(radius, group) %>%
  summarise(estimate = mean(mean_prop),
            std_error = sd(mean_prop)/sqrt(length(mean_prop))) %>%
  na.omit()

cnn_grouped2 <- cnn_grouped %>%
  group_by(radius, group) %>%
  summarise(text_label = round(mean(mean_prop), 3)) %>%
  mutate(mean_prop = text_label)

# get proportion plot:
plot_0 <- cnn_grouped %>%
  ggplot(aes(x = group, y = mean_prop, fill = group, color = group)) +
  geom_violin(alpha = 0.8) +
  geom_sina(alpha = 0.4) +
  geom_pointrange(data = stat_data_1, 
                  aes(x = group,
                      y = estimate,
                      ymin = estimate - 1.96*std_error,
                      ymax = estimate + 1.96*std_error),
                  color = "black", 
                  size = 0.7) +
  scale_y_continuous(
    name = "Mean proportion of mispredicted sites \nout of all local neighbors (per protein)",
    limits = c(0.00, 0.45),
    breaks = seq(from = 0.00, to = 0.45, by = 0.05),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "",
    labels = c("correct pred.", "mispredictions")) +
  facet_wrap(vars(fct_relevel(radius, "5A", "6A", "7A", "8A", "9A", "10A"))) +
  scale_fill_manual(values = c("#44b8a6", "#d94e9f"),
                    labels = c("correct pred.", "mispred.")) +
  scale_color_manual(values = c("#26665c", "#7d2d5c")) +
  geom_text(data = cnn_grouped2,
            aes(label = text_label, x = group, y = mean_prop),
            vjust = -6,
            hjust = 0.5,
            color = "black",
            position = position_dodge()) +
  theme_cowplot(16) +
  theme(plot.title = element_text(hjust = 0, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        panel.grid.minor.y = element_line(color = "grey92", size=0.5),
        legend.position = "none")
plot_0

ggsave(filename = "./analysis/figures/mean_prop_mispr.png", plot = plot_0, width = 11, height = 5)



#-----------------------------------------------------
#plotting the count of mispredictions

for_count <- joined_equal %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(gene, group, radius) %>%
  summarise(mean_count = mean(num_mispr))

stat_data_1 <- for_count %>%
  select(-gene) %>%
  group_by(radius, group) %>%
  summarise(estimate = mean(mean_count),
            std_error = sd(mean_count)/sqrt(length(mean_count))) %>%
  na.omit()

for_count2 <- for_count %>%
  group_by(radius, group) %>%
  summarise(text_label = round(mean(mean_count), 3)) %>%
  mutate(mean_count = text_label)

# get count plot:
plot_1 <- for_count %>%
  ggplot(aes(x = group, y = mean_count, fill = group, color = group)) +
  geom_violin(alpha = 0.8) +
  geom_sina(alpha = 0.4) +
  geom_pointrange(data = stat_data_1, 
                  aes(x = group,
                      y = estimate,
                      ymin = estimate - 1.96*std_error,
                      ymax = estimate + 1.96*std_error),
                  color = "black", 
                  size = 0.7) +
  scale_y_continuous(
    name = "Mean count of mispredicted sites \nout of all local neighbors (per protein)",
    limits = c(0, 5),
    breaks = seq(from = 0, to = 5, by = 1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "",
    labels = c("correct pred.", "mispredictions")) +
  facet_wrap(vars(fct_relevel(radius, "5A", "6A", "7A", "8A", "9A", "10A"))) +
  scale_fill_manual(values = c("#44b8a6", "#d94e9f"),
                    labels = c("correct pred.", "mispred.")) +
  scale_color_manual(values = c("#26665c", "#7d2d5c")) +
  geom_text(data = for_count2,
            aes(label = text_label, x = group, y = mean_count),
            vjust = -6,
            hjust = 0.5,
            color = "black",
            position = position_dodge()) +
  theme_cowplot(16) +
  theme(plot.title = element_text(hjust = 0, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        panel.grid.minor.y = element_line(color = "grey92", size=0.5),
        legend.position = "none")

plot_1

ggsave(filename = "./analysis/figures/mean_count_mispr.png", plot = plot_1, width = 11, height = 5)


#---------------------------------------------------------------------------------
#try doing binary (does it have at least one mispredicted neighbor or does it not)
# note that there are more correct predictions than mispredictions.

joined_counts <- joined_equal %>%
  group_by(group, radius) %>%
  count()

cnn_binary <- joined_equal %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  mutate(has_mispr_neigh = ifelse(num_mispr != 0, "has mispr. neighbor", "no mispr. neighbor")) %>%
  group_by(has_mispr_neigh, group, radius, gene) %>%
  count() %>%
  filter(has_mispr_neigh == "has mispr. neighbor") %>%
  group_by(has_mispr_neigh, group, radius) %>%
  summarise(mean_per_proteins = mean(n)) # mean number of positions (protein) that have at least one neighbor

plot <- cnn_binary %>%
  ggplot(aes(x = fct_relevel(radius, "5A", "6A", "7A", "8A", "9A", "10A"), y = mean_per_proteins, fill = group)) +
  geom_col(position = position_dodge()) +
  scale_y_continuous(
    name = "Mean number of residues per protein that \nhave at least one mispredicted site",
    limits = c(0, 32),
    breaks = seq(from = 0, to = 32, by = 2),
    expand = c(0, 0)) +
  labs(fill = "Focal residue") +
  scale_x_discrete(
    name = "Radius") +
  scale_fill_manual(values = c("#44b8a6", "#d94e9f"),
                    labels = c("correct pred.", "mispred.")) +
  theme_cowplot(16) +
  theme(plot.title = element_text(hjust = 0, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "right")
plot

ggsave(filename = "./analysis/figures/at_least_one_neighbor.png", plot = plot, width = 11, height = 5)



