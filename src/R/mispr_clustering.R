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

# get initial plot:
plot_0 <- cnn_grouped %>%
  ggplot(aes(x = group, y = mean_prop)) +
  geom_violin() +
  geom_sina() +
  facet_wrap(vars(fct_relevel(radius, "5A", "6A", "7A", "8A", "9A", "10A")))

plot_0

#-----------------------------------------------------
#plotting the count of mispredictions

for_count <- joined_equal %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(gene, group, radius) %>%
  summarise(mean_count = mean(num_mispr))

# get initial plot:
plot_0 <- for_count %>%
  ggplot(aes(x = group, y = mean_count)) +
  geom_violin() +
  geom_sina() +
  facet_wrap(vars(fct_relevel(radius, "5A", "6A", "7A", "8A", "9A", "10A")))

plot_0

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
    name = "Mean number of residues that \nhave at least one mispredicted site",
    limits = c(0, 32),
    breaks = seq(from = 0, to = 32, by = 2),
    expand = c(0, 0)) +
  labs(fill = "focal residue") +
  scale_x_discrete(
    name = "Radius") +
  scale_fill_manual(values = c("#44b8a6", "#d94e9f")) +
  theme_cowplot(16) +
  theme(plot.title = element_text(hjust = 0, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "right")
plot



