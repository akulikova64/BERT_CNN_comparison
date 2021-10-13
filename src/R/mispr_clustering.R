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

#lets count mispredictions vs correct predictions:
set.seed(12345)

group_count <- cnn_5A %>%
  group_by(group) %>%
  count()

#sample 5081 rows from correct predictions at random

correct_5A <- cnn_5A %>%
  filter(group == "corr")
  
filtered_cor_5A <- correct_5A[sample(nrow(correct_5A), 5081), ]

mispr_5A <- cnn_5A %>%
  filter(group == "mispr")

equal_5A <- rbind(mispr_5A, filtered_cor_5A)

#---------

correct_6A <- cnn_6A %>%
  filter(group == "corr")

filtered_cor_6A <- correct_5A[sample(nrow(correct_6A), 5081), ]

mispr_6A <- cnn_6A %>%
  filter(group == "mispr")

equal_6A <- rbind(mispr_6A, filtered_cor_6A)

#----------
correct_7A <- cnn_7A %>%
  filter(group == "corr")

filtered_cor_7A <- correct_7A[sample(nrow(correct_7A), 5081), ]

mispr_7A <- cnn_7A %>%
  filter(group == "mispr")

equal_7A <- rbind(mispr_7A, filtered_cor_7A)

#------------
correct_8A <- cnn_8A %>%
  filter(group == "corr")

filtered_cor_8A <- correct_8A[sample(nrow(correct_8A), 5081), ]

mispr_8A <- cnn_8A %>%
  filter(group == "mispr")

equal_8A <- rbind(mispr_8A, filtered_cor_8A)

#---------------
correct_9A <- cnn_9A %>%
  filter(group == "corr")

filtered_cor_9A <- correct_9A[sample(nrow(correct_9A), 5081), ]

mispr_9A <- cnn_9A %>%
  filter(group == "mispr")

equal_9A <- rbind(mispr_9A, filtered_cor_9A)


# lets just look at the means:
cnn_5A_means <- equal_5A %>%
  group_by(group) %>%
  summarise(mean_prop = mean(prop_mispr))

cnn_6A_means <- equal_6A %>%
  group_by(group) %>%
  summarise(mean_prop = mean(prop_mispr))

cnn_7A_means <- equal_7A %>%
  group_by(group) %>%
  summarise(mean_prop = mean(prop_mispr))

cnn_8A_means <- equal_8A %>%
  group_by(group) %>%
  summarise(mean_prop = mean(prop_mispr))

cnn_9A_means <- equal_9A %>%
  group_by(group) %>%
  summarise(mean_prop = mean(prop_mispr))

# perhaps we need to filter out the positions where 
# there is no other amino acid but the focal in the sphere.

cnn_5A_means <- equal_5A %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(group) %>%
  summarise(mean_prop = mean(prop_mispr))

cnn_6A_means <- equal_6A %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(group) %>%
  summarise(mean_prop = mean(prop_mispr))


#grouping by gene:

cnn_7A_grouped <- equal_7A %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(gene, group) %>%
  summarise(mean_prop = mean(prop_mispr))

# get initial plot:
plot_0 <- cnn_7A_grouped %>%
  ggplot(aes(x = group, y = mean_prop)) +
  geom_violin()

plot_0

#-----------------------------------------------------
#ploting the count of mispredictions

cnn_9A_grouped <- equal_9A %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  group_by(gene, group) %>%
  summarise(mean_count = mean(num_mispr))

# get initial plot:
plot_0 <- cnn_9A_grouped %>%
  ggplot(aes(x = group, y = mean_count)) +
  geom_violin()

plot_0

#---------------------------------------------------------------------------------
#try doing binary (does it have at least one mispredicted neighbor or does it not)
# note that there are more correct predictions than mispredictions.


cnn_7A_grouped <- equal_7A %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  mutate(has_mispr_neigh = ifelse(num_mispr != 0, "has mispr. neighbor", "no mispr. neighbor")) %>%
  group_by(has_mispr_neigh, group) %>%
  count()

cnn_8A_grouped <- equal_8A %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  mutate(has_mispr_neigh = ifelse(num_mispr != 0, "has mispr. neighbor", "no mispr. neighbor")) %>%
  group_by(has_mispr_neigh, group) %>%
  count()

cnn_9A_grouped <- equal_9A %>%
  mutate(status = ifelse(num_correct == 0 & num_mispr == 0, "empty", "filled")) %>%
  filter(status != "empty") %>%
  mutate(has_mispr_neigh = ifelse(num_mispr != 0, "has mispr. neighbor", "no mispr. neighbor")) %>%
  group_by(has_mispr_neigh, group) %>%
  count()

--------------------------------------------------------------------------------
# now lets make real plots. First, with counts of mispr around residues.
  



