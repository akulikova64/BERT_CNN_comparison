library(tidyverse)
library(cowplot)
library(data.table)
options(scipen = 999)

# get the 

data <- read.csv(file = "./output/cnn_wt_max_freq_3ogo.csv", header=TRUE, sep=",")

data2 <- data %>%
  filter(chain == 'G') %>%
  pivot_wider(names_from = group, values_from = c(aa, freq))

seq = data2$aa_wt
paste(seq, collapse = "")


write.csv(data2, "./output/cnn_3ogo.csv")

data_bert <- read.csv(file = "./output/bert_prediction_3ogo.csv", header=TRUE, sep=",")


#ploting data2 from above:

plot_a <- data2 %>%
  ggplot(aes(x = position, y = freq_wt)) +
  geom_line()

plot_a


#---------------cleaning data for Jupyter notebook-----------------------------

#cleaning the raw output:
data <- read.csv(file = "./output/3ogo_final_tot.csv", header=TRUE, sep=",")

count<- data %>%
  group_by(chain_id) %>%
  count()

clean_data <- data %>%
  filter(chain_id == "B") %>%
  select(c(prALA:prVAL)) %>%
  mutate_if(is.numeric, round, digits = 6)

names(clean_data) <-  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

write.csv(clean_data, "./output/clean_cnn_3ogo.csv")






