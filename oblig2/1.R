
library(tidyverse)
df = read.csv("data/qsar_aquatic_toxicity.csv",
         sep = ";",
         header = F) %>% 
  as_tibble()

p <- df %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point() 

ggsave(filename = "plots/randomplot.pdf", 
       plot = p, width = 12, height = 6)
