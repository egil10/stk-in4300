
library(tidyverse)
library(HistData)

library(MASS)
df <- HistData::PearsonLee %>% as_tibble()
lm1 <- lm(child~parent+par+chl, data = df, weights = frequency)
lm2 <- stepAIC(lm1, direction = "both", trace = FALSE)

lm1 %>% summary()
lm2 %>% summary()
