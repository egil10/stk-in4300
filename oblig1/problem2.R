
library(tidyverse)
library(HistData)
library(MASS)

df <- HistData::PearsonLee %>% as_tibble()

lm1 <- lm(child~parent+par+chl, data = df, weights = frequency)
lm1

lm2 <- stepAIC(lm1, direction = "both", trace = FALSE)
lm2