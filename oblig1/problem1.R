
library(tidyverse)
library(HistData)

df <- HistData::PearsonLee %>% as_tibble()
df %>% summary()

?PearsonLee

# a)

p <- ggplot(df, aes(x = child, weight = frequency)) +
  geom_histogram(bins = 3,
                 fill = "#D2B48C",
                 color = "black") +
  ggtitle("Distribution of child height") +
  coord_flip() +
  theme_minimal() 

ggsave(
  "plots/child_height_bad.pdf",
  plot = p,
  width = 12,
  height = 6,
  units = "in"
)


q <- ggplot(df, aes(x = par, weight = frequency, fill = par)) +
  geom_bar() +
  coord_polar(theta = "y") +   
  labs(title = "Parent sex distribution", 
       x = NULL, 
       y = NULL) +
  theme_minimal()

ggsave(
  "plots/parent_sex_distribution.pdf",
  plot = q,
  width = 12,
  height = 6,
  units = "in"
)

  

# c)

p2 <- ggplot(df, aes(x = child, weight = frequency)) +
  geom_histogram(bins = 29,
                 fill = "#D2B48C",
                 color = "black") +
  ggtitle("Distribution of child height") +
  xlab("Height of children in inches") +
  ylab("Count") +
  theme_minimal() 

ggsave(
  "plots/child_height_good.pdf",
  plot = p2,
  width = 12,
  height = 6,
  units = "in"
)

# -------------------------------------------------------------------------

q2 <- ggplot(df, aes(
  x = parent,
  weight = frequency,
  fill = par,
  colour = par
)) +
  geom_density(alpha = 0.35, adjust = 1.1) +
  scale_fill_manual(values = c(Father = "#4C78A8", Mother = "#F56598")) +
  scale_colour_manual(values = c(Father = "#4C78A8", Mother = "#F56598")) +
  labs(
    title = "Parent height distribution by sex",
    x = "Parent height in inches",
    y = "Density",
    fill = NULL,
    colour = NULL
  ) +
  theme_minimal()

ggsave(
  "plots/parent_sex_distribution_good.pdf",
  plot = q2,
  width = 12,
  height = 6,
  units = "in"
)
