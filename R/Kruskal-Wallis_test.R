# Kruskal-Wallis Test

library(tidyverse)

# create sample data
df <- data.frame(
  group = c(
    rep("A", 3), 
    rep("B", 4), 
    rep("C", 4)),
  
  value = c(
    22, 30, 42, 
    36, 40, 50, 53,
    25, 32, 45, 48)
)

df

# Kruskal-Wallis Test by hand
# 1. rank the values

df <-  df %>% 
  arrange(value) %>% 
  mutate(rank = row_number()) %>% 
  arrange(group, value)

df 

# 2. calculate the sum of ranks for each group
sum_ranks <- df %>% 
  group_by(group) %>% 
  summarise(sum_rank = sum(rank))

# convert to named vector
sum_rank_i <- as.numeric(sum_ranks$sum_rank)
names(sum_rank_i) <- sum_ranks$group

sum_rank_i

# 3. calculate the test statistic
N <- nrow(df) # total number of observations
n_i <- table(df$group) # number of observations in each group

# calculate H statistic
H <- (12 / (N * (N + 1))) * sum(sum_rank_i^2 / n_i) - 3 * (N + 1)

H

# also, calculate H statistic using matrix multiplication
# %*% is the matrix multiplication operator
# exsample: c(1, 2, 3) %*% c(4, 5, 6) = 1*4 + 2*5 + 3*6 = 32
H <- (12 / (N * (N + 1))) * (sum_rank_i %*% (sum_rank_i / n_i)) - 3 * (N + 1)

H

# or calculate H statistic using group mean ranks and median
mean_ranks <- df %>% 
  group_by(group) %>% 
  summarise(mean_rank = mean(rank))

mean_rank_i <- as.numeric(mean_ranks$mean_rank)
names(mean_rank_i) <- mean_ranks$group

H <- (12 / (N * (N + 1))) * sum(n_i * (mean_rank_i - (N + 1) / 2)^2)
H

# 4. calculate the p-value
dof <- length(unique(df$group)) - 1 # degrees of freedom
p_value <- 1 - pchisq(H, dof)
p_value

# compare with kruskal.test
kruskal.test(value ~ group, data = df)