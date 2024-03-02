# Rで積分計算をする
# 例として、正規分布の積分を計算してみる。

library(ggplot2)
library(plotly)
library(tidyverse)

# ---
# 1変数の積分
# ---

# まずは、正規分布の確率密度関数を定義する。
standard_normal_pdf <- function(x) {
  mean <- 0
  sd <- 1

  return(
    (1 / (sd * sqrt(2 * pi))) * exp(-((x - mean)^2) / (2 * sd^2))
         )
}

# 標準正規分布の確率密度関数をプロットする。
x <- seq(-4, 4, length.out = 100)
f_x <- standard_normal_pdf(x)

df <- data.frame(x, f_x)
df %>% head()

ggplot(df, aes(x, f_x)) +
  geom_line() +
  labs(title = "Standard Normal Distribution",
       x = "x",
       f_x = "f(x)")

# 次に、標準正規分布の確率密度関数を積分計算する。
result <- integrate(standard_normal_pdf, -Inf, Inf) # -InfからInfまでの積分
result <- integrate(standard_normal_pdf, -1.96, 1.96) # -1.96から1.96までの積分

# ---
# 2変数の積分
# ---

# まずは、2変量標準正規分布の確率密度関数を定義する。
f_xy <- function(x, y) {
  # parameters
  mean_x <- 0
  mean_y <- 0
  sd_x <- 1
  sd_y <- 1
  rho <- 0.3
  
  # z 
  z <- ((x - mean_x)^2 / sd_x^2) - (2 * rho * (x - mean_x) * (y - mean_y) / (sd_x * sd_y)) + ((y - mean_y)^2 / sd_y^2)
  return((1 / (2 * pi * sd_x * sd_y * sqrt(1 - rho^2))) * exp(-z / (2 * (1 - rho^2))))
}

# 2変量標準正規分布の確率密度関数をプロットする。
x <- rnorm(1000, mean = 0, sd = 1) * 4
y <- rnorm(1000, mean = 0, sd = 1) * 4
z <- f_xy(x, y)

df <- data.frame(x, y, z)

# plotlyで3Dプロットする
plot_ly(df, x = ~x, y = ~y, z = ~z, type = "mesh3d") %>%
  layout(scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "f(x, y)")))

# 次に、2変量標準正規分布の確率密度関数を積分計算する。
# obtain merginal distribution of x
# xに対して、そのxを固定したときのf_xyの積分を計算している。
f_x <- function(x) {
  integrate(function(y)f_xy(x=x, y=y), -1000, 1000)$value
}

# integrate f_x
# 数値的に誤差があるのは仕方ない
integrate(Vectorize(f_x), -1000, 1000, rel.tol = 1e-1)