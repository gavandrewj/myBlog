library(tictoc)
library(data.table)
library(brms)
library(progressr)
library(tidyverse)
library(plotly)
library(furrr)
library(stringr)
library(dials)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(bayesplot)
library(truncdist)


# where we saving
file_path =  "C:/Users/gavin/Documents/blog_power_files/linear_regression/problem_three/"




regress <- brm(bf(y ~ x ,sigma ~ x), dataframes[[1]])

model <- update(regress,newdata = dataframes[[1]],formula. = bf(y ~ x ,sigma ~ x))


# summary(regress)

saveRDS(
  model,
  paste0(file_path,'regression_model.rds')
)
