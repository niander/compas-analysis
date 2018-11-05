library(plyr)
library(tidyverse)

compas_data <- read.csv("./compas-scores-two-years.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))

compas_data_v <- read.csv("./compas-scores-two-years-violent.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))

cox_data <- read.csv("./cox-parsed.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))
