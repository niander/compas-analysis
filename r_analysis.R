library(plyr)
library(tidyverse)
library(purrrlyr)

compas_data <- read.csv("./compas-scores-two-years.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))

compas_data_v <- read.csv("./compas-scores-two-years-violent.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))

cox_data <- read.csv("./cox-parsed.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))

#library(causaleffect)
library(pcalg)
library(mediation)
library(igraph)

g <- graph_from_literal(S:A:R -+ NP:CR, NP -+ CR, S:A:R:NP:CR -+ SCP)

causal.effect(y = "SCP", x = "R", G = g)

cpss <- compas_data %>% select(id = id, S = sex, A = age_cat, R = race, 
                               NP = priors_count, CR = c_charge_degree, SCP = decile_score) %>%
  mutate(SCP = as.factor(SCP))

logical_gaps <- cpss %>%
  select(-id, -NP) %>% {
    matrix(FALSE, nrow = length(colnames(.)), ncol = length(colnames(.)))
  } %>% {
    colnames(.) <- colnames(select(cpss, -id, -NP))
    rownames(.) <- colnames(.)
    .
  }

logical_gaps["S", "A"] <- logical_gaps["A", "S"] <- T
logical_gaps["A", "R"] <- logical_gaps["R", "A"] <- T
logical_gaps["S", "R"] <- logical_gaps["R", "S"] <- T

skel.cpss <- cpss %>%
  select(-id, -NP) %>% {
  skeleton(list(dm = mutate_all(., function(a) as.integer(a) - 1), nlev = as.numeric(dmap(., nlevels)),
                adaptDF = FALSE), 
           indepTest = disCItest,
           labels = colnames(.), alpha = 0.01,
           fixedGaps = logical_gaps)
  }

plot(skel.cpss, main = "Estimated CPDAG")

###
# Demo pcalg
data("gmG")
suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x))
varNames <- gmG8$g@nodes
skel.gmG8 <- skeleton(suffStat, indepTest = gaussCItest,
                      labels = varNames, alpha = 0.01)

data("gmD")
V <- colnames(gmD$x)
## define sufficient statistics
suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)
## estimate CPDAG
pc.D <- pc(suffStat,
           ## independence test: G^2 statistic
           indepTest = disCItest, alpha = 0.01, labels = V, verbose = TRUE)
