#p_load(causaleffect)
p_load(plyr)
p_load(tidyverse)
#p_load(tidymodels)
p_load(purrrlyr)

compas_data <- read.csv("./compas-scores-two-years.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))

compas_data_v <- read.csv("./compas-scores-two-years-violent.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))

cox_data <- read.csv("./cox-parsed.csv") %>%
  as_tibble() %>%
  select(-ends_with(".1"))

##################
p_load(gRain)
p_load(bnlearn)


dag <- model2network("[S][A][R][NP|S:A:R][CR|S:A:R:NP][SCP|S:A:R:NP:CR]")

cmps <- compas_data %>% 
  select(id = id, S = sex, A = age_cat, R = race, 
         NP = priors_count, CR = c_charge_degree, SCP = decile_score) %>%
  mutate(NP = as.numeric(NP), SCP = as.factor(SCP >= 5))

bn <- bn.fit(dag, select(cmps, -id) %>% discretize())

gr.bn <- compile(as.grain(bn))

gr.bn.white <- setEvidence(gr.bn, nodes = "R", states = "Caucasian")
gr.bn.black <- setEvidence(gr.bn, nodes = "R", states = "African-American")

E_y_xz <- querygrain(gr.bn.white, c("SCP", "NP", "CR", "S", "A"), "conditional") %>%
  ar_slice(list(SCP = "TRUE"))
P_z1_x <- querygrain(gr.bn.white, c("CR", "NP", "S", "A"), "conditional")
P_z2_x <- querygrain(gr.bn.white, c("NP", "S", "A"), "conditional")

E_y_x1z <- querygrain(gr.bn.black, c("SCP", "NP", "CR", "S", "A"), "conditional") %>%
  ar_slice(list(SCP = "TRUE"))
P_z1_x1 <- querygrain(gr.bn.black, c("CR", "NP", "S", "A"), "conditional")
P_z2_x1 <- querygrain(gr.bn.black, c("NP", "S", "A"), "conditional")

P_c <- querygrain(gr.bn, c("S", "A"), "joint")

ADE_x1_x = sum((E_y_xz %a-% E_y_x1z) %a*% P_z1_x1  %a*% P_z2_x1 %a*% P_c)
ADE_x_x1 = sum((E_y_x1z %a-% E_y_xz) %a*% P_z_x)

#P_z_x1 <- querygrain(gr.bn.black, c("NP", "CR", "S", "A"), "joint")
#P_z_x <- querygrain(gr.bn.white, c("NP", "CR", "S", "A"), "joint")

##################
p_load(mediation)
p_unload(MASS)

cmps.data <- select(cmps, -id) %>%
  mutate(SCP = if_else(SCP == "TRUE", 1, 0)) %>%
  filter(R %in% c("African-American", "Caucasian")) %>%
  mutate(R = fct_drop(R))


med.fit <- lm(NP ~ S + A + R, cmps.data)
out.fit <- glm(SCP ~ S + A + CR + R * NP, cmps.data, family = binomial("probit"))
med.out <- mediate(model.m = med.fit, model.y = out.fit, 
                   treat = "R", outcome = "SCP", mediator = "NP",
                   covariates = c("S", "A", "CR"),
                   control.value = "Caucasian", treat.value = "African-American",
                   boot = T, sims = 100)

summary(med.out)
plot(med.out)

test.TMint(med.out, conf.level = .95)

med.out <- mediations(list("cmps" = cmps.data), 
                      treatment = c("R"), 
                      mediators = c("NP", "CR"), 
                      outcome = c("SCP"), 
                      covariates = "S + A",
                      families = c("gaussian", "binomial"))

summary(med.out)

################
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital + 
                nonwhite + educ + income, data = jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + 
                occp + marital + nonwhite + educ + income, data = jobs)
out.1 <- mediate(model.m, model.y, sims = 1000, boot = TRUE, treat = "treat", mediator = "job_seek")

summary(out.1)

#################
p_load(mma)

###################
p_load(igraph)
g <- graph_from_literal(S:A:R -+ NP:CR, NP -+ CR, S:A:R:NP:CR -+ SCP)

causal.effect(y = "SCP", x = "R", G = g)

cmps <- compas_data %>% 
  select(id = id, S = sex, A = age_cat, R = race, 
         NP = priors_count, CR = c_charge_degree, SCP = decile_score) %>%
  mutate(SCP = as.factor(SCP))



#######################
logical_gaps <- cmps %>%
  select(-id, -NP) %>% {
    matrix(FALSE, nrow = length(colnames(.)), ncol = length(colnames(.)))
  } %>% {
    colnames(.) <- colnames(select(cmps, -id, -NP))
    rownames(.) <- colnames(.)
    .
  }

logical_gaps["S", "A"] <- logical_gaps["A", "S"] <- T
logical_gaps["A", "R"] <- logical_gaps["R", "A"] <- T
logical_gaps["S", "R"] <- logical_gaps["R", "S"] <- T

skel.cmps <- cmps %>%
  select(-id, -NP) %>% {
  skeleton(list(dm = mutate_all(., function(a) as.integer(a) - 1), nlev = as.numeric(dmap(., nlevels)),
                adaptDF = FALSE), 
           indepTest = disCItest,
           labels = colnames(.), alpha = 0.01,
           fixedGaps = logical_gaps)
  }

plot(skel.cmps, main = "Estimated CPDAG")

###
#p_load(pcalg)
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
