library(ggmcmc)
library(ggthemes)
library(ggridges)
###################
# data processing #
###################
options(mc.cores = parallel::detectCores())
data <- read.csv("ipldata.csv", header=TRUE)
#exclude data of 2008, 2009 and 2010 due to different structure of playoffs
ipl.data <- data[!(data$Year==2008|data$Year==2009|data$Year==2010|
                     data$Year==2021),] 
head(ipl.data)
#ipl.data$End_Position <- as.factor(ipl.data$End_Position)
ipl.data$League_Position <- as.factor(ipl.data$League_Position)
ipltest <- ipl.data[,c(13,2,10,11)]
ipl.playoff <- ipltest
names(ipl.playoff) <- c("champion", "seed", "nrr", "last5")
head(ipl.playoff) #cleaned dataset 
ipl.2021 <- data[data$Year==2021,c(2,10,11,12,13)]
colnames(ipl.2021) <- c("seed", "nrr", "last5", "result","champion")
head(ipl.2021)
#########################################################################
data1 <- ipl.data[ipl.data$League_Position==1,] 
data2 <- ipl.data[ipl.data$League_Position==2,] 
data3 <- ipl.data[ipl.data$League_Position==3,]
data4 <- ipl.data[ipl.data$League_Position==4,]

################
## JAGS Model ##
################
library(rjags)
library(reshape2)
library(ggplot2)
library(plyr)
library(Rmisc)

jags_model = "
model {
  for (i in 1:n){
    y[i] ~ dcat(pos[])
  }
  pos ~ ddirich(a)
}"

##############################################
##for team placed first in the league stages##
##############################################
dat1 = list(n = length(data1$End_Position), y = data1$End_Position, a=rep(1,3))
m1 = jags.model(textConnection(jags_model), data=dat1, n.chains = 3)
res1 =coda.samples(m1, 'pos', 1000)
summary(res1)
plot(res1)

# plot to compare posterior of pi
dr1 = ldply(res1, function(x) {as.data.frame(x)})
m1 = melt(dr1[,grep("pos", names(dr1))], variable.name="estimate", 
          value.name="sample")
(ggplot(m1, aes(x=sample, color=estimate))+geom_density())

###############################################
##for team placed second in the league stages##
###############################################
dat2 = list(n = length(data2$End_Position), y = data2$End_Position, a=rep(1,3))
m2 = jags.model(textConnection(jags_model), data=dat2, n.chains = 3)
res2 =coda.samples(m2, 'pos', 1000)
summary(res2)
plot(res2)

# plot to compare posterior of pi
dr2 = ldply(res2, function(x) {as.data.frame(x)})
m2 = melt(dr2[,grep("pos", names(dr2))], variable.name="estimate", 
          value.name="sample")
(ggplot(m2, aes(x=sample, color=estimate))+geom_density())


###############################################
##for team placed third in the league stages##
###############################################
dat3 = list(n = length(data3$End_Position), y = data3$End_Position, a=rep(1,4))
m3 = jags.model(textConnection(jags_model), data=dat3, n.chains = 3)
res3 =coda.samples(m3, 'pos', 1000)
summary(res3)
plot(res3)

# plot to compare posterior of pi
dr3 = ldply(res3, function(x) {as.data.frame(x)})
m3 = melt(dr3[,grep("pos", names(dr3))], variable.name="estimate", 
          value.name="sample")
(ggplot(m3, aes(x=sample, color=estimate))+geom_density())


###############################################
##for team placed forth in the league stages##
###############################################
dat4 = list(n = length(data4$End_Position), y = data4$End_Position, a=rep(1,4))
m4 = jags.model(textConnection(jags_model), data=dat4, n.chains = 3)
res4 =coda.samples(m4, 'pos', 1000)
summary(res4)
plot(res4)

# plot to compare posterior of pi
dr4 = ldply(res4, function(x) {as.data.frame(x)})
m4 = melt(dr4[,grep("pos", names(dr4))], variable.name="estimate", 
          value.name="sample")
(ggplot(m4, aes(x=sample, color=estimate))+geom_density())


################################################
# Multiple plots in the same figure
################################################

library(bayesplot)
library(mcmcOutput)

pos1 <- mcmcOutput(res1); pos2 <- mcmcOutput(res2)
pos3 <- mcmcOutput(res3); pos4 <- mcmcOutput(res4)

data.seed.1 <- cbind(pos1$pos, rep(0,nrow(pos1)))
data.seed.2 <- cbind(pos2$pos, rep(0,nrow(pos2)))
data.seed.3 <- cbind(pos3$pos)
data.seed.4 <- cbind(pos4$pos)

#colnames(data.seed.1) <- colnames(data.seed.2) <- paste0("pos[", 1:4, "]")
#colnames(data.seed.3) <- colnames(data.seed.4) <- paste0("pos[", 1:4, "]")

colnames(data.seed.1) <- colnames(data.seed.2) <- 1:4
colnames(data.seed.3) <- colnames(data.seed.4) <- 1:4


combined <- rbind(mcmc_intervals_data(data.seed.1, prob = 0.8), 
                  mcmc_intervals_data(data.seed.2, prob = 0.8),
                  mcmc_intervals_data(data.seed.3, prob = 0.8), 
                  mcmc_intervals_data(data.seed.4, prob = 0.8))

combined$seed <- rep(c("Seed 1", "Seed 2",
                       "Seed 3", "Seed 4"), each = 4)

pos <- position_nudge(y = ifelse(combined$seed == "Seed 4", 0.4, 
                              ifelse(combined$seed == "Seed 3",0.3,
                                  ifelse(combined$seed == "Seed 2",0.2,0))))
ggplot(combined, aes(x = m, y = parameter, color = seed))+
  geom_linerange(aes(xmin = l, xmax = h), position = pos, size=2)+
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos)+
  geom_point(position = pos, color="black")+
  ylab("result") + xlab("probability") + 
  ggtitle("result and probabilities by seed")

################
## STAN Model ##
################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)
categorical.stan="
data {
  int<lower = 1> N;
  int<lower = 1> m;
  int<lower = 1, upper = m> y[N];
}
parameters {
  simplex[m] pi;
}
model {
  target += dirichlet_lpdf(pi | rep_vector(2, m));
  for(n in 1:N)
    target += categorical_lpmf(y[n] | pi);
}
generated quantities{
  int Y[N];
  for(n in 1:N)
    Y[n] = categorical_rng(pi);
}
"

##############################################
##for team placed first in the league stages##
##############################################
m = stan_model(model_code=categorical.stan)
dat1 = list(N = length(data1$End_Position), m = 3, y = data1$End_Position)
r1 = sampling(m, dat1, c("pi"), iter=10000,
              control = list(max_treedepth=50, adapt_delta = 0.99))
r1 #p1 = 0.29; pi2 = 0.41; pi3=0.29
plot(r1, main="Team placed 1 after league stage")
traceplot(r1, pars = c("pi")) #proper mixing

#posterior plots
posterior1 <- as.array(r1)
mcmc_hist(posterior1, pars = c("pi[1]","pi[2]","pi[3]"))
mcmc_hist_by_chain(posterior1, pars = c("pi[1]","pi[2]","pi[3]"))
pairs(r1,pars=c("pi"))

###############################################
##for team placed second in the league stages##
###############################################

dat2 = list(N = length(data2$End_Position), m = 3, y = data2$End_Position)
r2 = sampling(m, dat2, c("pi"), iter=10000,
              control = list(max_treedepth=50, adapt_delta = 0.99))
r2 #p1 = 0.53; pi2 = 0.35; pi3=0.12
plot(r2, main="Team placed 2 after league stage")
traceplot(r2, pars = c("pi")) #properly mixing

#posterior plots
posterior2 <- as.array(r2)
mcmc_hist(posterior2, pars = c("pi[1]","pi[2]","pi[3]"))
mcmc_hist_by_chain(posterior2, pars = c("pi[1]","pi[2]","pi[3]"))
pairs(r2,pars=c("pi"))

###############################################
##for team placed third in the league stages##
###############################################

dat3 = list(N = length(data3$End_Position), m = 4, y = data3$End_Position)
r3 = sampling(m, dat3, c("pi"), iter=10000,
              control = list(max_treedepth=50, adapt_delta = 0.99))
r3 #p1 = 0.16; pi2 = 0.10; pi3=0.42; pi4=0.32
plot(r3, main="Team placed 3 after league stage")
traceplot(r3, pars = c("pi")) #properly mixing

#posterior plots
posterior3 <- as.array(r3)
mcmc_hist(posterior3, pars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
mcmc_hist_by_chain(posterior3, pars = c("pi[1]","pi[2]","pi[3]"))
pairs(r3,pars=c("pi"))

###############################################
##for team placed fourth in the league stages##
###############################################

dat4 = list(N = length(data4$End_Position), m = 4, y = data4$End_Position)
r4 = sampling(m, dat4, c("pi"), iter=10000,
              control = list(max_treedepth=50, adapt_delta = 0.99))
r4 #p1 = 0.16; pi2 = 0.10; pi3=0.42; pi4=0.32
plot(r4, main="Team placed 3 after league stage")
traceplot(r4, pars = c("pi")) #properly mixing

#posterior plots
posterior4 <- as.array(r4)
mcmc_hist(posterior4, pars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
mcmc_hist_by_chain(posterior4, pars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
pairs(r4,pars=c("pi"))

############################################
# Logistic regression through brms package #
############################################
library(brms)
options(mc.cores = parallel::detectCores())
#ipl.playoff is the cleaned dataset to be used
head(ipl.playoff)
str(ipl.playoff)
################
## get priors ##
################
pr <- get_prior(champion ~ 1 + seed + nrr + last5,
                data = ipl.playoff,
                family = bernoulli("logit"))

pr

################
prior_ma <- prior(normal(0, 10), class = "b") +
  prior(normal(0, 10), class = "Intercept")

prior_hs <- set_prior(horseshoe(df = 1, par_ratio = 0.5))

fit_ma1 <- brm( champion ~ 1 + seed + nrr + last5,
                data = ipl.playoff,
                family = bernoulli("logit"),
                prior = prior_ma,
                iter = 4000)

summary(fit_ma1)

# mean          se_     mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# b_Intercept  -3.62    0.04 2.11  -8.20  -4.97  -3.51  -2.17   0.26  3412    1
# b_seed2       1.80    0.02 1.24  -0.54   0.95   1.74   2.59   4.42  2507    1
# b_seed3      -1.50    0.03 1.60  -5.12  -2.43  -1.39  -0.43   1.37  2311    1
# b_seed4      -8.13    0.15 5.72 -21.45 -11.31  -7.02  -3.87  -0.07  1482    1
# b_nrr         1.37    0.03 1.52  -1.58   0.31   1.36   2.35   4.47  2865    1
# b_last5       0.57    0.01 0.61  -0.61   0.16   0.56   0.98   1.80  3135    1

plot(fit_ma1)
mcmc_plot(fit_ma1, type = "trace")

library(bayesplot)
mcmc_trace(fit_ma1, pars = c("b_Intercept", "b_nrr", "b_last5",
                             "b_seed2", "b_seed3", "b_seed4"))

#check convergence
modelposterior <- as.mcmc(fit_ma1) 
gelman.diag(modelposterior)

#plot conditional effect of seed
conditional_effects(fit_ma1, "seed")
condition = data.frame(seed = c(1,2,3,4))
conditional_effects(fit_ma1, "nrr", conditions = condition)
conditional_effects(fit_ma1, "last5", conditions = condition)

#prediction on 2021 data
newdata <- data.frame(ipl.2021[,-4])
pred <- predict(fit_ma1, newdata = newdata, type="response")[,1]
predictions <- data.frame(ipl.2021, pred)
colnames(predictions) <- c(colnames(ipl.2021), "P(champ)")
predictions$cp <- round((predictions$`P(champ)`)/sum((predictions$`P(champ)`)),4)
predictions

# > predictions
# seed    nrr last5 result P(champ)     cp
# 53    1  0.481     3      3  0.25950 0.2956
# 54    2  0.455     2      1  0.48875 0.5568
# 55    3 -0.140     4      4  0.11025 0.1256
# 56    4  0.587     3      2  0.01925 0.0219

###########################################################
#horseshoe prior

fit_ma2 <- brm( champion ~ 1 + seed + nrr + last5,
                data = ipl.playoff,
                family = bernoulli("logit"),
                prior = prior_hs)

summary(fit_ma2)
plot(fit_ma2, "seed")
plot(fit_ma2, c("last5", "nrr", "Intercept"))


