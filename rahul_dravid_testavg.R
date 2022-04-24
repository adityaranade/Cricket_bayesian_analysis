#packages required
library(rjags)
library(rstan)
library(plyr)
library(reshape2)
library(Rmisc)
library(ggplot2)
set.seed(23)

#test match data
d <- read.csv("rahul.test.data.csv",header=TRUE)

# keep only runs, location, year, opposition and mode of dismissal
d.test <- d[,c(2,4,6,7,9)]

#keep only the records which have entry in runs column
data <- d.test[complete.cases(d.test$runs),]
data$year <- as.factor(data$year)

#create ID varible for STAN model
data$season <- as.numeric(factor(data$year)) 

###############################
## Exploratory Data Analysis ##
###############################
#boxplot for test match runs scored yearly
ggplot(data=data, aes(x=year,y=runs, color=year))+
  geom_boxplot(show.legend = FALSE)+
  labs(x="Year",y="Runs")+
  ggtitle("Boxplot of test match runs scored every year")

#box plot for test matche runs scored by opponent
ggplot(data=data, aes(x=opposition,y=runs,color=opposition))+
  geom_boxplot()+
  labs(x="Opposition",y="Runs")+
  ggtitle("Boxplot of test match runs scored against opposition")

##############################
# Test match model
##############################
data$y = data$runs
#mean(data$y);sd(data$y)
hier_dat = list(n=nrow(data), y=data$y, nc=nlevels(as.factor(data$year)), 
                year=as.numeric(as.factor(data$year)))

#model with separate variance for all years
hier_model = "
model {
for (i in 1:n) {
y[i] ~ dnorm(mu[year[i]], 1/sigma2[year[i]])
}

for (c in 1:nc) {
mu[c] ~ dnorm(theta,1/tau2)
sigma2[c] = (sigma[c])^2
sigma[c] ~ dt(0,20,1)T(0,)

}
theta ~ dnorm(52,0.05)   
tau2 ~ dunif(10,60)
tau = sqrt(tau2)
}

"

# Run the model with separate variance for each year
hier_m = jags.model(textConnection(hier_model), n.chains=3, hier_dat, quiet=TRUE)
hier_r = coda.samples(hier_m, 
                       c("mu","sigma","theta","tau"), 
                       10000)
summary(hier_r)

# plot to compare posterior Test batting average by season- densities (separate variance for every year)
dr = ldply(hier_r, function(x) {as.data.frame(x)})
m = melt(dr[,grep("mu", names(dr))], variable.name="year", value.name="sample")
(ggplot(m, aes(x=sample, color=year))+geom_density()+geom_vline(xintercept=46.47,col=2)+
    labs(x="Batting Average in tests")+ggtitle("Test batting average by year"))

#plot the posterior standard deviation
m2 = melt(dr[,grep("sigma", names(dr))], variable.name="year", value.name="sample.avg")
(ggplot(m2, aes(x=sample.avg, color=year))+geom_density()+
    labs(x="Batting Average in tests")+ggtitle("Test batting average by year"))

# ci for model with separate variance for each year
# credible intervals to compare posterior Test batting averages by year/season
ci = ddply(m, .(year), summarize, lcl=quantile(sample, .025), ucl=quantile(sample, .975))
ci$Year = as.character(1996:2012)
(ggplot(ci, aes(x=lcl,xend=ucl,y=Year,yend=Year,col=Year))+geom_segment(lwd=2)+
geom_vline(xintercept=46.47,col=2)+ coord_cartesian(xlim=c(15,85)) +
labs(x="Test Batting Average by year",y="Year")+ggtitle("Test batting average by year"))

# trace plot for year 1998, 2003 and 2010
par(mfrow=c(1,1))
plot(hier_r[,c("theta","tau")], trace=TRUE)
traceplot(hier_r[,c("mu[15]")], main="Year 2010")
traceplot(hier_r[,c("mu[8]")],main="Year 2003")
traceplot(hier_r[,c("mu[3]")],main="Year 1998")

#posterior predictive p-value to check if the model is decent
ttest = t.test(runs~location,data)
ttest
obs_test_statistic_test = diff(ttest$estimate)

mu = hier_r[,grep("mu",varnames(hier_r))][[1]]
sigma = (hier_r[,grep("sigma",varnames(hier_r))][[1]])
tmp = data
tmp$location = as.factor(tmp$location)
n = nrow(data)

yrep = rdply(nrow(mu), {
  tmp$runs = rnorm(n, mu[hier_dat$year], sigma[hier_dat$year])
  data.frame(test_statistic_test=diff(t.test(runs~location, tmp)$estimate))
})

sum(yrep$test_statistic_test>obs_test_statistic_test)/max(yrep$.n) 
#posterior predictive p value is 0.4321 which is not too less or too low

#posterior predictive plot test career
qplot(test_statistic_test, data=yrep, geom="histogram", bins=100) +
  geom_vline(xintercept=obs_test_statistic_test,col=2)+
  ggtitle("Posterior predictive plot")

# Posterior predictive plot indicates the t test value from replicated data 
# is very close to the test statistic from the observed statistic which 
# indicates the model is working well.

################
## STAN Code  ##
################
library(rstan)
hier_model2 = "
data {
int<lower=1> n; //number of years
int<lower=1> N; //number of observations
real y[N];
int season[N]; //season corresponding to observation
}

parameters {
real mu[n];
real<lower=0> sigma[n];
real theta;
real<lower=0> tau2;
}

transformed parameters {
real<lower=0> tau;
tau = sqrt(tau2);
}

model {
//Priors
target += normal_lpdf(theta |50,16);
target += cauchy_lpdf(tau2 |10,60);

//Likelihood
for(k in 1:n) {
target += normal_lpdf(mu[k] |theta,tau);
target += cauchy_lpdf(sigma[k] |10,50);
}
for(j in 1:N) {
target += normal_lpdf(y[j] |mu[season[j]],sigma[season[j]]);
}
}
"
stan_m = stan_model(model_code=hier_model2)
dat = list(y = data$runs, n = max(data$season), season = data$season, N = nrow(data))
r = sampling(stan_m, dat, c("mu", "sigma", "theta"), iter=10000, 
              control = list(max_treedepth=50, adapt_delta = 0.99))
r
plot(r, pars = "mu")
plot(r, pars = "sigma")
###############################################################################