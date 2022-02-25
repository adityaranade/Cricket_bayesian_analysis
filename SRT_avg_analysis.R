#About the project

##Bayesian data analysis to compare yearly 
##batting average of Sachin Tendulkar


#load data and select the required columns
d <- read.csv(file.choose(),header=TRUE)
data <- d[,c(1,10,20,21,22,23)]
data$year <- as.character(data$year)

data.test <- data[1:338,]; as.numeric(as.factor(data.test$year))
data.odi <- data[339:790,]; as.numeric(as.factor(data.odi$year))

######################
## Test match model ##
######################
library(rjags)
library(ggplot2)
library(plyr)
library(reshape2)
library(Rmisc)
set.seed(23)
data.test$y = data.test$Runs
mean(data.test$y);sd(data.test$y)
hier_dat = list(n=nrow(data.test), y=data.test$y, 
                nc=nlevels(as.factor(data.test$year)), 
                year=as.numeric(as.factor(data.test$year)))

hier_model1 = "
model {
for (i in 1:n) {
y[i] ~ dnorm(mu[year[i]], 1/sigma2)
}

for (c in 1:nc) {
mu[c] ~ dnorm(theta,1/tau2)
}

sigma2 = sigma*sigma
sigma ~ dt(0,20,1)T(0,)
theta ~ dnorm(50,0.0625)   
tau2 = tau*tau
tau ~ dunif(10,50)
}

"

hier_m1 = jags.model(textConnection(hier_model1), n.chains=3, 
                     hier_dat, quiet=TRUE)
hier_r1 = coda.samples(hier_m1, 
                       c("mu","sigma","theta","tau"), 
                       10000)
summary(hier_r1)


# plot to compare posterior Test batting average by season- densities
dr = ldply(hier_r1, function(x) {as.data.frame(x)})
m = melt(dr[,grep("mu", names(dr))], variable.name="year", value.name="sample")
#par(mfrow=c(1,2))
(ggplot(m, aes(x=sample, color=year))+geom_density()+
    geom_vline(xintercept=47.92,col=2)+
    labs(x="Batting Average in tests")+ggtitle("Test batting average by year"))

# credible intervals to compare posterior Test batting averages by year/season
ci = ddply(m, .(year), summarize, lcl=quantile(sample, .025), 
           ucl=quantile(sample, .975))
ci$year = as.character(1989:2013)
(ggplot(ci, aes(x=lcl,xend=ucl,y=year,yend=year,col=year))+geom_segment(lwd=2)+
geom_vline(xintercept=47.92,col=2)+ coord_cartesian(xlim=c(15,80)) +
labs(x="Test Batting Average by year")+ggtitle("Test batting average by year"))

par(mfrow=c(1,1))
plot(hier_r1[,c("theta","sigma","tau")], trace=TRUE)
traceplot(hier_r1[,c("mu[22]")], main="year 2010")
traceplot(hier_r1[,c("mu[12]")],main="year 2000")
traceplot(hier_r1[,c("mu[3]")],main="year 1991")

ttest = t.test(Runs~location,data.test)
ttest
obs_test_statistic_test = diff(ttest$estimate)

mu = hier_r1[,grep("mu",varnames(hier_r1))][[1]]
sigma = (hier_r1[,grep("sigma",varnames(hier_r1))][[1]])
tmp = data.test
tmp$location = as.factor(tmp$location)
n = nrow(data.test)

yrep = rdply(nrow(mu), {
  tmp$Runs = rnorm(n, mu[hier_dat$year], sigma)
  data.frame(test_statistic_test=diff(t.test(Runs~location, tmp)$estimate))
})

qplot(test_statistic_test, data=yrep, geom="histogram", bins=100) +
  geom_vline(xintercept=obs_test_statistic_test,col=2)

#####################
## ODI match model ##
#####################
set.seed(23)
data.odi$y = data.odi$Runs
par(mfrow=c(1,2))

mean(data.odi$y);var(data.odi$y)
hier_dat2 = list(n=nrow(data.odi), y=data.odi$y, 
                 nc=nlevels(as.factor(data.odi$year)), 
                year=as.numeric(as.factor(data.odi$year)))

hier_model2 = "
model {
for (i in 1:n) {
y[i] ~ dnorm(mu[year[i]], 1/sigma2)
}

for (c in 1:nc) {
mu[c] ~ dnorm(theta,1/tau2)
}

sigma2 = sigma*sigma
sigma ~ dt(0,40,1)T(0,)
theta ~ dnorm(50,0.0625)   
tau2 = tau*tau
tau ~ dunif(10,50)
}
"
hier_m2 = jags.model(textConnection(hier_model2), 
                     n.chains=3, hier_dat2, quiet=TRUE)
hier_r2 = coda.samples(hier_m2, 
                       c("mu","sigma","theta","tau"), 
                       10000)
summary(hier_r2)

# plot to compare posterior ODI batting average by season- densities
dr2 = ldply(hier_r2, function(x) {as.data.frame(x)})
m2 = melt(dr2[,grep("mu", names(dr2))], 
          variable.name="year", value.name="sample")
(ggplot(m2, aes(x=sample, color=year))+geom_vline(xintercept=40.77,col=2)+
    labs(x="Batting Average in ODIs")+geom_density()+
    ggtitle("ODI Batting average by year"))

# credible intervals to compare posterior ODI batting averages by year/season
ci2 = ddply(m2, .(year), summarize, lcl=quantile(sample, .025), 
            ucl=quantile(sample, .975))
ci2$year = as.character(1989:2012)
(ggplot(ci2, aes(x=lcl,xend=ucl,y=year,yend=year,col=year))+geom_segment(lwd=2)
+labs(x="ODI Batting Average by year")+ geom_vline(xintercept=40.77,col=2)+
 ggtitle("ODI Batting average by year"))

par(mfrow=c(1,1))
plot(hier_r2[,c("theta","sigma","tau")], trace=TRUE)
traceplot(hier_r2[,c("mu[10]")], main="year 1998")
traceplot(hier_r2[,c("mu[17]")],main="year 2005")
traceplot(hier_r2[,c("mu[24]")],main="year 2012")


ttest2 = t.test(Runs~location,data.odi)
ttest2
obs_test_statistic_odi = diff(ttest2$estimate)

mu2 = hier_r2[,grep("mu",varnames(hier_r2))][[1]]
sigma2 = (hier_r2[,grep("sigma",varnames(hier_r2))][[1]])
tmp2 = data.odi
tmp2$location = as.factor(tmp2$location)
n2 = nrow(data.odi)

yrep2 = rdply(nrow(mu2), {
  tmp2$Runs = rnorm(n2, mu2[hier_dat2$year], sigma2)
  data.frame(test_statistic_odi=diff(t.test(Runs~location, tmp2)$estimate))
})

qplot(test_statistic_odi, data=yrep2, geom="histogram", bins=100) +
  geom_vline(xintercept=obs_test_statistic_odi,col=2)
