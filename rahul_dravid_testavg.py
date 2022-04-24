#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pystan
import arviz # to plot posterior and traceplots


# In[12]:


data = pd.read_csv("D:/Study Iowa State University/Projects/Rahul Dravid dataset/rahul.test.data.csv")
print(data.head())


# In[24]:


data['season'] = data['season'].astype('int')
print(data.head())


# In[15]:


df = data[['runs', 'year', 'season' ]].copy() # select only required columns
df.dropna(axis = 0, how = 'any', inplace = True) #delete rows with NA
print(df.head())


# In[17]:


#STAN model
hier_model = """
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
"""


# In[29]:


dat_hm = {
    "N": len(df.index),
    "n": 17,
    "y": df["runs"],
    "season": df["season"]
}
fit = pystan.stan(model_code = hier_model, data=dat_hm, iter=10000, chains=1)


# In[30]:


print(fit)  # parameter estimates
arviz.plot_trace(fit) # posterior density and traceplots


# In[ ]:




