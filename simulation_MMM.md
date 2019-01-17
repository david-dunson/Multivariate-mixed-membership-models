

# Sample from MMM model


This tutorial shows how to sample from the posterior of the proposed MMM model, following the steps highlighted in Algorithm 1 of the paper [**Multivariate mixed membership modeling: Inferring domain-specific risk profiles**](https://arxiv.org/pdf/1901.05191.pdf).

* First we generated data from different simulation scenarios following the tutorial in [`simulated_data.md`](simulated_data.md), and having the correspondent R data in the folder `simulated_data/`.


```r
library(MMM)
## loaded as required 
## library(mvnfast)
## library(MCMCpack)
```

##  Number of profiles correctly specified 

We consider the 4 simulation scenarios proposed in Section 6 of the paper.

### Scenario 1


```r
## Load data from 
load('simulated_data/Sim1_data_and_true_parameters.Rdata')


## seed for reproducibility
set.seed(1040)

## set prior and startig values
param        = list(
                  omega  = matrix(NA,data$n,2),
                  psi    = matrix(0,data$n,2),
                  ki     = matrix(NA,data$n,2),
                  mu     = c(0,0),
                  lambda = matrix(0.5,data$n,2),
                  kern   = array(NA,dim  = c(4,data$p,2)),
                  zi     = true_parameter$sbj_profiles,
                  Sigma  =  diag(1,2)
)





prior = list(
           mu     = c(0,0),
           nu     = 2,
           Psi    = diag(1,2),
           Sigma  = diag(1,2)
)




Scenario1= estimate_MMM(data=data ,## list containing: 
                                    # n = number of subjects; 
                                    # p = number of variables;
                                    # gj = group indicators;
                                    # data nxp  data-matrix
                 prior = prior,    ## (optional) default as in paper
                 param = param,  ## (optional) list of parameter initialization
                 nrep     = 5000) ## number of MCMC samples 

save(Scenario1,file = 'results/Sim1_MMM_results.Rdata')
rm(list = ls());gc();
```


### Scenario 2


```r
## Load data from 
load('simulated_data/Sim2_data_and_true_parameters.Rdata')

## seed for reproducibility
set.seed(10300)

## set prior and startig values
param        = list(
                  omega  = matrix(NA,data$n,2),
                  psi    = matrix(0,data$n,2),
                  ki     = matrix(NA,data$n,2),
                  mu     = c(0,0),
                  lambda = matrix(0.5,data$n,2),
                  kern   = array(NA,dim  = c(4,data$p,2)),
                  zi     = true_parameter$sbj_profiles,
                  Sigma  =  diag(1,2)
)





prior = list(
           mu     = c(0,0),
           nu     = 2,
           Psi    = diag(1,2),
           Sigma  = diag(1,2)
)




Scenario2= estimate_MMM(data=data ,## list containing: 
                                    # n = number of subjects; 
                                    # p = number of variables;
                                    # gj = group indicators;
                                    # data nxp  data-matrix
                 prior = prior,    ## (optional) default as in paper
                 param = param,  ## (optional) list of parameter initialization
                 nrep     = 5000) ## number of MCMC samples 

save(Scenario2,file = 'results/Sim2_MMM_results.Rdata')
rm(list = ls());gc();
```

### Scenario 3


```r
## Load data from 
load('simulated_data/Sim3_data_and_true_parameters.Rdata')

## seed for reproducibility
set.seed(1320)

## set prior and startig values
param        = list(
                  omega  = matrix(NA,data$n,2),
                  psi    = matrix(0,data$n,2),
                  ki     = matrix(NA,data$n,2),
                  mu     = c(0,0),
                  lambda = matrix(0.5,data$n,2),
                  kern   = array(NA,dim  = c(4,data$p,2)),
                  zi     = true_parameter$sbj_profiles,
                  Sigma  =  diag(1,2)
)





prior = list(
           mu     = c(0,0),
           nu     = 2,
           Psi    = diag(1,2),
           Sigma  = diag(1,2)
)




Scenario3= estimate_MMM(data=data ,## list containing: 
                                    # n = number of subjects; 
                                    # p = number of variables;
                                    # gj = group indicators;
                                    # data nxp  data-matrix
                 prior = prior,    ## (optional) default as in paper
                 param = param,    ## (optional) list of parameter initialization
                 nrep     = 5000) ## number of MCMC samples 

save(Scenario3,file = 'results/Sim3_MMM_results.Rdata')
rm(list = ls());gc();
```

### Scenario 4


```r
## Load data from 
load('simulated_data/Sim4_data_and_true_parameters.Rdata')

## seed for reproducibility
set.seed(1984)
## set prior and startig values
param        = list(
                  omega  = matrix(NA,data$n,2),
                  psi    = matrix(0,data$n,2),
                  ki     = matrix(NA,data$n,2),
                  mu     = c(0,0),
                  lambda = matrix(0.5,data$n,2),
                  kern   = array(NA,dim  = c(4,data$p,2)),
                  zi     = true_parameter$sbj_profiles,
                  Sigma  = diag(1,2)
)





prior = list(
           mu     = c(0,0),
           nu     = 2,
           Psi    = diag(1,2),
           Sigma  = diag(1,2)
)





Scenario4= estimate_MMM(data=data ,## list containing: 
                                    # n = number of subjects; 
                                    # p = number of variables;
                                    # gj = group indicators;
                                    # data nxp  data-matrix
                 prior = prior,    ## (optional) default as in paper
                 param = param,  ## (optional) list of parameter initialization
                 nrep     = 5000) ## number of MCMC samples 

save(Scenario4,file = 'results/Sim4_MMM_results.Rdata')
rm(list = ls());gc();
```

## Misspecification: more than two pure types


```r
load('simulated_data/data_and_true_parameters_misspecified.Rdata')

## seed for reproducibility
set.seed(1984)
## set prior and startig values
param        = list(
                  omega  = matrix(NA,data$n,2),
                  psi    = matrix(0,data$n,2),
                  ki     = matrix(NA,data$n,2),
                  mu     = c(0,0),
                  lambda = matrix(0.5,data$n,2),
                  kern   = array(NA,dim  = c(4,data$p,2)),
                  zi     = matrix(sample(1:2,size = data$n * data$p,replace = T),data$n,data$p)

                  Sigma  = diag(1,2)
)


prior = list(
           mu     = c(0,0),
           nu     = 2,
           Psi    = diag(1,2),
           Sigma  = diag(1,2)
)


Scenario_miss = estimate_MMM(data=data ,## list containing: 
                                    # n = number of subjects; 
                                    # p = number of variables;
                                    # gj = group indicators;
                                    # data nxp  data-matrix
                 prior = prior,    ## (optional) default as in paper
                 param = param,  ## (optional) list of parameter initialization
                 nrep     = 5000) ## number of MCMC samples 

save(Scenario_miss,file = 'results/Sim_miss_MMM_results.Rdata')
rm(list = ls());gc();
```

Once results have been saved we can load and analyze them following the conde in [`plots_and_tables.md`](plots_and_tables.md)
