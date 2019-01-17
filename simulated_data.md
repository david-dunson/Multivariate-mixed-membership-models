

# Simulate Data 


This tutorial shows how to reproduce the simulated data analyzed in Section 6 of the paper [**Multivariate mixed membership modeling: Inferring domain-specific risk profiles**](https://arxiv.org/pdf/1901.05191.pdf).


```r
## Libraries
library(mvnfast)
library(MCMCpack)
```


# Number of profiles correctly specified 

We consider the 4 simulation scenarios proposed in Section 6 of the paper.

## Kernels
Kernels are shared across the 4 simulation scenarios and are defined as


```r
## seed for reproducibility
set.seed(12425)

## data specific quantities
data          = list()
data$n        = 1000
data$p1       = 5 
data$p2       = 5 
data$p        = data$p1+data$p2
data$gj       = c(rep(1,data$p1),rep(2,data$p2)) # group indicator


true_parameter = list()
## G = 1; d=4; p=10; H=2 
true_parameter$kern = array(NA,dim = c(4,data$p,2))


## kernels for g=1
for(j in 1:data$p1)
{
 true_parameter$kern[,j,1] = rdirichlet(1,c(10,3,2,1))
 true_parameter$kern[,j,2] = rdirichlet(1,c(1,1,1,11))
}


## kernels for g=2
for(j in (data$p1+1):data$p)
{
 true_parameter$kern[,j,1] = c(rdirichlet(1,c(5,5,1,0)))
 true_parameter$kern[,j,2] = c(rdirichlet(1,c(1,1,1,8)))
}
```

## Generate data 

We need to generate the true score vectors for each subject, and generate the data according to model (4) in Section 3 of the paper.

### Scenario 1

In the first Scenario profiles distribution P is a normal truncated in the unit square with parameters
 **μ**=(1/2,1/2)<sup>T </sup> and **Ʃ** = (Ʃ<sub>11</sub>, Ʃ<sub>21</sub>, Ʃ<sub>12</sub>, Ʃ<sub>22</sub>)<sup>T</sup> = (0.05,0.02,0.02,0.05)<sup>T</sup>.  This formulation induces positive dependence between the two scores with their distribution having ellipsoid contours truncated at the borders.
 

We can generate profiles and data using the following code


```r
## Truncated normal distribution
## hit or miss simulation
Sig = matrix(c(0.05,0.02,0.02,0.05),2,2)
mu  = c(0.5,0.5)

true_parameter$sbj_prob      = matrix(NA,data$n,2)
i = 0
while(i <= data$n)
{
 tmp = rmvn(1,mu,Sig)
 if((tmp[1] < 1) & (tmp[1]>0) & (tmp[2] < 1) & (tmp[2]>0)){
 true_parameter$sbj_prob[i,] = tmp
 i = i+1
 }
}
rm(tmp);gc(FALSE);


## Assign each subject to a profile with probability given by the scores
true_parameter$sbj_profiles  = matrix(NA,data$n,data$p)

for(i in 1:data$n)
{
 for(j in 1:data$p){
  true_parameter$sbj_profiles[i,j] = sample(1:2,
                                            size = 1,
                                            prob = c(1-true_parameter$sbj_prob[i,data$gj[j]],true_parameter$sbj_prob[i,data$gj[j]]))
 }
}


## Generate the data
data$data = matrix(NA,data$n,data$p)

for(i in 1:data$n)
{
 for(j in 1:data$p)
 {
 data$data[i,j]  = sample(1:4,
                          size = 1,
                          prob = true_parameter$kern[,j,true_parameter$sbj_profiles[i,j]] )
 }
}

## Save data & parameters
save(data,true_parameter, file = 'simulated_data/Sim1_data_and_true_parameters.Rdata')
```

### Scenario 2
 In this scenario we consider the multivariate logistic normal distribution proposed in section (4) of the paper, having **μ** = (-1.2,1.0)<sup>T</sup> and vec(**Ʃ**) = (3.0,-2.4,-2.4,3.5)<sup>T</sup>. 




```r
true_parameter$mu            = c(-1.2,1)
true_parameter$Sigma         = rbind(c(3.0,-2.4),c(-2.4,3.5))

true_parameter$sbj_prob      = apply(rmvn(data$n, 
                                          mu = true_parameter$mu,
                                          sigma = true_parameter$Sigma),2,plogis)


## Assign each subject to a profile with probability given by the scores
true_parameter$sbj_profiles  = matrix(NA,data$n,data$p)

for(i in 1:data$n)
{
 for(j in 1:data$p)
 {
 true_parameter$sbj_profiles[i,j] =
    sample(1:2,
           size = 1,
           prob  = c(1-true_parameter$sbj_prob[i,data$gj[j]],true_parameter$sbj_prob[i,data$gj[j]]))
 }
}


## Generate the data
data$data = matrix(NA,data$n,data$p)
for(i in 1:data$n)
{
 for(j in 1:data$p)
 {
  data$data[i,j]  = sample(1:4,
                           size = 1,
                           prob = true_parameter$kern[,j,true_parameter$sbj_profiles[i,j]])
 }
}

## Save data & parameters
save(data,true_parameter, file = 'simulated_data/Sim2_data_and_true_parameters.Rdata')
```

### Scenario 3

In the third scenario, we rely on the generate the profiles from a single MM model, having profile distribution shared by all variables; we generate this profile from a uniform distribution.




```r
true_parameter$sbj_prob      = matrix(NA,data$n,2) 
true_parameter$sbj_prob[,1]  = true_parameter$sbj_prob[,2] = runif(data$n)


## Assign each subject to a profile with probability given by the scores
true_parameter$sbj_profiles  = matrix(NA,data$n,data$p)

for(i in 1:data$n)
{
 for(j in 1:data$p)
 {
  true_parameter$sbj_profiles[i,j] = sample(1:2,
                                            size = 1,
                                            prob = c(1-true_parameter$sbj_prob[i,data$gj[j]], true_parameter$sbj_prob[i,data$gj[j]]))
 }
}

## Generate the data
data$data = matrix(NA,data$n,data$p)

for(i in 1:data$n)
{
 for(j in 1:data$p)
 {
  data$data[i,j]  = sample(1:4,
                           size = 1,
                           prob = true_parameter$kern[,j,true_parameter$sbj_profiles[i,j]] )
 }
}

## Save data & parameters
save(data,true_parameter, file = 'simulated_data/Sim3_data_and_true_parameters.Rdata')
```


### Scenario 4

Finally, in the fourth simulation scenario, we consider P to be the product of two independent uniforms, forcing independence in the variables belonging to different groups, which translates into the case in which two separate models for the groups represents the correctly specified model. 


```r
true_parameter$sbj_prob     = matrix(NA,data$n,2)
true_parameter$sbj_prob[,1] = runif(data$n)
true_parameter$sbj_prob[,2] = runif(data$n)


## Assign each subject to a profile with probability given by the scores
true_parameter$sbj_profiles  = matrix(NA,data$n,data$p)

for(i in 1:data$n)
{
 for(j in 1:data$p)
 {
  true_parameter$sbj_profiles[i,j] = sample(1:2,
                                            size = 1,
                                            prob = c(1-true_parameter$sbj_prob[i,data$gj[j]],true_parameter$sbj_prob[i,data$gj[j]]))
 }
}


## Generate the data
data$data = matrix(NA,data$n,data$p)

for(i in 1:data$n)
{
 for(j in 1:data$p)
 {
  data$data[i,j]  = sample(1:4,
                           size = 1,
                           prob = true_parameter$kern[,j,true_parameter$sbj_profiles[i,j]] )
 }
}


## Save data & parameters
save(data,true_parameter, file = 'simulated_data/Sim4_data_and_true_parameters.Rdata')
```
## Misspecification: more than two pure types

## Kernels


```r
## seed for reproducibility
set.seed(12425)

data          = list()
data$n        = 1000
data$p1       = 5 
data$p2       = 5 
data$p        = data$p1+data$p2
data$gj       = c(rep(1,data$p1),rep(2,data$p2)) # group indicator
data$H        = 4

true_parameter = list()
## G=1; d=4; p=10; H=2 
true_parameter$kern = array(NA,dim = c(4,data$p,data$H))


## kernels for g=1
for(j in 1:data$p1)
{
  true_parameter$kern[,j,1] = c(0.85,0.05,0.05,0.05)
  true_parameter$kern[,j,2] = c(0.05,0.85,0.05,0.05)
  true_parameter$kern[,j,3] = c(0.05,0.05,0.85,0.05)
  true_parameter$kern[,j,4] = c(0.05,0.05,0.05,0.85)
}

## kernels for g=2
for(j in (data$p1+1):data$p){
 true_parameter$kern[,j,1] = c(rdirichlet(1,c(5,5,1,0)))
 true_parameter$kern[,j,2] = c(rdirichlet(1,c(1,1,1,8)))
}
```

## Generate Data


```r
true_parameter$sbj_prob            = list()
true_parameter$sbj_prob[[1]]       = rdirichlet(data$n,c(0.25,0.25,0.25,0.25))
true_parameter$sbj_prob[[2]]       = runif(data$n)



## Assign each subject to a profile with probability given by the scores
true_parameter$sbj_profiles  = matrix(NA,data$n,data$p)

for(i in 1:data$n)
{
 for(j in 1:data$p1)
 {
  true_parameter$sbj_profiles[i,j] = sample(1:data$H,
                                            size = 1,
                                            prob = true_parameter$sbj_prob[[1]][i,])
 }
 for(j in (data$p1 + 1):data$p)
 {
  true_parameter$sbj_profiles[i,j] = sample(1:2,
                                              size = 1,
                                              prob = c(1-true_parameter$sbj_prob[[2]][i],true_parameter$sbj_prob[[2]][i]))
 }
}



## Generate the data
data$data = matrix(NA,data$n,data$p)
for(i in 1:data$n){
 for(j in 1:data$p)
 {
 data$data[i,j]  = sample(1:4,size = 1,prob = true_parameter$kern[,j,true_parameter$sbj_profiles[i,j]] )
 }
}

## Save data & parameters
save(data,true_parameter, file = 'simulated_data/data_and_true_parameters_misspecified.Rdata')
rm(list =ls());gc();
```
