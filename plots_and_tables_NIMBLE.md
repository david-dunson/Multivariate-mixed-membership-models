Plots & Tables: Supplementary Material
======================================

We first some libraries

``` r
library(tidyverse)
library(reshape2)
library(ggplot2)
library(mvnfast)
library(gridExtra)
library(nimble)
library(MCMCpack)
library(knitr)
```

We also need the package <tt>mvtnorm</tt> to be installed, as we use the
function <tt>mvtnorm::pmvnorm</tt>, and to load a custom function to
embed <tt>.Rdata</tt> objects in a named list. The function is taken
from
[stackoverflow](https://stackoverflow.com/questions/14757668/combine-multiple-rdata-files-containing-objects-with-the-same-name-into-one-sin)

``` r
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}
```

Load results & simulated data
-----------------------------

``` r
Scenario1 = c(loadRData('results/Sim1_MMM_results.Rdata'),
              loadRData('results/Sim1_MM_results.Rdata'))

Scenario1$Nimble1 = readRDS("results/Sim1_nimble_model_1.rds")
Scenario1$Nimble2 = readRDS("results/Sim1_nimble_model_2.rds")


load('simulated_data/Sim1_data_and_true_parameters.Rdata')
Scenario1$data  = data
Scenario1$true_parameter  = true_parameter
rm(data,true_parameter)

Scenario2 = c(loadRData('results/Sim2_MMM_results.Rdata'),
              loadRData('results/Sim2_MM_results.Rdata'))

Scenario2$Nimble1 = readRDS("results/Sim2_nimble_model_1.rds")
Scenario2$Nimble2 = readRDS("results/Sim2_nimble_model_2.rds")

load('simulated_data/Sim2_data_and_true_parameters.Rdata')
Scenario2$data  = data
Scenario2$true_parameter  = true_parameter
rm(data,true_parameter)

Scenario3 = c(loadRData('results/Sim3_MMM_results.Rdata'),
              loadRData('results/Sim3_MM_results.Rdata'))


Scenario3$Nimble1 = readRDS("results/Sim3_nimble_model_1.rds")
Scenario3$Nimble2 = readRDS("results/Sim3_nimble_model_2.rds")
        
load('simulated_data/Sim3_data_and_true_parameters.Rdata')
Scenario3$data  = data
Scenario3$true_parameter  = true_parameter
rm(data,true_parameter)
              


Scenario4 = c(loadRData('results/Sim4_MMM_results.Rdata'),
              loadRData('results/Sim4_MM_results.Rdata'))


Scenario4$Nimble1 = readRDS("results/Sim4_nimble_model_1.rds")
Scenario4$Nimble2 = readRDS("results/Sim4_nimble_model_2.rds")

load('simulated_data/Sim4_data_and_true_parameters.Rdata')
Scenario4$data  = data
Scenario4$true_parameter  = true_parameter
rm(data,true_parameter)
gc()
```

### Plot Scenario 1

``` r
num_of_samp = 1000 ## number of samples

membership_scores = matrix(0,num_of_samp,2)
for(r in (5000 - num_of_samp + 1):5000)
{
 membership_scores[r -(5000 - num_of_samp) , ] = plogis(rmvn(1, mu = Scenario1$mu[r,],sigma = Scenario1$Sigma[r,,]))
}


## bivariate truncated normal density function
dTNORM = function(x,mu,Sig) 
{
    A = c(t(x-mu) %*% solve(Sig) %*% (x-mu))
   return( exp( -0.5*A)/mvtnorm::pmvnorm( lower = c(0,0), upper = c(1,1), mean = mu ,sigma = Sig)[1])
}


xx = seq(0.01,0.99, l = 30)
yy = seq(0.01,0.99, l = 30)
eval_points = expand.grid(xx,yy)
eval_points = cbind(eval_points, value = 
apply(eval_points,1, function(x) dTNORM(x,mu = c(0.5,0.5),
                                       Sig = matrix(c(0.05,0.02,0.02,0.05),2,2))))



## vs nimble MCMC estimates
tmp3 = cbind(
             apply(Scenario1$Nimble1[3001:4000,str_detect(colnames(Scenario1$Nimble1),'tilde_gamma')], 1, function(x) rdirichlet(1,x)[,1]),
             apply(Scenario1$Nimble2[3001:4000,str_detect(colnames(Scenario1$Nimble2),'tilde_gamma')], 1, function(x) rdirichlet(1,x)[,2]))
       

ggdata3 = data.frame(m1= membership_scores[,1],
                     m2 = membership_scores[,2],
                     l1 = tmp3[,1],
                     l2 = tmp3[,2]) 


plot_scenario1_NIMBLE = 
                       ggplot(ggdata3) +
                       geom_contour(aes(Var1,Var2,z = value),bins = 13,data = eval_points,col='grey50',size = 1.0)+
                       geom_point(aes(m1,m2), alpha = 0.6,size = 0.8) +
                       geom_point(aes(l1,l2),col ="#00249C",shape = 4, alpha = 0.6,size = 0.8) +
                       scale_x_continuous(breaks = seq(0,1,0.25),labels = sprintf("%.2f" ,seq(0,1,0.25)),expand=c(0.009,0.02))+
                       scale_y_continuous(breaks = seq(0.25,1,0.25),labels = sprintf("%.2f" ,seq(0.25,1,0.25)),expand=c(0.001,0.001))+
                       facet_wrap(~I('SCENARIO 1'))+
                       labs(x  = expression(hat(lambda)^{(1)}),y =  expression(hat(lambda)^{(2)})) + 
                       theme_bw()
```

### Plot Scenario 2

``` r
membership_scores = matrix(0,num_of_samp,2)

for(r in (5000 - num_of_samp + 1):5000)
{
 membership_scores[r -(5000 - num_of_samp + 1) , ] = plogis(rmvn(1, mu = Scenario2$mu[r,],sigma = Scenario2$Sigma[r,,]))
}



## density function for bivariate logistic normal
dMLND = function(x,mu,Sigma)
{
  logitx = log(x/(1-x))
  
out = 
  1/(sqrt(2*pi*det(Sigma))*prod(x)*prod(1-x))*
  exp( -0.5*( t(logitx  - mu) %*%solve(Sigma) %*% (logitx  - mu)))
return(out)
}

xx = seq(0.01,0.99, l = 100)
yy = seq(0.01,0.99, l = 100)
eval_points = expand.grid(xx,yy)
eval_points = cbind(eval_points, value = 
apply(eval_points,1, function(x) dMLND(x,mu = Scenario2$true_parameter$mu,
                                       Sigma = Scenario2$true_parameter$Sigma)))



## vs nimble MCMC estimates
tmp3 = cbind(
             apply(Scenario2$Nimble1[3001:4000,str_detect(colnames(Scenario2$Nimble1),'tilde_gamma')], 1, function(x) rdirichlet(1,x)[,2]),
             apply(Scenario2$Nimble2[3001:4000,str_detect(colnames(Scenario2$Nimble2),'tilde_gamma')], 1, function(x) rdirichlet(1,x)[,2]))
       

ggdata3 = data.frame(m1= membership_scores[,1],
                     m2 = membership_scores[,2],
                     l1 = tmp3[,1],
                     l2 = tmp3[,2]) 


plot_scenario2_NIMBLE = 
                       ggplot(ggdata3) +
                       geom_contour(aes(Var1,Var2,z = value),bins = 150,data = eval_points,col='grey50',size = 1.0)+
                       geom_point(aes(m1,m2), alpha = 0.6,size = 0.8) +
                       geom_point(aes(l1,l2),col = '#00249C',shape = 4,  alpha = 0.6,size = 0.8) +
                       scale_x_continuous(breaks = seq(0,1,0.25),labels = sprintf("%.2f" ,seq(0,1,0.25)),expand=c(0.009,0.02))+
                       scale_y_continuous(breaks = seq(0.25,1,0.25),labels = sprintf("%.2f" ,seq(0.25,1,0.25)),expand=c(0.002,0.002))+
                       facet_wrap(~I('SCENARIO 2'))+
                       labs(x  = expression(hat(lambda)^{(1)}),y =  expression(hat(lambda)^{(2)})) + 
                       theme_bw()
```

### Plot Scenario 3

``` r
membership_scores = matrix(0,num_of_samp,2)
for(r in (5000 - num_of_samp + 1):5000)
{
 membership_scores[r - (5000 - num_of_samp), ] = plogis(rmvn(1, mu = Scenario3$mu[r,],sigma = Scenario3$Sigma[r,,]))

}




## vs nimble MCMC estimates
tmp3 = cbind(
             apply(Scenario3$Nimble1[3001:4000,str_detect(colnames(Scenario3$Nimble1),'tilde_gamma')], 1, function(x) rdirichlet(1,x)[,1]),
             apply(Scenario3$Nimble2[3001:4000,str_detect(colnames(Scenario3$Nimble2),'tilde_gamma')], 1, function(x) rdirichlet(1,x)[,1]))


ggdata3 = data.frame(m1= membership_scores[,1],
                     m2 = membership_scores[,2],
                     l1 = tmp3[,1],
                     l2 = tmp3[,2]) 


plot_scenario3_NIMBLE = 
                        ggplot(ggdata3) +
                        geom_abline(intercept = 0, slope = 1,col = "grey50",size=1.0 ) +
                        geom_point(aes(m1,m2), alpha = 0.6,size = 0.8) +
                        geom_point(aes(l1,l2),col = '#00249C',shape = 4, alpha = 0.6,size = 0.8) +
                        scale_x_continuous(breaks = seq(0,1,0.25),labels = sprintf("%.2f" ,seq(0,1,0.25)),expand=c(0.009,0.02))+
                        scale_y_continuous(breaks = seq(0.25,1,0.25),labels = sprintf("%.2f" ,seq(0.25,1,0.25)),expand=c(0.001,0.001))+
                        facet_wrap(~I('SCENARIO 3'))+
                        labs(x  = expression(hat(lambda)^{(1)}),y =  expression(hat(lambda)^{(2)})) + 
                        theme_bw()
```

### Plot Scenario 4

``` r
membership_scores = matrix(0,num_of_samp,2)
for(r in  (5000 - num_of_samp + 1):5000)
{
 membership_scores[r -(5000 - num_of_samp) , ] = plogis(rmvn(1, mu = Scenario4$mu[r,],sigma = Scenario4$Sigma[r,,]))
}


## vs nimble MCMC estimates
tmp3 = cbind(
             apply(Scenario4$Nimble1[3001:4000,str_detect(colnames(Scenario4$Nimble1),'tilde_gamma')], 1, function(x) rdirichlet(1,x)[,2]),
             apply(Scenario4$Nimble2[3001:4000,str_detect(colnames(Scenario4$Nimble2),'tilde_gamma')], 1, function(x) rdirichlet(1,x)[,2]))


ggdata3 = data.frame(m1= membership_scores[,1],
                     m2 = membership_scores[,2],
                     l1 = tmp3[,1],
                     l2 = tmp3[,2]) 


plot_scenario4_NIMBLE = 
                       ggplot(ggdata3) +
                       geom_polygon(aes(x,y),data = data.frame(x=c(0,1,1,0),y=c(0,0,1,1))
                       ,fill = "grey50",alpha = 0.3 ) + 
                       geom_point(aes(m1,m2), alpha = 0.6,size = 0.8) +
                       geom_point(aes(l1,l2),col = '#00249C',shape = 4, alpha = 0.6,size = 0.8) +
                       scale_x_continuous(breaks = seq(0,1,0.25),labels = sprintf("%.2f" ,seq(0,1,0.25)),expand=c(0.009,0.02))+
                       scale_y_continuous(breaks = seq(0.25,1,0.25),labels = sprintf("%.2f" ,seq(0.25,1,0.25)),expand=c(0.001,0.001))+
                       facet_wrap(~I('SCENARIO 4'))+
                       labs(x  = expression(hat(lambda)^{(1)}),y =  expression(hat(lambda)^{(2)})) + 
                       theme_bw()
```

### Merge plots

Following code produces Figure S1 of the Supplementary Material
displayed below

``` r
tot_plot_NIMBLE = arrangeGrob(plot_scenario1_NIMBLE,plot_scenario2_NIMBLE,plot_scenario3_NIMBLE,plot_scenario4_NIMBLE)

ggsave(plot = tot_plot_NIMBLE,file = "plots/All_membership_distr_nimble.png",height = 15,width =20,units = 'cm',dpi = 300)
```

![](plots/All_membership_distr_nimble.png)

MSE table
---------

``` r
MSE = tibble()

## SCENARIO 1
membership_scores = apply(Scenario1$lambda[2501:5000,,],c(2,3), mean)
nimble_scores =  cbind(colMeans(Scenario1$Nimble1[-c(1:1500), str_detect(colnames(Scenario1$Nimble1),"lambda.[0-9]{1,4}[, ] 2")]),  
colMeans(Scenario1$Nimble2[-c(1:1500), str_detect(colnames(Scenario1$Nimble2),"lambda.[0-9]{1,4}[, ] 2")]))

err2 =  (membership_scores[,1] - Scenario1$true_parameter$sbj_prob[,1])^2
MSE[1,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 =  (membership_scores[,2] - Scenario1$true_parameter$sbj_prob[,2])^2
MSE[2,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))
         
err2 = (Scenario1$lambda.point1[,2] - Scenario1$true_parameter$sbj_prob[,1])^2 
MSE[3,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = (Scenario1$lambda.point2[,2] - Scenario1$true_parameter$sbj_prob[,2])^2 
MSE[4,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = (nimble_scores[,1] - Scenario1$true_parameter$sbj_prob[,1])^2
MSE[5,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = ( nimble_scores[,2] - Scenario1$true_parameter$sbj_prob[,2] )^2
MSE[6,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


MSE = cbind(` `=c('MMM g = 1','MMM g = 2','mixedMem g = 1','mixedMem g = 2',
                   'MCMC g = 1', ' MCMC g = 2'),MSE)

names(MSE)[2]  = "SCENARIO 1"

## SCENARIO 2

membership_scores = apply(Scenario2$lambda[2501:5000,,],c(2,3),mean)
nimble_scores =  cbind(colMeans(Scenario2$Nimble1[-c(1:1500), str_detect(colnames(Scenario2$Nimble1),"lambda.[0-9]{1,4}[, ] 2")]),  
colMeans(Scenario2$Nimble2[-c(1:1500), str_detect(colnames(Scenario2$Nimble2),"lambda.[0-9]{1,4}[, ] 2")])) 



err2 =  (membership_scores[,1] - Scenario2$true_parameter$sbj_prob[,1])^2
MSE[1,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 =  (membership_scores[,2] - Scenario2$true_parameter$sbj_prob[,2])^2
MSE[2,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))
         
err2 = (Scenario2$lambda.point1[,2] - Scenario2$true_parameter$sbj_prob[,1])^2 
MSE[3,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = (Scenario2$lambda.point2[,2] - Scenario2$true_parameter$sbj_prob[,2])^2 
MSE[4,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = (nimble_scores[,1]- Scenario2$true_parameter$sbj_prob[,1])^2
MSE[5,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = (nimble_scores[,2]- Scenario2$true_parameter$sbj_prob[,2])^2
MSE[6,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

names(MSE)[3]  = "SCENARIO 2"

## SCENARIO 3

membership_scores = apply(Scenario3$lambda[2501:5000,,],c(2,3),mean)
err2 =  (membership_scores[,1] - Scenario3$true_parameter$sbj_prob[,1])^2

nimble_scores =  cbind(colMeans(Scenario3$Nimble1[-c(1:1500), str_detect(colnames(Scenario3$Nimble1),"lambda.[0-9]{1,4}[, ] 2")]),  
colMeans(Scenario3$Nimble2[-c(1:1500), str_detect(colnames(Scenario3$Nimble2),"lambda.[0-9]{1,4}[, ] 2")]))




MSE[1,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 =  (membership_scores[,2] - Scenario3$true_parameter$sbj_prob[,2])^2
MSE[2,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))
         
err2 = (Scenario3$lambda.point1[,2] - Scenario3$true_parameter$sbj_prob[,1])^2 
MSE[3,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = (Scenario3$lambda.point2[,2] - Scenario3$true_parameter$sbj_prob[,2])^2 
MSE[4,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))



err2 = (nimble_scores[,1] - Scenario3$true_parameter$sbj_prob[,1])^2 
MSE[5,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


err2 = (nimble_scores[,2] - Scenario3$true_parameter$sbj_prob[,2])^2
MSE[6,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


names(MSE)[4]  = "SCENARIO 3"

## SCENARIO 4


membership_scores = apply(Scenario4$lambda[2501:5000,,],c(2,3),mean)
nimble_scores =  cbind(colMeans(Scenario4$Nimble1[-c(1:1500), str_detect(colnames(Scenario4$Nimble1),"lambda.[0-9]{1,4}[, ] 2")]),  
colMeans(Scenario4$Nimble2[-c(1:1500), str_detect(colnames(Scenario4$Nimble2),"lambda.[0-9]{1,4}[, ] 2")]))


err2 =  (membership_scores[,1] - Scenario4$true_parameter$sbj_prob[,1])^2
MSE[1,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 =  (membership_scores[,2] - Scenario4$true_parameter$sbj_prob[,2])^2
MSE[2,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))
         
err2 = (Scenario4$lambda.point1[,2] - Scenario4$true_parameter$sbj_prob[,1])^2 
MSE[3,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = (Scenario4$lambda.point2[,2] - Scenario4$true_parameter$sbj_prob[,2])^2 
MSE[4,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


err2 = (nimble_scores[,1] - Scenario4$true_parameter$sbj_prob[,1])^2 
MSE[5,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = (nimble_scores[,2] - Scenario4$true_parameter$sbj_prob[,2])^2 
MSE[6,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


names(MSE)[5]  = "SCENARIO 4"


## generate latex table
MSE  %>% knitr::kable(.,booktab = 3,format = 'latex') %>% cat(.,file = 'plots/MSE_sim_nimble.tex')
```

The resulting table is

|                | SCENARIO 1   | SCENARIO 2   | SCENARIO 3   | SCENARIO 4   |
|:---------------|:-------------|:-------------|:-------------|:-------------|
| MMM g = 1      | 0.027(0.035) | 0.025(0.040) | 0.023(0.032) | 0.037(0.045) |
| MMM g = 2      | 0.026(0.034) | 0.029(0.042) | 0.023(0.036) | 0.030(0.044) |
| mixedMem g = 1 | 0.029(0.035) | 0.034(0.051) | 0.043(0.050) | 0.043(0.049) |
| mixedMem g = 2 | 0.030(0.041) | 0.040(0.054) | 0.033(0.047) | 0.031(0.042) |
| MCMC g = 1     | 0.078(0.085) | 0.029(0.046) | 0.036(0.048) | 0.037(0.047) |
| MCMC g = 2     | 0.032(0.043) | 0.035(0.054) | 0.035(0.055) | 0.031(0.046) |

L1-Norm table
-------------

``` r
L1 = tibble()

## SCENARIO 1
membership_scores = apply(Scenario1$lambda[2501:5000,,],c(2,3), mean)
err2 =  abs(membership_scores[,1] - Scenario1$true_parameter$sbj_prob[,1])
L1[1,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 =  abs(membership_scores[,2] - Scenario1$true_parameter$sbj_prob[,2])
L1[2,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))
         
err2 = abs(Scenario1$lambda.point1[,2] - Scenario1$true_parameter$sbj_prob[,1])
L1[3,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(Scenario1$lambda.point2[,2] - Scenario1$true_parameter$sbj_prob[,2])
L1[4,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


err2 = abs(colMeans(Scenario1$Nimble1[-c(1:2500), str_detect(colnames(Scenario1$Nimble1),"lambda.[0-9]{1,4}[, ] 2")]) 
           - Scenario1$true_parameter$sbj_prob[,1] )
L1[5,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(colMeans(Scenario1$Nimble2[-c(1:2500), str_detect(colnames(Scenario1$Nimble2),"lambda.[0-9]{1,4}[, ] 2")])
           - Scenario1$true_parameter$sbj_prob[,2] )
L1[6,1] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))



L1 = cbind(` `=c('MMM g = 1','MMM g = 2','mixedMem g = 1','mixedMem g = 2', 'MCMC g = 1', ' MCMC g = 2'),L1)

names(L1)[2]  = "SCENARIO 1"
## SCENARIO 2

membership_scores = apply(Scenario2$lambda[2501:5000,,],c(2,3),mean)

err2 =  abs(membership_scores[,1] - Scenario2$true_parameter$sbj_prob[,1])
L1[1,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 =  abs(membership_scores[,2] - Scenario2$true_parameter$sbj_prob[,2])
L1[2,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))
         
err2 = abs(Scenario2$lambda.point1[,2] - Scenario2$true_parameter$sbj_prob[,1])
L1[3,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(Scenario2$lambda.point2[,2] - Scenario2$true_parameter$sbj_prob[,2])
L1[4,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(colMeans(Scenario2$Nimble1[-c(1:2500), str_detect(colnames(Scenario2$Nimble1),"lambda.[0-9]{1,4}[, ] 2")]) - Scenario2$true_parameter$sbj_prob[,1] )
L1[5,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(colMeans(Scenario2$Nimble2[-c(1:2500), str_detect(colnames(Scenario2$Nimble2),"lambda.[0-9]{1,4}[, ] 2")])
          - Scenario2$true_parameter$sbj_prob[,2] 
           )
L1[6,3] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


names(L1)[3]  = "SCENARIO 2"

## SCENARIO 3

membership_scores = apply(Scenario3$lambda[2501:5000,,],c(2,3),mean)
err2 =  abs(membership_scores[,1] - Scenario3$true_parameter$sbj_prob[,1])
L1[1,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 =  abs(membership_scores[,2] - Scenario3$true_parameter$sbj_prob[,2])
L1[2,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))
         
err2 = abs(Scenario3$lambda.point1[,2] - Scenario3$true_parameter$sbj_prob[,1])
L1[3,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(Scenario3$lambda.point2[,2] - Scenario3$true_parameter$sbj_prob[,2])
L1[4,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


err2 = abs(colMeans(Scenario3$Nimble1[-c(1:2500), str_detect(colnames(Scenario3$Nimble1),"lambda.[0-9]{1,4}[, ] 2")])  - Scenario3$true_parameter$sbj_prob[,1])
L1[5,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(colMeans(Scenario3$Nimble2[-c(1:2500), str_detect(colnames(Scenario3$Nimble2),"lambda.[0-9]{1,4}[, ] 2")]) - Scenario3$true_parameter$sbj_prob[,2])
L1[6,4] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

names(L1)[4]  = "SCENARIO 3"

## SCENARIO 4

membership_scores = apply(Scenario4$lambda[2501:5000,,],c(2,3),mean)
err2 =  abs(membership_scores[,1] - Scenario4$true_parameter$sbj_prob[,1])
L1[1,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 =  abs(membership_scores[,2] - Scenario4$true_parameter$sbj_prob[,2])
L1[2,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))
         
err2 = abs(Scenario4$lambda.point1[,2] - Scenario4$true_parameter$sbj_prob[,1])
L1[3,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(Scenario4$lambda.point2[,2] - Scenario4$true_parameter$sbj_prob[,2])
L1[4,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))


err2 = abs(colMeans(Scenario4$Nimble1[-c(1:2500), str_detect(colnames(Scenario4$Nimble1),"lambda.[0-9]{1,4}[, ] 2")]) - Scenario4$true_parameter$sbj_prob[,1])
L1[5,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

err2 = abs(colMeans(Scenario4$Nimble2[-c(1:2500), str_detect(colnames(Scenario4$Nimble2),"lambda.[0-9]{1,4}[, ] 2")])- Scenario4$true_parameter$sbj_prob[,2])
L1[6,5] = sprintf("%.3f(%.3f)", mean(err2),sd(err2))

names(L1)[5]  = "SCENARIO 4"


## generate latex table
L1 %>% knitr::kable(.,booktab = 3,format = 'latex') %>% cat(.,file = 'plots/L1_sim.tex')
```

The resulting table is

|                | SCENARIO 1   | SCENARIO 2   | SCENARIO 3   | SCENARIO 4   |
|:---------------|:-------------|:-------------|:-------------|:-------------|
| MMM g = 1      | 0.132(0.096) | 0.126(0.097) | 0.122(0.090) | 0.162(0.106) |
| MMM g = 2      | 0.130(0.094) | 0.134(0.103) | 0.117(0.095) | 0.138(0.105) |
| mixedMem g = 1 | 0.139(0.096) | 0.148(0.110) | 0.174(0.113) | 0.174(0.113) |
| mixedMem g = 2 | 0.140(0.104) | 0.162(0.119) | 0.147(0.108) | 0.141(0.103) |
| MCMC g = 1     | 0.233(0.150) | 0.131(0.106) | 0.156(0.111) | 0.157(0.110) |
| MCMC g = 2     | 0.153(0.111) | 0.147(0.118) | 0.147(0.117) | 0.139(0.109) |

Kernel plots
------------

### Plot summaries

#### Scenario 1

``` r
pdata1 = 
Scenario1$kern[2501:5000,,,] %>% melt %>% rename(samp = Var1,d = Var2, j = Var3, h = Var4) %>% group_by(d,j,h) %>% summarize(MMM = median(value),
                                         post_low = quantile(value,0.1),
                                         post_up = quantile(value,0.9)) 
pdata1 = pdata1 %>%
left_join(
Scenario1$true_parameter$kern %>% melt %>% rename( d = Var1, j = Var2, h = Var3,true = value))




tmp_nimble1 = Scenario1$Nimble1[-c(1:1500), str_subset(colnames(Scenario1$Nimble1),"^theta")] %>% melt 

tmp_nimble1$Var2 = str_extract(tmp_nimble1$Var2,"[1-4], [1-5], [1-2]")
tmp2 = do.call("rbind", str_split(tmp_nimble1$Var2,","))
tmp_nimble1$d =  as.numeric(tmp2[,1])
tmp_nimble1$j =  as.numeric(tmp2[,2])
tmp_nimble1$h =  as.numeric(tmp2[,3]) 


tmp_nimble2 = Scenario1$Nimble2[-c(1:1500), str_subset(colnames(Scenario1$Nimble2),"^theta")] %>% melt 
tmp_nimble2$Var2 = str_extract(tmp_nimble2$Var2,"[1-4], [1-5], [1-2]")

tmp2 = do.call("rbind", str_split(tmp_nimble2$Var2,","))
tmp_nimble2$d =  as.numeric(tmp2[,1])
tmp_nimble2$j =  as.numeric(tmp2[,2]) + 5
tmp_nimble2$h =  as.numeric(tmp2[,3]) 


pdata1 = 
pdata1 %>% left_join(
                    tmp_nimble1 %>%
                      group_by(d,j,h) %>%
                      summarize(MCMC     = median(value),
                                MCMC_low = quantile(value,prob=0.1),
                                MCMC_up  = quantile(value,prob=0.9)) %>%
                      bind_rows(
                                tmp_nimble2 %>%
                                  group_by(d,j,h) %>%
                                  summarize(MCMC     = median(value),
                                            MCMC_low = quantile(value,prob=0.1),
                                            MCMC_up  = quantile(value,prob=0.9)))
)
```

#### Scenario 2

``` r
pdata2 = 
Scenario2$kern[2501:5000,,,] %>% melt %>% rename(samp = Var1,d = Var2, j = Var3, h = Var4) %>% group_by(d,j,h) %>% summarize(MMM = median(value),
                                         post_low = quantile(value,0.1),
                                         post_up = quantile(value,0.9)) 
pdata2 = pdata2 %>%
left_join(
Scenario2$true_parameter$kern %>% melt %>% rename( d = Var1, j = Var2, h = Var3,true = value))




tmp_nimble1 = Scenario2$Nimble1[-c(1:1500), str_subset(colnames(Scenario2$Nimble1),"^theta")] %>% melt 

tmp_nimble1$Var2 = str_extract(tmp_nimble1$Var2,"[1-4], [1-5], [1-2]")
tmp2 = do.call("rbind", str_split(tmp_nimble1$Var2,","))
tmp_nimble1$d =  as.numeric(tmp2[,1])
tmp_nimble1$j =  as.numeric(tmp2[,2])
tmp_nimble1$h =  as.numeric(tmp2[,3]) 


tmp_nimble2 = Scenario2$Nimble2[-c(1:1500), str_subset(colnames(Scenario2$Nimble2),"^theta")] %>% melt 
tmp_nimble2$Var2 = str_extract(tmp_nimble2$Var2,"[1-4], [1-5], [1-2]")

tmp2 = do.call("rbind", str_split(tmp_nimble2$Var2,","))
tmp_nimble2$d =  as.numeric(tmp2[,1])
tmp_nimble2$j =  as.numeric(tmp2[,2]) + 5
tmp_nimble2$h =  as.numeric(tmp2[,3]) 


pdata2 = 
pdata2 %>% left_join(
                    tmp_nimble1 %>%
                      group_by(d,j,h) %>%
                      summarize(MCMC     = median(value),
                                MCMC_low = quantile(value,prob=0.1),
                                MCMC_up  = quantile(value,prob=0.9)) %>%
                      bind_rows(
                                tmp_nimble2 %>%
                                  group_by(d,j,h) %>%
                                  summarize(MCMC     = median(value),
                                            MCMC_low = quantile(value,prob=0.1),
                                            MCMC_up  = quantile(value,prob=0.9)))
)
```

#### Scenario 3

``` r
pdata3 = 
Scenario3$kern[2501:5000,,,] %>% melt %>% rename(samp = Var1,d = Var2, j = Var3, h = Var4) %>% group_by(d,j,h) %>% summarize(MMM = median(value),
                                         post_low = quantile(value,0.1),
                                         post_up = quantile(value,0.9)) 
pdata3 = pdata3 %>%
left_join(
Scenario3$true_parameter$kern %>% melt %>% rename( d = Var1, j = Var2, h = Var3,true = value))





tmp_nimble1 = Scenario3$Nimble1[-c(1:1500), str_subset(colnames(Scenario3$Nimble1),"^theta")] %>% melt 

tmp_nimble1$Var2 = str_extract(tmp_nimble1$Var2,"[1-4], [1-5], [1-2]")
tmp2 = do.call("rbind", str_split(tmp_nimble1$Var2,","))
tmp_nimble1$d =  as.numeric(tmp2[,1])
tmp_nimble1$j =  as.numeric(tmp2[,2])
tmp_nimble1$h =  as.numeric(tmp2[,3]) 


tmp_nimble2 = Scenario3$Nimble2[-c(1:1500), str_subset(colnames(Scenario3$Nimble2),"^theta")] %>% melt 
tmp_nimble2$Var2 = str_extract(tmp_nimble2$Var2,"[1-4], [1-5], [1-2]")

tmp2 = do.call("rbind", str_split(tmp_nimble2$Var2,","))
tmp_nimble2$d =  as.numeric(tmp2[,1])
tmp_nimble2$j =  as.numeric(tmp2[,2]) + 5
tmp_nimble2$h =  as.numeric(tmp2[,3]) 


pdata3 = 
pdata3 %>% left_join(
                    tmp_nimble1 %>%
                      group_by(d,j,h) %>%
                      summarize(MCMC     = median(value),
                                MCMC_low = quantile(value,prob=0.1),
                                MCMC_up  = quantile(value,prob=0.9)) %>%
                      bind_rows(
                                tmp_nimble2 %>%
                                  group_by(d,j,h) %>%
                                  summarize(MCMC     = median(value),
                                            MCMC_low = quantile(value,prob=0.1),
                                            MCMC_up  = quantile(value,prob=0.9)))
)
```

#### Scenario 4

``` r
pdata4 = 
Scenario4$kern[2501:5000,,,] %>% melt %>% rename(samp = Var1,d = Var2, j = Var3, h = Var4) %>% group_by(d,j,h) %>% summarize(MMM = median(value),
                                         post_low = quantile(value,0.1),
                                         post_up = quantile(value,0.9)) 
pdata4 = pdata4 %>%
left_join(
Scenario4$true_parameter$kern %>% melt %>% rename( d = Var1, j = Var2, h = Var3,true = value))



tmp_nimble1 = Scenario4$Nimble1[-c(1:1500), str_subset(colnames(Scenario4$Nimble1),"^theta")] %>% melt 

tmp_nimble1$Var2 = str_extract(tmp_nimble1$Var2,"[1-4], [1-5], [1-2]")
tmp2 = do.call("rbind", str_split(tmp_nimble1$Var2,","))
tmp_nimble1$d =  as.numeric(tmp2[,1])
tmp_nimble1$j =  as.numeric(tmp2[,2])
tmp_nimble1$h =  as.numeric(tmp2[,3]) 


tmp_nimble2 = Scenario4$Nimble2[-c(1:1500), str_subset(colnames(Scenario4$Nimble2),"^theta")] %>% melt 
tmp_nimble2$Var2 = str_extract(tmp_nimble2$Var2,"[1-4], [1-5], [1-2]")

tmp2 = do.call("rbind", str_split(tmp_nimble2$Var2,","))
tmp_nimble2$d =  as.numeric(tmp2[,1])
tmp_nimble2$j =  as.numeric(tmp2[,2]) + 5
tmp_nimble2$h =  as.numeric(tmp2[,3]) 


pdata4 = 
pdata4 %>% left_join(
                    tmp_nimble1 %>%
                      group_by(d,j,h) %>%
                      summarize(MCMC     = median(value),
                                MCMC_low = quantile(value,prob=0.1),
                                MCMC_up  = quantile(value,prob=0.9)) %>%
                      bind_rows(
                                tmp_nimble2 %>%
                                  group_by(d,j,h) %>%
                                  summarize(MCMC     = median(value),
                                            MCMC_low = quantile(value,prob=0.1),
                                            MCMC_up  = quantile(value,prob=0.9)))
)
```

#### Total plots

``` r
pdata_tot = bind_rows(cbind.data.frame(pdata1, type = 'SCENARIO 1'),
                      cbind.data.frame(pdata2, type = 'SCENARIO 2'),
                      cbind.data.frame(pdata3, type = 'SCENARIO 3'),
                      cbind.data.frame(pdata4, type = 'SCENARIO 4')
                  )

pdata_tot$h =factor(ifelse(pdata_tot$h ==1,'PROFILE 1', 'PROFILE 2'))

tmp_data = pdata_tot  %>%
           filter(j==4) %>%
           dplyr::select(-j) %>% 
           gather(method,value,-d,-h,-post_low,-post_up,-type) 



tmp_data[tmp_data$method=='MCMC',"post_up"] = tmp_data[tmp_data$method=='MCMC_up',"value"]
tmp_data[tmp_data$method=='MCMC',"post_low"] = tmp_data[tmp_data$method=='MCMC_low',"value"]


tmp_data[tmp_data$method=='MM',"post_low"] = tmp_data[tmp_data$method=='lowMM',"value"]
tmp_data[tmp_data$method=='MM',"post_up"] = tmp_data[tmp_data$method=='upMM',"value"]




tmp_data = tmp_data %>% dplyr::filter(method %in% c('MMM','MM',"true","MCMC"))
tmp_data$method = tmp_data$method %>% toupper %>% factor(., levels = c('MMM','TRUE','MM','MCMC'))


###tmp_data$variable = tmp_data$variable %>% as.character %>% toupper %>% factor(., levels = c('MMM','TRUE','MM'))
### 
tmp_data$post_low[tmp_data$method=='TRUE'] = NA
tmp_data$post_up[tmp_data$method=='TRUE']  = NA
```

Finally we can plot the results

``` r
## comparison with NIMBLE MCMC
tmp_data %>%
 filter(method!='MM')%>%
 ggplot +
 geom_bar(aes(x = d, y = value,fill = method),position = position_dodge(),stat = 'identity',col = 'grey70') +
 geom_errorbar(aes(x = d, ymin = post_low, ymax = post_up,group=method), width=.5,
                 position=position_dodge(.9) ) +
scale_fill_grey(start = 0.5, end = 0.8)+
facet_grid(factor(h)~type)+ 
labs(x = '',y = '', fill = '') + 
theme_bw()
ggsave(file = 'plots/kern_group1_nimble.png',height = 15,width = 30,units = 'cm')


tmp_data = pdata_tot  %>%
           filter(j==4) %>%
           dplyr::select(-j) %>% 
           gather(method,value,-d,-h,-post_low,-post_up,-type) 
```

Previous code produces the figure

![](plots/kern_group1_nimble.png)

We do the same for a variable in group 2

``` r
### GROUP 2 variable
tmp_data = pdata_tot  %>%
           filter(j==8) %>%
           dplyr::select(-j) %>% 
           gather(method,value,-d,-h,-post_low,-post_up,-type) 


tmp_data[tmp_data$method=='MCMC',"post_up"] = tmp_data[tmp_data$method=='MCMC_up',"value"]
tmp_data[tmp_data$method=='MCMC',"post_low"] = tmp_data[tmp_data$method=='MCMC_low',"value"]

tmp_data[tmp_data$method=='MM',"post_low"] = tmp_data[tmp_data$method=='lowMM',"value"]
tmp_data[tmp_data$method=='MM',"post_up"] = tmp_data[tmp_data$method=='upMM',"value"]


tmp_data = tmp_data %>% dplyr::filter(method %in% c('MMM','MM',"true","MCMC"))
tmp_data$method = tmp_data$method %>% toupper %>% factor(., levels = c('MMM','TRUE','MM','MCMC'))
tmp_data$post_low[tmp_data$method=='TRUE'] = NA
tmp_data$post_up[tmp_data$method=='TRUE']  = NA

tmp_data %>% 
 filter(method!='MM')%>%
 ggplot +
 geom_bar(aes(x = d, y = value,fill = method),position = position_dodge(),stat = 'identity',col = 'grey70') +
 geom_errorbar(aes(x = d, ymin = post_low, ymax = post_up,group=method), width=.5,
                 position=position_dodge(.9) ) +
scale_fill_grey(start = 0.5, end = 0.8)+
facet_grid(factor(h)~type)+ 
labs(x = '',y = '', fill = '') + 
theme_bw()
ggsave(file = 'plots/kern_group2_nimble.png',height = 15,width = 30,units = 'cm')
```

Having as result

![](plots/kern_group2_nimble.png)

Misspecification: more than two pure types
------------------------------------------

``` r
Scenario_miss = c(loadRData('results/Sim_miss_MMM_results.Rdata'),
                  loadRData('results/Sim_miss_MM_results.Rdata'))

Scenario_miss$Nimble =  readRDS("results/Sim_nimble_model_miss_model_1.rds")

load('simulated_data/data_and_true_parameters_misspecified.Rdata')
Scenario_miss$data  = data
Scenario_miss$true_parameter  = true_parameter
rm(data,true_parameter)


membership_score = apply(Scenario_miss$lambda[2501:5000,,],c(2,3), mean)
nimble_score = cbind(colMeans(Scenario_miss$Nimble[-c(1:1500), str_detect(colnames(Scenario_miss$Nimble),"lambda.[0-9]{1,4}[, ] 2")]),  
colMeans(Scenario_miss$Nimble[-c(1:1500), str_detect(colnames(Scenario_miss$Nimble),"lambda.[0-9]{1,4}[, ] 2")]))

pdata = cbind.data.frame(nimble_score[,1],membership_score[,1], Scenario_miss$true_parameter$sbj_prob[[1]]) 

names(pdata) = c('MM-MCMC','MMM','PROFILE 1' ,'PROFILE 2','PROFILE 3', 'PROFILE 4')


 pdata %>% melt(c(1,2)) %>% melt(c(3,4)) %>% setNames(c('truth','y','model','x')) %>%  
    ggplot +
    geom_point(aes(x = x, y = y,shape = model,col = model)) + 
    scale_shape(solid = FALSE) + 
    scale_colour_brewer(palette='Set1')+
    scale_x_continuous(limits = c(0,1),expand = c(-0.1,0)) + 
    scale_y_continuous(limits = c(0,1)) +
    facet_wrap(~truth) + 
    geom_hline(yintercept = 0.5, lty = 'dashed') +
    labs(x = expression(hat(lambda)^{(1)}),y = '',col ='',shape = '')+
    theme_bw()
ggsave(file = 'plots/multi_prof_nimble.png',height = 15,width = 30,units = 'cm')
```

![](plots/multi_prof_nimble.png)

### Mixed profiles

We can obtain the elements of Table S2 in the Supplementary Material
with the following code

``` r
## MMM variable 1 profile 1 kernel
Scenario_miss$kern[2501:5000,,1,1] %>%  apply(2, function(x) sprintf("%.3f (%.3f;%.3f)", mean(x),quantile(x,0.1),quantile(x,0.9)))

## MMM variable 1 profile 2 kernel
Scenario_miss$kern[2501:5000,,1,2] %>%  apply(2, function(x) sprintf("%.3f (%.3f;%.3f)", mean(x),quantile(x,0.1),quantile(x,0.9)))

## MM variable 1 profile 1 kernel
Scenario_miss$Nimble[-c(1:1500),   str_subset(colnames(Scenario_miss$Nimble),"^theta\\[[1-4], 1, 1")]                  %>% apply(2, function(x) sprintf("%.3f (%.3f;%.3f)", mean(x),quantile(x,0.1),quantile(x,0.9)))

# MM variable 1 profile 2 kernel
Scenario_miss$Nimble[-c(1:1500),   str_subset(colnames(Scenario_miss$Nimble),"^theta\\[[1-4], 1, 2")]                  %>% apply(2, function(x) sprintf("%.3f (%.3f;%.3f)", mean(x),quantile(x,0.1),quantile(x,0.9)))
```
