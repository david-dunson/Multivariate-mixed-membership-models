estimate_MMM  = function(data,    ## list containing: n = number of subjects; p = number of variables; gj = group indicators; data nxp matrix
                 prior=list(),    ## list having mu; Sigma; nu; Psi
                 param = list(),  ##optional list of parameter initialization
                 nrep     = 5000) ## number of MCMC samples 
{
##########################################
##########################################
## REQUIRED LIBRARIES 
##########################################
##########################################
suppressMessages(require(mvnfast))
suppressMessages(require(MCMCpack))

##########################################
##########################################
## Model  parameters 
##########################################
##########################################


if(length(param) ==0)
{
 param$omega  = matrix(NA,data$n,2) ## data augmented polyagamma values
 param$psi    = matrix(0,data$n,2)
 param$ki     = matrix(NA,data$n,2)
 param$mu     = c(0,0)
 param$lambda = matrix(0.5,data$n,2)
 param$kern   = array(NA,dim  = c(4,data$p,2))
 param$zi     = matrix(sample(1:2,replace = T, size = data$n*data$p),ncol = data$p)
 param$Sigma  =  diag(1,2)
}
##########################################
##########################################
## Prior
##########################################
##########################################

if(length(prior) ==0)
{
  prior = list()
  prior$mu     = c(0,0)
  prior$nu     = 2
  prior$Psi    = diag(1,2)
  prior$Sigma  = diag(1,2)
}


## initial values for variance and convariance matrix; set prior
  Sigma    = prior$Sigma
  invSigma = solve(Sigma)




##########################################
##########################################
# Result place holder
##########################################
##########################################

results  =list( 
              mu      = matrix(NA,nrep,2),
              psi     = array(NA,dim = c(nrep,data$n,2)),
              lambda  = array(NA,dim = c(nrep,data$n,2)),
              Sigma   = array(NA,dim = c(nrep,2,2)),
              kern    = array(NA,dim = c(nrep,4,data$p,2)),
              zi      = array(NA,dim = c(nrep,data$n,data$p)))



##########################################
##########################################
## ALGORITHM
##########################################
##########################################

#define a progress bar
progress <- txtProgressBar(min = 0, max = nrep, style = 3)

for(r in 1:nrep){

##########################################
##########################################
# update the kernels
##########################################
##########################################
for(j in 1:data$p)
{
 for(h in 1:2)
 {
     param$kern[,j,h] = c(rdirichlet(1,table(factor(data$data[param$zi[,j]==h,j],1:4)) +1/4))
 }
}
    
##########################################
##########################################
## individual - latent analysis analysis
##########################################
##########################################

## some parameters are related to time varing analysis and not acltually needed here
## 

param$eta = t(sapply(1:data$n,function(i)  param$mu))

## reshape kernels in list form
tmp_kern = list()
for(j in 1:data$p){
tmp_kern[[j]] = param$kern[,j,]
}

attr(invSigma,'dim') =c(2,2,1) ## convert to array

tmp_latent = .update_subject_latent( 
                                    data                    = list(data$data),
                                    Kern                   = tmp_kern,
                                    ZI                     = list(param$zi),
                                    ki                     = param$ki,
                                    eta                    = param$eta,
                                    psi                    = param$psi,
                                    lambda                 = param$lambda,
                                    omega                  = param$omega,
                                    ISigma_psi             = invSigma,
                                    subject_year_indicator = rep(1,data$n),
                                    g_j                    = list(data$gj),
                                    i_t                    = 1:data$n,
                                    inv_j_t                = list(1:data$p))

param$ki     = tmp_latent[['ki']]
param$psi    = tmp_latent[['psi']]
param$lambda = tmp_latent[['lambda']]
param$omega  = tmp_latent[['omega']]
param$zi     = tmp_latent[['ZI']][[1]] + 1


##########################################
##########################################
## update general mean
##########################################
##########################################

ytilde  = t(sapply(1:data$n,function(i) param$ki[i,]/param$omega[i,]))
 
ISigma_data = apply(
 simplify2array(lapply(1:data$n,function(i) solve(diag(1/param$omega[i,]) + param$Sigma))),c(1,2),sum)

W = apply(simplify2array(lapply(1:data$n,function(i) solve(diag(1/param$omega[i,]) + param$Sigma)%*%ytilde[i,] ) ),c(1,2),sum)
V_omega = solve(solve(prior$Sigma) + ISigma_data)
m_omega = c(V_omega%*%(solve(prior$Sigma)%*%prior$mu + W))

param$mu  = c(rmvn(1, m_omega, V_omega))


##########################################
##########################################
## Sample sigma
##########################################
##########################################

bar_psi          = colMeans(param$psi)
V                = crossprod(param$psi,param$psi) + data$n*tcrossprod(param$mu,param$mu) - data$n*tcrossprod(param$mu,bar_psi) -data$n*tcrossprod(bar_psi,param$mu)


param$Sigma      = riwish(data$n+ prior$nu ,  V + prior$Psi)
invSigma         = solve(param$Sigma)


##########################################
##########################################
## save the results across each iteration
##########################################
##########################################

results$mu[r,]        = param$mu
results$Sigma[r,,]    = param$Sigma
results$psi[r,,]      = param$psi
results$lambda[r,,]   = param$lambda
results$kern[r,,,]    = param$kern
results$zi[r,,]       = param$zi

## update progress-bar
setTxtProgressBar(progress, r)
}
## start a new line in the prompt 
cat('\n')
## return results
results$psi = NULL
return(results)
}
