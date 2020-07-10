## MM model using MCMC to approximate the posterior
## see for example:
## Erosheva, Elena A.; Fienberg, Stephen E.; Joutard, Cyrille. Describing disability through individual-level mixture models for multivariate binary data. Ann. Appl. Stat. 1 (2007), no. 2, 502--537. doi:10.1214/07-AOAS126. https://projecteuclid.org/euclid.aoas/1196438029



library(nimble)
#'@param d: number of categories for each variable
#'@param N: number of subjects
#'@param p: number of variables
#'@param H: dimension of the latent space
#'@param a_gamma, a_gamma: hyperprior for the gamma prior on Dirichlet concentration parameters

MM_nimble_code = nimbleCode(
{
 
  for(i in 1:N)
  {
    for(j in 1:p)
    {
     X[i,j] ~ dcat(theta[1:d,j,Z[i,j]])
     Z[i,j] ~ dcat(lambda[i,1:H])
    }
    lambda[i,1:H] ~ ddirch(tilde_gamma[1:H])
  }

  xi ~ dgamma(a_gamma,rate = b_gamma)
  gamma[1:H] ~ ddirch(alpha_gamma[1:H]) 
  tilde_gamma[1:H]  <- gamma[1:H] * xi

  for(j in 1:p)
  {
    for(h in 1:H)
    {
     theta[1:d,j,h] ~ ddirch(alpha[1:d])
    }
}
}
)

