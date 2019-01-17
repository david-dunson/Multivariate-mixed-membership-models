r2MLND = function(n,mu,Sigma)
{
 res = apply(rmvn(n, mu = mu, sigma = Sigma),2,plogis)
 return(res)
}


