d2MLND = function(x,mu,Sigma)
{
  logitx = log(x/(1-x))
  
out = 
  1/(sqrt(2*pi*det(Sigma))*prod(x)*prod(1-x))*
  exp( -0.5*( t(logitx  - mu) %*%solve(Sigma) %*% (logitx  - mu)))
return(out)
}
