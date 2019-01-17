#include "RcppArmadillo.h"
#include <R_ext/Utils.h>
#include <iostream>
#include <exception>
#include "RNG.h"
#include "PolyaGamma.h"
// Rcpp::plugins(cpp11)
// Rcpp::depends(RcppArmadillo)

using namespace Rcpp;


/*
 *
 * This code implement the subject specific part on Algorithm 1 leveraging C++
 * To include polyagamma data augmentation needed from the algorithm we use Scott package skeleton publicly available at
 * https://github.com/jgscott/helloPG
 */


//////////////////////////////////////////
//////////////////////////////////////////
// Utilities (for internal use)
//////////////////////////////////////////
//////////////////////////////////////////



/*
 * Equivalent to R plogis function
 */
arma::vec plogis(arma::vec x){
    arma::vec retval(x.n_elem);
    for(arma::uword i=0;i<x.n_elem;i++){
    retval(i) = exp(x(i))/(1 + exp(x(i)));
    }
    return retval;
}


/*
 *  This function "normalize" the log-vector (i.e.  exp(arma::vec)/sum(exp(arma::vec)) ), leveraging logarithms
 */

arma::vec log_vec_renorm(arma::vec x)
{
  arma::vec  retval(x.n_elem);
    for(arma::uword i =0; i<x.n_elem;i++)
    {
      retval(i) =  exp(-log(sum(exp(x - x(i)))));
    }
  return retval;
}


/*
 * multivariate normal generation. Original code from Ahmadou Dicko available at  
 *http://gallery.rcpp.org/articles/simulate-multivariate-normal/
*/
 arma::mat rmvnorm(arma::vec mu, arma::mat sigma)
{
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return  mu.t() + Y * arma::chol(sigma);
}

/*
 * This function generate a Dirichlet distribution using Gamma functions in R generator 
 */


arma::vec rdirichlet(const arma::vec alpha)
{
  arma::vec res(alpha.n_elem);
  for(arma::uword i=0;i<alpha.n_elem;i++)
  {
    res[i] = R::rgamma(alpha[i],1.0);
  }  
  return res/sum(res);
}


/*
 * This function compute the counts for a categorical vector given the vector (x) and the possible (unique) values it can assume (xu)
 */

NumericVector table(const NumericVector x, const NumericVector xu)
{
  std::map<double, int> Elt;
  
  Elt.clear();
  
  for(auto k=0; k < xu.size(); k++)
  {
    Elt[xu[k]] =0;
  }
  
  // Count each element
  for (auto i = 0; i != x.size(); ++i) 
  {
    Elt[ x[i] ] += 1;
  }
  
  return wrap(Elt);
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Function to update latent part of the model for all subject
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// note that we assume the input starting from 1,... NOT 0,... (as in R) !!!
// We combine all individual step in Algorithm 1 for computational speed


// [[Rcpp::export(".update_subject_latent")]]


List update_subject_latent(
  arma::field<arma::umat> data,              // observed data for t=1,..,T
  arma::field<arma::mat>  Kern,              // Kernels
  arma::field<arma::umat> ZI,                // profiles for each group
  arma::mat ki,                        // data from polya-gamma data augmentation
  arma::mat eta,                       // linear predictor for eta_i for each subject
  arma::mat psi,                       // profiles odds as in Polson & Scott
  arma::mat lambda,                    // profiles probabilities
  arma::mat omega,                     // polya-gamma latent variables
  arma::cube  ISigma_psi,              // inverse of the variance and covarianve matrix
  arma::vec subject_year_indicator,    // which year is subject i  
  arma::field<arma::uvec> g_j,               // groups indicator
  arma::vec i_t,                       // map from i=1,..,N to equivalent i_t
  arma::field<arma::vec> inv_j_t             // map from (t,j_t) to j=1,...,p
){
//////////////////////////////////////////
//////////////////////////////////////////
// define & compute useful quantities
//////////////////////////////////////////
//////////////////////////////////////////

//Poliagamma generator
#ifdef USE_R
GetRNGstate();
#endif

RNG PGr;
PolyaGamma PG;

// discrete generator
std::random_device rd;
std::mt19937 gen(rd());

int N = ki.n_rows;
int T = data.n_rows;

arma::uvec p_t(T);
arma::uvec n_t(T);
arma::umat pg_t(T,2);

for(int t =0;t<T;t++){
  p_t(t) = data(t).n_cols;
  n_t(t) = data(t).n_rows; 
  // just for 2 groups 
  pg_t(t,0) =  sum(g_j(t) == 1);
  pg_t(t,1) =  sum(g_j(t) == 2);
}

//int p = sum(p_t);
arma::vec pg;

// temporary quantities
int posg;
int posit;
int t;
arma::vec tmp_prob(2);
arma::mat sig_psi(2,2);
arma::vec mu_psi(2);


// 

for(auto i=0 ; i<N ; i++){
 // define year of subject i
 t  =  subject_year_indicator(i) -1;
 posit = i_t(i)-1;

// update the latent variable indicator 
 for(unsigned int j=0;j<p_t(t);j++){
  posg  = (g_j(t))(j) -1; 
  
  
  tmp_prob(0) = log(1-lambda(i,posg)) + log((Kern(inv_j_t(t)(j) -1))((data(t))(posit,j) -1,0));      
  tmp_prob(1) = log(lambda(i,posg))   + log((Kern(inv_j_t(t)(j) -1))((data(t))(posit,j) -1,1));   
  tmp_prob    =   log_vec_renorm(tmp_prob);
   
  // set probability for the generator 
  std::discrete_distribution<> zi_disc(tmp_prob.begin(),tmp_prob.end());
  // update latent class for the subject i
   (ZI(t))(posit,j) = zi_disc(gen);
    
 }


    for(auto g=0;g<2;g++)
    { 
    // update polya-gamma data
     arma::umat tmp_ZI = ((ZI(t)).row(posit));
     arma::umat tmp_ZI_sel  =  tmp_ZI.cols(find(g_j(t) == int(g +1)));
     ki(i,g) =  arma::as_scalar(sum(tmp_ZI_sel,1))- pg_t(t,g)/2.0;

    //polyagamma data augmentation 
     omega(i,g) = PG.draw(pg_t(t,g),psi(i,g),PGr);
    }


    //update the multivariate random effect psi
    sig_psi = arma::inv((ISigma_psi.slice(t)) + arma::diagmat(omega.row(i)));
    mu_psi = sig_psi*((ISigma_psi.slice(t))*eta.row(i).t() +  ki.row(i).t());
    psi.row(i) = rmvnorm(mu_psi, sig_psi);
    lambda.row(i) = plogis(psi.row(i).t()).t();

}

// return results to R
return List::create(
        _["eta"]    = eta, // this is not updated but is the input one
        _["omega"]  = omega,
        _["psi"]    = psi,
        _["lambda"] = lambda,
        _["ZI"]     = ZI,
        _["ki"]     = ki);
}




