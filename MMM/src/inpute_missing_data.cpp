#include "RcppArmadillo.h"
// Rcpp::plugins(cpp11)
// Rcpp::depends(RcppArmadillo)

using namespace Rcpp;
using namespace arma;


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Function to inpute missing data 
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// note that we assume the input to have 1,... NOT 0,... (as in R) !!!

// [[Rcpp::export(".inpute_missing_data")]]



arma::field<arma::umat>   inpute_missing_data(
                           const arma::field<arma::umat> data,              // observed data for t=1,..,T
                           arma::field<arma::mat>  Kern,              // Kernels
                           arma::field<arma::mat>  missing_indices,    // row and columns of the missing indices for t=1,,T (!! FROM 1 not 0!!!)
                           arma::mat lambda,                    // profiles probabilities
                           arma::vec subject_year_indicator,    // which year is subject i  
                           arma::field<arma::uvec> g_j,               // groups indicator
                           arma::vec i_t,                       // map from i=1,..,N to equivalent i_t
                           arma::field<arma::vec> inv_j_t,             // map from (t,j_t) to j=1,...,p
                           arma::field<arma::vec> inv_i_t             // map from (t,j_t) to j=1,...,p
)
{
 
// discrete generator
std::random_device rd;
std::mt19937 gen(rd());
 
// temporary quantities
int T = data.n_rows;
int posj;
int posi;
int posg;
unsigned int j;
arma::vec tmp_prob(2);

arma::field<umat> results = data;

for( auto t= 0 ; t<T ; t++)
{ 
 for(unsigned int i= 0; i< (missing_indices(t)).n_rows ; i++ )
 {
   
   // set indicies  
   posi = (missing_indices(t))(i,0) -1;
   posj = (missing_indices(t))(i,1) -1 ;
   posg = (g_j(t)(posj)) -1;
   j    = inv_j_t(t)(posj) -1;

   //////////////////////////////////////////
   //////////////////////////////////////////
   // this is tought for binary variables
   //////////////////////////////////////////
   //////////////////////////////////////////

  tmp_prob(0) = exp( 
          log(1-lambda(inv_i_t(t)(i) -1 ,posg)) + log((Kern(j))(0,0)) ) + 
                 exp(log(lambda(inv_i_t(t)(i) -1,posg)) + log((Kern(j)(0,1))));      


  tmp_prob(1) = exp(log(1-lambda(inv_i_t(t)(i) -1,posg)) + log((Kern(j))(1,0))) + 
                 exp(log(lambda(inv_i_t(t)(i) -1,posg))  + log((Kern(j))(1,1)));      

  // not necessay
  //tmp_prob    =   log_vec_renorm(tmp_prob);


  // set probability for the generator 
  std::discrete_distribution<> new_data_prob(tmp_prob.begin(),tmp_prob.end());
  // update the data
  (results(t))(posi,posj) = new_data_prob(gen)  + 1;
 }
}
return results;
}
