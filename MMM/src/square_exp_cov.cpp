#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::plugins(cpp11)]]

// [[Rcpp::export(".square_exp_cov")]]
NumericMatrix square_exp_cov(NumericMatrix Coords,NumericVector length_scale,double eta2,double nugget){
int N = Coords.nrow();
//int d = Coords.ncol();
NumericMatrix results(N,N);
for(int i=0; i < N ; i++){
 for(int j=0; j<=i ;  j++){
  results(i,j) =  eta2 * exp( -0.5 *  sum(length_scale * pow(Coords(i,_) - Coords(j,_),2)));
 if(i==j){
     results(i,i) += nugget;
}
else{
    results(j,i) = results(i,j);
}
 }
} 
return(results);
}

