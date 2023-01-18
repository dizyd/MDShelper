#include <Rcpp.h>
using namespace Rcpp;


//' Calculate pairwise distance between points
//'
//' @param x numeric matrix with dimensions in the columns and stimuli in the rows
//' @param r Minkowski metric (1 = city block, 2 = euclidean)
//'

// [[Rcpp::export]]
NumericMatrix pair_dist_cpp(NumericMatrix x, double r) {

  int nrow = x.nrow();
  int ndim = x.ncol();
  NumericMatrix out(nrow,nrow);


  for(int i = 0; i < nrow; i ++){

    for( int j = 0; j < nrow; j++){

      double t_p = 0;

      for(int k = 0; k < ndim; k++){

       t_p += pow(x(i,k)-x(j,k),r);

      }

      out(i,j) = pow(t_p,1/r);

    }

  }


  return out;
}


