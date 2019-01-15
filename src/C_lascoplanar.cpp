#include <RcppArmadillo.h>
#include "Progress.h"
#include "QuadTree.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector C_lascoplanar(S4 las, int k, double th1, double th2)
{
  DataFrame data = as<Rcpp::DataFrame>(las.slot("data"));

  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];

  unsigned int n = X.length();

  LogicalVector output(n);

  QuadTree qtree(las);
  arma::mat A(k,3);
  arma::mat coeff;
  arma::mat score;
  arma::vec latent;

  Progress pb(n, "Eigenvalues computation: ");

  bool colinear_mode = (th2 == 0) ? true : false;

  for (unsigned int i = 0 ; i < n ; i++)
  {
    pb.check_abort();
    pb.increment();

    PointXYZ p(X[i], Y[i], Z[i]);

    std::vector<PointXYZ> pts;
    qtree.knn(p, k, pts);

    for (unsigned int j = 0 ; j < pts.size() ; j++)
    {
      A(j,0) = pts[j].x;
      A(j,1) = pts[j].y;
      A(j,2) = pts[j].z;
    }

    arma::princomp(coeff, score, latent, A);

    if (!colinear_mode && latent[1] > th1*latent[2] && th2*latent[1] > latent[0])
      output[i] = true;
    else if (colinear_mode && th1*latent[2] < latent[0] && th1*latent[1] < latent[0])
      output[i] = true;
  }

  return(output);
}