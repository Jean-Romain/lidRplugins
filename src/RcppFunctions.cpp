#include <Rcpp.h>
#include <SpatialIndexes.h>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
IntegerVector C_count_in_disc(NumericVector X, NumericVector Y, NumericVector x, NumericVector y, double radius, int ncpu)
{
  unsigned int n = x.length();
  IntegerVector output(n);

  lidR::GridPartition tree(X,Y);

  #pragma omp parallel for num_threads(ncpu)
  for(unsigned int i = 0 ; i < n ; i++)
  {
    lidR::Circle disc(x[i], y[i], radius);
    std::vector<lidR::PointXYZ> pts;
    tree.lookup(disc, pts);

    #pragma omp critical
    {
      output[i] = pts.size();
    }
  }

  return output;
}
