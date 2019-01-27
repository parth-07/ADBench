#ifndef GMM_B_LOADED
#define GMM_B_LOADED
/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.10 (r5717) - 30 Jul 2015 16:03
*/
#include "../defs.h"
#include "gmm.h"

// d: dim
// k: number of gaussians
// n: number of points
// alphas: k logs of mixture weights (unnormalized), so
//			weights = exp(log_alphas) / sum(exp(log_alphas))
// means: d*k component means
// icf: (d*(d+1)/2)*k inverse covariance factors 
//					every icf entry stores firstly log of diagonal and then 
//          columnwise other entris
//          To generate icf in MATLAB given covariance C :
//              L = inv(chol(C, 'lower'));
//              inv_cov_factor = [log(diag(L)); L(au_tril_indices(d, -1))]
// wishart: wishart distribution parameters
// x: d*n points
// err: 1 output
void gmm_objective(int d, int k, int n, double *alphas, double *means, double 
    *icf, double *x, Wishart wishart, double *err);
void gmm_objective_b(int d, int k, int n, double *alphas, double *alphasb, 
    double *means, double *meansb, double *icf, double *icfb, double *x, 
    Wishart wishart, double *err, double *errb);


void gmm_objective_split_inner_b(int d, int k,
  double *alphas,
  double *alphasb,
  double *means,
  double *meansb,
  double *icf,
  double *icfb,
  double *x,
  double *err,
  double *errb);
void gmm_objective_split_other_b(int d, int k, int n,
  double *alphas,
  double *alphasb,
  double *icf,
  double *icfb,
  Wishart wishart,
  double *err,
  double *errb);

#endif