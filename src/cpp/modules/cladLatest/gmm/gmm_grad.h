#include "gmm.h"

void gmm_objective_pullback(int d, int k, int n, const double *const alphas, const double *const means, const double *const icf, const double *const x, Wishart wishart, double *err, int *_d_d, int *_d_k, int *_d_n, double *_d_alphas, double *_d_means, double *_d_icf, double *_d_x, Wishart *_d_wishart, double *_d_err);