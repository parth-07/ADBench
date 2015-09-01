/*        Generated by TAPENADE     (INRIA, Tropics team)
Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
*/
#include "ba.h"
#include "adStack.h"

#define NB_DIRS_REPROJ_BV 2

/*
Differentiation of sqsum in reverse (adjoint) mode:
gradient     of useful results: *x sqsum
with respect to varying inputs: *x
Plus diff mem management of: x:in
*/
void sqsum_bv(int n,
  double *x,
  double(*xb)[NBDirsMaxReproj_BV],
  double sqsumb[NBDirsMaxReproj_BV],
  int nbdirs);

/*
Differentiation of cross in reverse (adjoint) mode:
gradient     of useful results: *out *a *b
with respect to varying inputs: *a *b
Plus diff mem management of: out:in a:in b:in
*/
void cross_bv(
  double *a,
  double(*ab)[NBDirsMaxReproj_BV],
  double *b,
  double(*bb)[NBDirsMaxReproj_BV],
  double *out,
  double(*outb)[NBDirsMaxReproj_BV],
  int nbdirs);

// rot 3 rotation parameters
// pt 3 point to be rotated
// rotatedPt 3 rotated point
// this is an efficient evaluation (part of
// the Ceres implementation)
// easy to understand calculation in matlab:
//	theta = sqrt(sum(w. ^ 2));
//	n = w / theta;
//	n_x = au_cross_matrix(n);
//	R = eye(3) + n_x*sin(theta) + n_x*n_x*(1 - cos(theta));
/*
Differentiation of rodrigues_rotate_point in reverse (adjoint) mode:
gradient     of useful results: *rot *rotatedPt
with respect to varying inputs: *rot *pt
Plus diff mem management of: rot:in rotatedPt:in pt:in
*/
void rodrigues_rotate_point_bv(
  double *rot,
  double(*rotb)[NBDirsMaxReproj_BV],
  double *pt,
  double(*ptb)[NBDirsMaxReproj_BV],
  double *rotatedPt,
  double(*rotatedPtb)[NBDirsMaxReproj_BV],
  int nbdirs);

// rad_params 2 radial distortion parameters
// proj 2 projection to be distorted
/*
Differentiation of radial_distort in reverse (adjoint) mode:
gradient     of useful results: *rad_params *proj
with respect to varying inputs: *rad_params *proj
Plus diff mem management of: rad_params:in proj:in
*/
void radial_distort_bv(
  double *rad_params,
  double(*rad_paramsb)[NBDirsMaxReproj_BV],
  double *proj,
  double(*projb)[NBDirsMaxReproj_BV],
  int nbdirs);

// cam 11 cameras in format [r1 r2 r3 C1 C2 C3 f u0 v0 k1 k2]
//            r1, r2, r3 are angle - axis rotation parameters(Rodrigues)
//			  [C1 C2 C3]' is the camera center
//            f is the focal length in pixels
//			  [u0 v0]' is the principal point
//            k1, k2 are radial distortion parameters
// X 3 point
// proj 2 projection
// projection: 
// Xcam = R * (X - C)
// distorted = radial_distort(projective2euclidean(Xcam), radial_parameters)
// proj = distorted * f + principal_point
// err = sqsum(proj - measurement)
/*
Differentiation of project in reverse (adjoint) mode:
gradient     of useful results: *cam *X *proj
with respect to varying inputs: *cam *X
Plus diff mem management of: cam:in X:in proj:in-out
*/
void project_bv(
  double *cam,
  double(*camb)[NBDirsMaxReproj_BV],
  double *X,
  double(*Xb)[NBDirsMaxReproj_BV],
  double *proj,
  double(*projb)[NBDirsMaxReproj_BV],
  int nbdirs);

/*
Differentiation of computeReprojError in reverse (adjoint) mode:
gradient     of useful results: *err
with respect to varying inputs: *err *w *cam *X
RW status of diff variables: *err:in-out *w:out *cam:out *X:out
Plus diff mem management of: err:in w:in cam:in X:in
*/
void computeReprojError_bv(
  double *cam,
  double(*camb)[NBDirsMaxReproj_BV],
  double *X,
  double(*Xb)[NBDirsMaxReproj_BV],
  double *w,
  double(*wb)[NBDirsMaxReproj_BV],
  double feat_x, double feat_y,
  double *err,
  double(*errb)[NBDirsMaxReproj_BV],
  int nbdirs);

void computeZachWeightError_b(double *w, double *wb, double *err, double *errb);