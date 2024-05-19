// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
#ifndef CLADLATEST_BA_H
#define CLADLATEST_BA_H
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "defs.h"
#include "matrix.h"

////////////////////////////////////////////////////////////
//////////////////// Declarations //////////////////////////
////////////////////////////////////////////////////////////

// cam: 11 camera in format [r1 r2 r3 C1 C2 C3 f u0 v0 k1 k2]
//            r1, r2, r3 are angle - axis rotation parameters(Rodrigues)
//            [C1 C2 C3]' is the camera center
//            f is the focal length in pixels
//            [u0 v0]' is the principal point
//            k1, k2 are radial distortion parameters
// X: 3 point
// feats: 2 feature (x,y coordinates)
// reproj_err: 2
// projection function:
// Xcam = R * (X - C)
// distorted = radial_distort(projective2euclidean(Xcam), radial_parameters)
// proj = distorted * f + principal_point
// err = sqsum(proj - measurement)
template <typename T>
void computeReprojError(const T *const cam, const T *const X, const T *const w,
                        const double *const feat, T *err);

// w: 1
// w_err: 1
template <typename T> void computeZachWeightError(const T *const w, T *err);

// n number of cameras
// m number of points
// p number of observations
// cams: 11*n cameras in format [r1 r2 r3 C1 C2 C3 f u0 v0 k1 k2]
//            r1, r2, r3 are angle - axis rotation parameters(Rodrigues)
//            [C1 C2 C3]' is the camera center
//            f is the focal length in pixels
//            [u0 v0]' is the principal point
//            k1, k2 are radial distortion parameters
// X: 3*m points
// obs: 2*p observations (pairs cameraIdx, pointIdx)
// feats: 2*p features (x,y coordinates corresponding to observations)
// reproj_err: 2*p errors of observations
// w_err: p weight "error" terms
// projection function:
// Xcam = R * (X - C)
// distorted = radial_distort(projective2euclidean(Xcam), radial_parameters)
// proj = distorted * f + principal_point
// err = sqsum(proj - measurement)
template <typename T>
void ba_objective(int n, int m, int p, const T *const cams, const T *const X,
                  const T *const w, const int *const obs,
                  const double *const feats, T *reproj_err, T *w_err);

// rot: 3 rotation parameters
// pt: 3 point to be rotated
// rotatedPt: 3 rotated point
// this is an efficient evaluation (part of
// the Ceres implementation)
// easy to understand calculation in matlab:
//  theta = sqrt(sum(w. ^ 2));
//  n = w / theta;
//  n_x = au_cross_matrix(n);
//  R = eye(3) + n_x*sin(theta) + n_x*n_x*(1 - cos(theta));
template <typename T>
void rodrigues_rotate_point(const T *const rot, const T *const pt,
                            T *rotatedPt);

////////////////////////////////////////////////////////////
//////////////////// Definitions ///////////////////////////
////////////////////////////////////////////////////////////

template <typename T> T sqsum(int n, const T *const x) {
  T res = 0;
  int i = 0;
  for (; i < n; i++)
    res = res + x[i] * x[i];
  return res;
}

template <typename T>
void rodrigues_rotate_point(const T * rot, const T * pt,
                            T *rotatedPt) {
  T sqtheta = sqsum(3, rot);
  int i0 = 0;
  int i1 = 0;
  int i2 = 0;
  T theta, costheta, sintheta, theta_inverse, tmp;
  T w[3], w_cross_pt[3];
  T rot_cross_pt[3];
  if (sqtheta != 0) {
    theta = sqrt(sqtheta);
    costheta = cos(theta);
    sintheta = sin(theta);
    theta_inverse = 1.0 / theta;

    for (; i0 < 3; i0++)
      w[i0] = rot[i0] * theta_inverse;

    cross(w, pt, w_cross_pt);

    tmp = (w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * (1. - costheta);

    for (; i1 < 3; i1++)
      rotatedPt[i1] = pt[i1] * costheta + w_cross_pt[i1] * sintheta + w[i1] * tmp;
  } else {
    cross(rot, pt, rot_cross_pt);

    for (; i2 < 3; i2++)
      rotatedPt[i2] = pt[i2] + rot_cross_pt[i2];
  }
}

template <typename T> void radial_distort(const T *const rad_params, T *proj) {
  T rsq, L;
  rsq = sqsum(2, proj);
  L = 1. + rad_params[0] * rsq + rad_params[1] * rsq * rsq;
  proj[0] = proj[0] * L;
  proj[1] = proj[1] * L;
}

template <typename T>
void project(const T * cam, const T * X, T *proj) {
  const T *const C = &cam[3];
  T Xo[3], Xcam[3];

  Xo[0] = X[0] - C[0];
  Xo[1] = X[1] - C[1];
  Xo[2] = X[2] - C[2];

  rodrigues_rotate_point(cam + 0, Xo, Xcam);

  proj[0] = Xcam[0] / Xcam[2];
  proj[1] = Xcam[1] / Xcam[2];

  radial_distort(cam + 9, proj);

  proj[0] = proj[0] * cam[6] + cam[7];
  proj[1] = proj[1] * cam[6] + cam[8];
}

template <typename T>
void computeReprojError(const T * cam, const T * X, const T * w,
                        const double * feat, T *err) {
  T proj[2];
  project(cam, X, proj);

  err[0] = (*w) * (proj[0] - feat[0]);
  err[1] = (*w) * (proj[1] - feat[1]);
}

template <typename T> void computeZachWeightError(const T *const w, T *err) {
  *err = 1 - (*w) * (*w);
}

template <typename T>
void ba_objective(int n, int m, int p, const T *const cams, const T *const X,
                  const T *const w, const int *const obs,
                  const double *const feats, T *reproj_err, T *w_err) {
  for (int i = 0; i < p; i++) {
    int camIdx = obs[i * 2 + 0];
    int ptIdx = obs[i * 2 + 1];
    computeReprojError(&cams[camIdx * BA_NCAMPARAMS], &X[ptIdx * 3], &w[i],
                       &feats[i * 2], &reproj_err[2 * i]);
  }

  for (int i = 0; i < p; i++) {
    computeZachWeightError(&w[i], &w_err[i]);
  }
}

template <typename T>
T computeReprojError_wrapper(const T * cam, const T * X,
                             const T * w, const double * feat,
                             T *err) {
  computeReprojError(cam, X, w, feat, err);
  return err[0];
}

template <typename T>
T computeZachWeightError_wrapper(const T * w, T *err) {
    computeZachWeightError(w, err);
    return err[0];
}
#endif