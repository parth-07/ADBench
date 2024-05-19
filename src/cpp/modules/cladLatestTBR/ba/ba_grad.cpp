#include "ba.h"
#include "clad/Differentiator/Differentiator.h"

void sqsum_pullback(int n, const double *const x, double _d_y,
                    clad::array_ref<int> _d_n, clad::array_ref<double> _d_x) {
  double _d_res = 0;
  int _d_i = 0;
  unsigned long _t0;
  double res = 0;
  int i = 0;
  _t0 = 0;
  for (; i < n; i++) {
    _t0++;
    res = res + x[i] * x[i];
  }
  goto _label0;
_label0:
  _d_res += _d_y;
  for (; _t0; _t0--) {
    i--;
    double _r_d0 = _d_res;
    _d_res += _r_d0;
    _d_x[i] += _r_d0 * x[i];
    _d_x[i] += x[i] * _r_d0;
    _d_res -= _r_d0;
  }
}
void cross_pullback(const double *const a, const double *const b, double *out,
                    clad::array_ref<double> _d_a, clad::array_ref<double> _d_b,
                    clad::array_ref<double> _d_out) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
  {
    double _r_d2 = _d_out[2];
    _d_a[0] += _r_d2 * b[1];
    _d_b[1] += a[0] * _r_d2;
    _d_a[1] += -_r_d2 * b[0];
    _d_b[0] += a[1] * -_r_d2;
    _d_out[2] -= _r_d2;
    _d_out[2];
  }
  {
    double _r_d1 = _d_out[1];
    _d_a[2] += _r_d1 * b[0];
    _d_b[0] += a[2] * _r_d1;
    _d_a[0] += -_r_d1 * b[2];
    _d_b[2] += a[0] * -_r_d1;
    _d_out[1] -= _r_d1;
    _d_out[1];
  }
  {
    double _r_d0 = _d_out[0];
    _d_a[1] += _r_d0 * b[2];
    _d_b[2] += a[1] * _r_d0;
    _d_a[2] += -_r_d0 * b[1];
    _d_b[1] += a[2] * -_r_d0;
    _d_out[0] -= _r_d0;
    _d_out[0];
  }
}
void rodrigues_rotate_point_pullback(const double *rot, const double *pt,
                                     double *rotatedPt,
                                     clad::array_ref<double> _d_rot,
                                     clad::array_ref<double> _d_pt,
                                     clad::array_ref<double> _d_rotatedPt) {
  const double *_t0;
  double _d_sqtheta = 0;
  int _d_i0 = 0;
  int _d_i1 = 0;
  int _d_i2 = 0;
  double _d_theta = 0, _d_costheta = 0, _d_sintheta = 0, _d_theta_inverse = 0,
         _d_tmp = 0;
  clad::array<double> _d_w(3UL), _d_w_cross_pt(3UL);
  clad::array<double> _d_rot_cross_pt(3UL);
  bool _cond0;
  unsigned long _t1;
  const double *_t2;
  unsigned long _t3;
  const double *_t4;
  const double *_t5;
  unsigned long _t6;
  _t0 = rot;
  double sqtheta = sqsum(3, rot);
  int i0 = 0;
  int i1 = 0;
  int i2 = 0;
  double theta, costheta, sintheta, theta_inverse, tmp;
  double w[3], w_cross_pt[3];
  double rot_cross_pt[3];
  _cond0 = sqtheta != 0;
  if (_cond0) {
    theta = sqrt(sqtheta);
    costheta = cos(theta);
    sintheta = sin(theta);
    theta_inverse = 1. / theta;
    _t1 = 0;
    for (; i0 < 3; i0++) {
      _t1++;
      w[i0] = rot[i0] * theta_inverse;
    }
    _t2 = pt;
    cross(w, pt, w_cross_pt);
    tmp = (w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * (1. - costheta);
    _t3 = 0;
    for (; i1 < 3; i1++) {
      _t3++;
      rotatedPt[i1] =
          pt[i1] * costheta + w_cross_pt[i1] * sintheta + w[i1] * tmp;
    }
  } else {
    _t4 = rot;
    _t5 = pt;
    cross(rot, pt, rot_cross_pt);
    _t6 = 0;
    for (; i2 < 3; i2++) {
      _t6++;
      rotatedPt[i2] = pt[i2] + rot_cross_pt[i2];
    }
  }
  if (_cond0) {
    for (; _t3; _t3--) {
      i1--;
      double _r_d6 = _d_rotatedPt[i1];
      _d_pt[i1] += _r_d6 * costheta;
      _d_costheta += pt[i1] * _r_d6;
      _d_w_cross_pt[i1] += _r_d6 * sintheta;
      _d_sintheta += w_cross_pt[i1] * _r_d6;
      _d_w[i1] += _r_d6 * tmp;
      _d_tmp += w[i1] * _r_d6;
      _d_rotatedPt[i1] -= _r_d6;
      _d_rotatedPt[i1];
    }
    {
      double _r_d5 = _d_tmp;
      _d_w[0] += _r_d5 * (1. - costheta) * pt[0];
      _d_pt[0] += w[0] * _r_d5 * (1. - costheta);
      _d_w[1] += _r_d5 * (1. - costheta) * pt[1];
      _d_pt[1] += w[1] * _r_d5 * (1. - costheta);
      _d_w[2] += _r_d5 * (1. - costheta) * pt[2];
      _d_pt[2] += w[2] * _r_d5 * (1. - costheta);
      _d_costheta += -(w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * _r_d5;
      _d_tmp -= _r_d5;
    }
    {
      pt = _t2;
      cross_pullback(w, _t2, w_cross_pt, _d_w, _d_pt, _d_w_cross_pt);
    }
    for (; _t1; _t1--) {
      i0--;
      double _r_d4 = _d_w[i0];
      _d_rot[i0] += _r_d4 * theta_inverse;
      _d_theta_inverse += rot[i0] * _r_d4;
      _d_w[i0] -= _r_d4;
      _d_w[i0];
    }
    {
      double _r_d3 = _d_theta_inverse;
      double _r4 = _r_d3 * -1. / (theta * theta);
      _d_theta += _r4;
      _d_theta_inverse -= _r_d3;
    }
    {
      double _r_d2 = _d_sintheta;
      double _r3 =
          _r_d2 *
          clad::custom_derivatives::sin_pushforward(theta, 1.).pushforward;
      _d_theta += _r3;
      _d_sintheta -= _r_d2;
    }
    {
      double _r_d1 = _d_costheta;
      double _r2 =
          _r_d1 *
          clad::custom_derivatives::cos_pushforward(theta, 1.).pushforward;
      _d_theta += _r2;
      _d_costheta -= _r_d1;
    }
    {
      double _r_d0 = _d_theta;
      double _r1 =
          _r_d0 *
          clad::custom_derivatives::sqrt_pushforward(sqtheta, 1.).pushforward;
      _d_sqtheta += _r1;
      _d_theta -= _r_d0;
    }
  } else {
    for (; _t6; _t6--) {
      i2--;
      double _r_d7 = _d_rotatedPt[i2];
      _d_pt[i2] += _r_d7;
      _d_rot_cross_pt[i2] += _r_d7;
      _d_rotatedPt[i2] -= _r_d7;
      _d_rotatedPt[i2];
    }
    {
      rot = _t4;
      pt = _t5;
      cross_pullback(_t4, _t5, rot_cross_pt, _d_rot, _d_pt, _d_rot_cross_pt);
    }
  }
  {
    rot = _t0;
    int _grad0 = 0;
    sqsum_pullback(3, _t0, _d_sqtheta, &_grad0, _d_rot);
    int _r0 = _grad0;
  }
}
void radial_distort_pullback(const double *const rad_params, double *proj,
                             clad::array_ref<double> _d_rad_params,
                             clad::array_ref<double> _d_proj) {
  double _d_rsq = 0, _d_L = 0;
  double *_t0;
  double _t1;
  double _t2;
  double rsq, L;
  _t0 = proj;
  rsq = sqsum(2, proj);
  L = 1. + rad_params[0] * rsq + rad_params[1] * rsq * rsq;
  _t1 = proj[0];
  proj[0] = proj[0] * L;
  _t2 = proj[1];
  proj[1] = proj[1] * L;
  {
    proj[1] = _t2;
    double _r_d3 = _d_proj[1];
    _d_proj[1] += _r_d3 * L;
    _d_L += proj[1] * _r_d3;
    _d_proj[1] -= _r_d3;
    _d_proj[1];
  }
  {
    proj[0] = _t1;
    double _r_d2 = _d_proj[0];
    _d_proj[0] += _r_d2 * L;
    _d_L += proj[0] * _r_d2;
    _d_proj[0] -= _r_d2;
    _d_proj[0];
  }
  {
    double _r_d1 = _d_L;
    _d_rad_params[0] += _r_d1 * rsq;
    _d_rsq += rad_params[0] * _r_d1;
    _d_rad_params[1] += _r_d1 * rsq * rsq;
    _d_rsq += rad_params[1] * _r_d1 * rsq;
    _d_rsq += rad_params[1] * rsq * _r_d1;
    _d_L -= _r_d1;
  }
  {
    double _r_d0 = _d_rsq;
    proj = _t0;
    int _grad0 = 0;
    sqsum_pullback(2, _t0, _r_d0, &_grad0, _d_proj);
    int _r0 = _grad0;
    _d_rsq -= _r_d0;
  }
}
void project_pullback(const double *cam, const double *X, double *proj,
                      clad::array_ref<double> _d_cam,
                      clad::array_ref<double> _d_X,
                      clad::array_ref<double> _d_proj) {
  double *_d_C = 0;
  clad::array<double> _d_Xo(3UL), _d_Xcam(3UL);
  const double *_t0;
  const double *_t1;
  double *_t2;
  double _t3;
  double _t4;
  _d_C = &_d_cam[3];
  const double *const C = &cam[3];
  double Xo[3], Xcam[3];
  Xo[0] = X[0] - C[0];
  Xo[1] = X[1] - C[1];
  Xo[2] = X[2] - C[2];
  _t0 = cam + 0;
  rodrigues_rotate_point(cam + 0, Xo, Xcam);
  proj[0] = Xcam[0] / Xcam[2];
  proj[1] = Xcam[1] / Xcam[2];
  _t1 = cam + 9;
  _t2 = proj;
  // CFIX-start
  double tproj0 = proj[0];
  double tproj1 = proj[1];
  // CFIX-end
  radial_distort(cam + 9, proj);
  _t3 = proj[0];
  proj[0] = proj[0] * cam[6] + cam[7];
  _t4 = proj[1];
  proj[1] = proj[1] * cam[6] + cam[8];
  {
    proj[1] = _t4;
    double _r_d6 = _d_proj[1];
    _d_proj[1] += _r_d6 * cam[6];
    _d_cam[6] += proj[1] * _r_d6;
    _d_cam[8] += _r_d6;
    _d_proj[1] -= _r_d6;
    _d_proj[1];
  }
  {
    proj[0] = _t3;
    double _r_d5 = _d_proj[0];
    _d_proj[0] += _r_d5 * cam[6];
    _d_cam[6] += proj[0] * _r_d5;
    _d_cam[7] += _r_d5;
    _d_proj[0] -= _r_d5;
    _d_proj[0];
  }
  {
    proj = _t2;
    proj[0] = tproj0;
    proj[1] = tproj1;
    radial_distort_pullback(_t1, _t2, _d_cam.ptr_ref() + 9, _d_proj);
  }
  {
    double _r_d4 = _d_proj[1];
    _d_Xcam[1] += _r_d4 / Xcam[2];
    double _r1 = _r_d4 * -Xcam[1] / (Xcam[2] * Xcam[2]);
    _d_Xcam[2] += _r1;
    _d_proj[1] -= _r_d4;
    _d_proj[1];
  }
  {
    double _r_d3 = _d_proj[0];
    _d_Xcam[0] += _r_d3 / Xcam[2];
    double _r0 = _r_d3 * -Xcam[0] / (Xcam[2] * Xcam[2]);
    _d_Xcam[2] += _r0;
    _d_proj[0] -= _r_d3;
    _d_proj[0];
  }
  rodrigues_rotate_point_pullback(_t0, Xo, Xcam, _d_cam.ptr_ref() + 0, _d_Xo,
                                  _d_Xcam);
  {
    double _r_d2 = _d_Xo[2];
    _d_X[2] += _r_d2;
    _d_C[2] += -_r_d2;
    _d_Xo[2] -= _r_d2;
    _d_Xo[2];
  }
  {
    double _r_d1 = _d_Xo[1];
    _d_X[1] += _r_d1;
    _d_C[1] += -_r_d1;
    _d_Xo[1] -= _r_d1;
    _d_Xo[1];
  }
  {
    double _r_d0 = _d_Xo[0];
    _d_X[0] += _r_d0;
    _d_C[0] += -_r_d0;
    _d_Xo[0] -= _r_d0;
    _d_Xo[0];
  }
}
void computeReprojError_pullback(const double *cam, const double *X,
                                 const double *w, const double *feat,
                                 double *err, clad::array_ref<double> _d_cam,
                                 clad::array_ref<double> _d_X,
                                 clad::array_ref<double> _d_w,
                                 clad::array_ref<double> _d_feat,
                                 clad::array_ref<double> _d_err) {
  clad::array<double> _d_proj(2UL);
  const double *_t0;
  const double *_t1;
  double proj[2];
  _t0 = cam;
  _t1 = X;
  project(cam, X, proj);
  err[0] = *w * (proj[0] - feat[0]);
  err[1] = *w * (proj[1] - feat[1]);
  {
    double _r_d1 = _d_err[1];
    *_d_w += _r_d1 * (proj[1] - feat[1]);
    _d_proj[1] += *w * _r_d1;
    _d_feat[1] += -*w * _r_d1;
    _d_err[1] -= _r_d1;
    _d_err[1];
  }
  {
    double _r_d0 = _d_err[0];
    *_d_w += _r_d0 * (proj[0] - feat[0]);
    _d_proj[0] += *w * _r_d0;
    _d_feat[0] += -*w * _r_d0;
    _d_err[0] -= _r_d0;
    _d_err[0];
  }
  {
    cam = _t0;
    X = _t1;
    project_pullback(_t0, _t1, proj, _d_cam, _d_X, _d_proj);
  }
}
void computeReprojError_wrapper_grad(
    const double *cam, const double *X, const double *w, const double *feat,
    double *err, clad::array_ref<double> _d_cam, clad::array_ref<double> _d_X,
    clad::array_ref<double> _d_w, clad::array_ref<double> _d_feat,
    clad::array_ref<double> _d_err) {
  const double *_t0;
  const double *_t1;
  const double *_t2;
  const double *_t3;
  double *_t4;
  _t0 = cam;
  _t1 = X;
  _t2 = w;
  _t3 = feat;
  _t4 = err;
  computeReprojError(cam, X, w, feat, err);
  goto _label0;
_label0:
  _d_err[0] += 1;
  {
    cam = _t0;
    X = _t1;
    w = _t2;
    feat = _t3;
    err = _t4;
    computeReprojError_pullback(_t0, _t1, _t2, _t3, _t4, _d_cam, _d_X, _d_w,
                                _d_feat, _d_err);
  }
}
void computeZachWeightError_pullback(const double *const w, double *err,
                                     clad::array_ref<double> _d_w,
                                     clad::array_ref<double> _d_err) {
  double _t0;
  _t0 = *err;
  *err = 1 - *w * (*w);
  {
    *err = _t0;
    double _r_d0 = *_d_err;
    *_d_w += -_r_d0 * (*w);
    *_d_w += *w * -_r_d0;
    *_d_err -= _r_d0;
    *_d_err;
  }
}
void computeZachWeightError_wrapper_grad(const double *w, double *err,
                                         clad::array_ref<double> _d_w,
                                         clad::array_ref<double> _d_err) {
  const double *_t0;
  double *_t1;
  _t0 = w;
  _t1 = err;
  computeZachWeightError(w, err);
  goto _label0;
_label0:
  _d_err[0] += 1;
  {
    w = _t0;
    err = _t1;
    computeZachWeightError_pullback(_t0, _t1, _d_w, _d_err);
  }
}