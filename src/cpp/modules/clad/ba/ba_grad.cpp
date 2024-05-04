#include "ba.cpp"
#include "ba_grad.h"
#include "clad/Differentiator/Differentiator.h"

void sqsum_pullback(int n, const double *const x, double _d_y,
                    clad::array_ref<int> _d_n, clad::array_ref<double> _d_x) {
  double _d_res = 0;
  unsigned long _t0;
  int _d_i = 0;
  clad::tape<double> _t1 = {};
  clad::tape<int> _t2 = {};
  clad::tape<double> _t4 = {};
  clad::tape<int> _t5 = {};
  double res = 0;
  _t0 = 0;
  for (int i = 0; i < n; i++) {
    _t0++;
    res = res + clad::push(_t4, x[clad::push(_t2, i)]) *
                    clad::push(_t1, x[clad::push(_t5, i)]);
  }
  goto _label0;
_label0:
  _d_res += _d_y;
  for (; _t0; _t0--) {
    double _r_d0 = _d_res;
    _d_res += _r_d0;
    double _r0 = _r_d0 * clad::pop(_t1);
    int _t3 = clad::pop(_t2);
    _d_x[_t3] += _r0;
    double _r1 = clad::pop(_t4) * _r_d0;
    int _t6 = clad::pop(_t5);
    _d_x[_t6] += _r1;
    _d_res -= _r_d0;
  }
}
void cross_pullback(const double *const a, const double *const b, double *out,
                    clad::array_ref<double> _d_a, clad::array_ref<double> _d_b,
                    clad::array_ref<double> _d_out) {
  double _t0;
  double _t1;
  double _t2;
  double _t3;
  double _t4;
  double _t5;
  double _t6;
  double _t7;
  double _t8;
  double _t9;
  double _t10;
  double _t11;
  _t1 = a[1];
  _t0 = b[2];
  _t3 = a[2];
  _t2 = b[1];
  out[0] = _t1 * _t0 - _t3 * _t2;
  _t5 = a[2];
  _t4 = b[0];
  _t7 = a[0];
  _t6 = b[2];
  out[1] = _t5 * _t4 - _t7 * _t6;
  _t9 = a[0];
  _t8 = b[1];
  _t11 = a[1];
  _t10 = b[0];
  out[2] = _t9 * _t8 - _t11 * _t10;
  {
    double _r_d2 = _d_out[2];
    double _r8 = _r_d2 * _t8;
    _d_a[0] += _r8;
    double _r9 = _t9 * _r_d2;
    _d_b[1] += _r9;
    double _r10 = -_r_d2 * _t10;
    _d_a[1] += _r10;
    double _r11 = _t11 * -_r_d2;
    _d_b[0] += _r11;
    _d_out[2] -= _r_d2;
    _d_out[2];
  }
  {
    double _r_d1 = _d_out[1];
    double _r4 = _r_d1 * _t4;
    _d_a[2] += _r4;
    double _r5 = _t5 * _r_d1;
    _d_b[0] += _r5;
    double _r6 = -_r_d1 * _t6;
    _d_a[0] += _r6;
    double _r7 = _t7 * -_r_d1;
    _d_b[2] += _r7;
    _d_out[1] -= _r_d1;
    _d_out[1];
  }
  {
    double _r_d0 = _d_out[0];
    double _r0 = _r_d0 * _t0;
    _d_a[1] += _r0;
    double _r1 = _t1 * _r_d0;
    _d_b[2] += _r1;
    double _r2 = -_r_d0 * _t2;
    _d_a[2] += _r2;
    double _r3 = _t3 * -_r_d0;
    _d_b[1] += _r3;
    _d_out[0] -= _r_d0;
    _d_out[0];
  }
}
void rodrigues_rotate_point_pullback(const double rot[], const double pt[],
                                     double rotatedPt[],
                                     clad::array_ref<double> _d_rot,
                                     clad::array_ref<double> _d_pt,
                                     clad::array_ref<double> _d_rotatedPt) {
  const double *_t0;
  double _d_sqtheta = 0;
  bool _cond0;
  double _d_theta = 0, _d_costheta = 0, _d_sintheta = 0, _d_theta_inverse = 0,
         _d_tmp = 0;
  clad::array<double> _d_w(3UL), _d_w_cross_pt(3UL);
  double _t1;
  double _t2;
  double _t3;
  double _t4;
  unsigned long _t5;
  int _d_i = 0;
  clad::tape<int> _t6 = {};
  clad::tape<double> _t8 = {};
  clad::tape<int> _t9 = {};
  clad::tape<double> _t11 = {};
  clad::array<double> _t12(3UL);
  const double *_t13;
  clad::array<double> _t14(3UL);
  double _t15;
  double _t16;
  double _t17;
  double _t18;
  double _t19;
  double _t20;
  double _t21;
  double _t22;
  unsigned long _t23;
  int _d_i2 = 0;
  clad::tape<int> _t24 = {};
  clad::tape<double> _t26 = {};
  clad::tape<int> _t27 = {};
  clad::tape<double> _t29 = {};
  clad::tape<double> _t30 = {};
  clad::tape<int> _t31 = {};
  clad::tape<double> _t33 = {};
  clad::tape<double> _t34 = {};
  clad::tape<int> _t35 = {};
  clad::tape<double> _t37 = {};
  clad::array<double> _d_rot_cross_pt(3UL);
  const double *_t38;
  const double *_t39;
  clad::array<double> _t40(3UL);
  unsigned long _t41;
  int _d_i3 = 0;
  clad::tape<int> _t42 = {};
  clad::tape<int> _t44 = {};
  clad::tape<int> _t46 = {};
  _t0 = rot;
  double sqtheta = sqsum(3, rot);
  _cond0 = sqtheta != 0;
  if (_cond0) {
    double theta, costheta, sintheta, theta_inverse, tmp;
    double w[3], w_cross_pt[3];
    _t1 = sqtheta;
    theta = sqrt(_t1);
    _t2 = theta;
    costheta = cos(_t2);
    _t3 = theta;
    sintheta = sin(_t3);
    _t4 = theta;
    theta_inverse = 1. / _t4;
    _t5 = 0;
    for (int i = 0; i < 3; i++) {
      _t5++;
      w[clad::push(_t6, i)] = clad::push(_t11, rot[clad::push(_t9, i)]) *
                              clad::push(_t8, theta_inverse);
    }
    _t12 = w;
    _t13 = pt;
    _t14 = w_cross_pt;
    cross(w, pt, w_cross_pt);
    _t17 = w[0];
    _t16 = pt[0];
    _t19 = w[1];
    _t18 = pt[1];
    _t21 = w[2];
    _t20 = pt[2];
    _t22 = (_t17 * _t16 + _t19 * _t18 + _t21 * _t20);
    _t15 = (1. - costheta);
    tmp = _t22 * _t15;
    _t23 = 0;
    for (int i2 = 0; i2 < 3; i2++) {
      _t23++;
      rotatedPt[clad::push(_t24, i2)] =
          clad::push(_t29, pt[clad::push(_t27, i2)]) *
              clad::push(_t26, costheta) +
          clad::push(_t33, w_cross_pt[clad::push(_t31, i2)]) *
              clad::push(_t30, sintheta) +
          clad::push(_t37, w[clad::push(_t35, i2)]) * clad::push(_t34, tmp);
    }
  } else {
    double rot_cross_pt[3];
    _t38 = rot;
    _t39 = pt;
    _t40 = rot_cross_pt;
    cross(rot, pt, rot_cross_pt);
    _t41 = 0;
    for (int i3 = 0; i3 < 3; i3++) {
      _t41++;
      rotatedPt[clad::push(_t42, i3)] =
          pt[clad::push(_t44, i3)] + rot_cross_pt[clad::push(_t46, i3)];
    }
  }
  if (_cond0) {
    for (; _t23; _t23--) {
      int _t25 = clad::pop(_t24);
      double _r_d6 = _d_rotatedPt[_t25];
      double _r20 = _r_d6 * clad::pop(_t26);
      int _t28 = clad::pop(_t27);
      _d_pt[_t28] += _r20;
      double _r21 = clad::pop(_t29) * _r_d6;
      _d_costheta += _r21;
      double _r22 = _r_d6 * clad::pop(_t30);
      int _t32 = clad::pop(_t31);
      _d_w_cross_pt[_t32] += _r22;
      double _r23 = clad::pop(_t33) * _r_d6;
      _d_sintheta += _r23;
      double _r24 = _r_d6 * clad::pop(_t34);
      int _t36 = clad::pop(_t35);
      _d_w[_t36] += _r24;
      double _r25 = clad::pop(_t37) * _r_d6;
      _d_tmp += _r25;
      _d_rotatedPt[_t25] -= _r_d6;
      _d_rotatedPt[_t25];
    }
    {
      double _r_d5 = _d_tmp;
      double _r12 = _r_d5 * _t15;
      double _r13 = _r12 * _t16;
      _d_w[0] += _r13;
      double _r14 = _t17 * _r12;
      _d_pt[0] += _r14;
      double _r15 = _r12 * _t18;
      _d_w[1] += _r15;
      double _r16 = _t19 * _r12;
      _d_pt[1] += _r16;
      double _r17 = _r12 * _t20;
      _d_w[2] += _r17;
      double _r18 = _t21 * _r12;
      _d_pt[2] += _r18;
      double _r19 = _t22 * _r_d5;
      _d_costheta += -_r19;
      _d_tmp -= _r_d5;
    }
    {
      cross_pullback(_t12, _t13, _t14, _d_w, _d_pt, _d_w_cross_pt);
      clad::array<double> _r9(_d_w);
      clad::array<double> _r10(_d_pt);
      clad::array<double> _r11(_d_w_cross_pt);
    }
    for (; _t5; _t5--) {
      int _t7 = clad::pop(_t6);
      double _r_d4 = _d_w[_t7];
      double _r7 = _r_d4 * clad::pop(_t8);
      int _t10 = clad::pop(_t9);
      _d_rot[_t10] += _r7;
      double _r8 = clad::pop(_t11) * _r_d4;
      _d_theta_inverse += _r8;
      _d_w[_t7] -= _r_d4;
      _d_w[_t7];
    }
    {
      double _r_d3 = _d_theta_inverse;
      double _r5 = _r_d3 / _t4;
      double _r6 = _r_d3 * -1. / (_t4 * _t4);
      _d_theta += _r6;
      _d_theta_inverse -= _r_d3;
    }
    {
      double _r_d2 = _d_sintheta;
      double _r4 =
          _r_d2 *
          clad::custom_derivatives::sin_pushforward(_t3, 1.).pushforward;
      _d_theta += _r4;
      _d_sintheta -= _r_d2;
    }
    {
      double _r_d1 = _d_costheta;
      double _r3 =
          _r_d1 *
          clad::custom_derivatives::cos_pushforward(_t2, 1.).pushforward;
      _d_theta += _r3;
      _d_costheta -= _r_d1;
    }
    {
      double _r_d0 = _d_theta;
      double _r2 =
          _r_d0 *
          clad::custom_derivatives::sqrt_pushforward(_t1, 1.).pushforward;
      _d_sqtheta += _r2;
      _d_theta -= _r_d0;
    }
  } else {
    for (; _t41; _t41--) {
      int _t43 = clad::pop(_t42);
      double _r_d7 = _d_rotatedPt[_t43];
      int _t45 = clad::pop(_t44);
      _d_pt[_t45] += _r_d7;
      int _t47 = clad::pop(_t46);
      _d_rot_cross_pt[_t47] += _r_d7;
      _d_rotatedPt[_t43] -= _r_d7;
      _d_rotatedPt[_t43];
    }
    {
      cross_pullback(_t38, _t39, _t40, _d_rot, _d_pt, _d_rot_cross_pt);
      clad::array<double> _r26(_d_rot);
      clad::array<double> _r27(_d_pt);
      clad::array<double> _r28(_d_rot_cross_pt);
    }
  }
  {
    int _grad0 = 0;
    sqsum_pullback(3, _t0, _d_sqtheta, &_grad0, _d_rot);
    int _r0 = _grad0;
    clad::array<double> _r1(_d_rot);
  }
}
void radial_distort_pullback(const double rad_params[], double proj[],
                             clad::array_ref<double> _d_rad_params,
                             clad::array_ref<double> _d_proj) {
  double _d_rsq = 0, _d_L = 0;
  double *_t0;
  double _t1;
  double _t2;
  double _t3;
  double _t4;
  double _t5;
  double _t6;
  double _t7;
  double _t8;
  double _t9;
  double _t10;
  double rsq, L;
  _t0 = proj;
  rsq = sqsum(2, proj);
  _t2 = rad_params[0];
  _t1 = rsq;
  _t5 = rad_params[1];
  _t4 = rsq;
  _t6 = _t5 * _t4;
  _t3 = rsq;
  L = 1. + _t2 * _t1 + _t6 * _t3;
  _t8 = proj[0];
  _t7 = L;
  proj[0] = _t8 * _t7;
  _t10 = proj[1];
  _t9 = L;
  proj[1] = _t10 * _t9;
  {
    double _r_d3 = _d_proj[1];
    double _r10 = _r_d3 * _t9;
    _d_proj[1] += _r10;
    double _r11 = _t10 * _r_d3;
    _d_L += _r11;
    _d_proj[1] -= _r_d3;
    _d_proj[1];
  }
  {
    double _r_d2 = _d_proj[0];
    double _r8 = _r_d2 * _t7;
    _d_proj[0] += _r8;
    double _r9 = _t8 * _r_d2;
    _d_L += _r9;
    _d_proj[0] -= _r_d2;
    _d_proj[0];
  }
  {
    double _r_d1 = _d_L;
    double _r2 = _r_d1 * _t1;
    _d_rad_params[0] += _r2;
    double _r3 = _t2 * _r_d1;
    _d_rsq += _r3;
    double _r4 = _r_d1 * _t3;
    double _r5 = _r4 * _t4;
    _d_rad_params[1] += _r5;
    double _r6 = _t5 * _r4;
    _d_rsq += _r6;
    double _r7 = _t6 * _r_d1;
    _d_rsq += _r7;
    _d_L -= _r_d1;
  }
  {
    double _r_d0 = _d_rsq;
    int _grad0 = 0;
    sqsum_pullback(2, _t0, _r_d0, &_grad0, _d_proj);
    int _r0 = _grad0;
    clad::array<double> _r1(_d_proj);
    _d_rsq -= _r_d0;
  }
}
void project_pullback(const double cam[], const double X[], double proj[],
                      clad::array_ref<double> _d_cam,
                      clad::array_ref<double> _d_X,
                      clad::array_ref<double> _d_proj) {
  clad::array<double> _d_Xo(3UL), _d_Xcam(3UL);
  int _t0;
  int _t1;
  int _t2;
  const double *_t3;
  clad::array<double> _t4(3UL);
  clad::array<double> _t5(3UL);
  double _t6;
  double _t7;
  double _t8;
  double _t9;
  const double *_t10;
  double *_t11;
  double _t12;
  double _t13;
  double _t14;
  double _t15;
  double Xo[3], Xcam[3];
  _t0 = 3 + 0;
  Xo[0] = X[0] - cam[_t0];
  _t1 = 3 + 1;
  Xo[1] = X[1] - cam[_t1];
  _t2 = 3 + 2;
  Xo[2] = X[2] - cam[_t2];
  _t3 = cam;
  _t4 = Xo;
  _t5 = Xcam;
  rodrigues_rotate_point(cam, Xo, Xcam);
  _t7 = Xcam[0];
  _t6 = Xcam[2];
  proj[0] = _t7 / _t6;
  _t9 = Xcam[1];
  _t8 = Xcam[2];
  proj[1] = _t9 / _t8;
  _t10 = cam + 9;
  _t11 = proj;
  radial_distort(cam + 9, proj);
  _t13 = proj[0];
  _t12 = cam[6];
  proj[0] = _t13 * _t12 + cam[7];
  _t15 = proj[1];
  _t14 = cam[6];
  proj[1] = _t15 * _t14 + cam[8];
  {
    double _r_d6 = _d_proj[1];
    double _r11 = _r_d6 * _t14;
    _d_proj[1] += _r11;
    double _r12 = _t15 * _r_d6;
    _d_cam[6] += _r12;
    _d_cam[8] += _r_d6;
    _d_proj[1] -= _r_d6;
    _d_proj[1];
  }
  {
    double _r_d5 = _d_proj[0];
    double _r9 = _r_d5 * _t12;
    _d_proj[0] += _r9;
    double _r10 = _t13 * _r_d5;
    _d_cam[6] += _r10;
    _d_cam[7] += _r_d5;
    _d_proj[0] -= _r_d5;
    _d_proj[0];
  }
  {
    radial_distort_pullback(_t10, _t11, _d_cam.slice(9), _d_proj);
    clad::array<double> _r7(_d_cam.slice(9));
    clad::array<double> _r8(_d_proj);
  }
  {
    double _r_d4 = _d_proj[1];
    double _r5 = _r_d4 / _t8;
    _d_Xcam[1] += _r5;
    double _r6 = _r_d4 * -_t9 / (_t8 * _t8);
    _d_Xcam[2] += _r6;
    _d_proj[1] -= _r_d4;
    _d_proj[1];
  }
  {
    double _r_d3 = _d_proj[0];
    double _r3 = _r_d3 / _t6;
    _d_Xcam[0] += _r3;
    double _r4 = _r_d3 * -_t7 / (_t6 * _t6);
    _d_Xcam[2] += _r4;
    _d_proj[0] -= _r_d3;
    _d_proj[0];
  }
  {
    rodrigues_rotate_point_pullback(_t3, _t4, _t5, _d_cam, _d_Xo, _d_Xcam);
    clad::array<double> _r0(_d_cam);
    clad::array<double> _r1(_d_Xo);
    clad::array<double> _r2(_d_Xcam);
  }
  {
    double _r_d2 = _d_Xo[2];
    _d_X[2] += _r_d2;
    _d_cam[_t2] += -_r_d2;
    _d_Xo[2] -= _r_d2;
    _d_Xo[2];
  }
  {
    double _r_d1 = _d_Xo[1];
    _d_X[1] += _r_d1;
    _d_cam[_t1] += -_r_d1;
    _d_Xo[1] -= _r_d1;
    _d_Xo[1];
  }
  {
    double _r_d0 = _d_Xo[0];
    _d_X[0] += _r_d0;
    _d_cam[_t0] += -_r_d0;
    _d_Xo[0] -= _r_d0;
    _d_Xo[0];
  }
}
void computeReprojError_pullback(const double cam[], const double X[],
                                 const double w[], const double feat[],
                                 double err[], clad::array_ref<double> _d_cam,
                                 clad::array_ref<double> _d_X,
                                 clad::array_ref<double> _d_w,
                                 clad::array_ref<double> _d_feat,
                                 clad::array_ref<double> _d_err) {
  clad::array<double> _d_proj(2UL);
  const double *_t0;
  const double *_t1;
  clad::array<double> _t2(2UL);
  double _t3;
  double _t4;
  double _t5;
  double _t6;
  double proj[2];
  _t0 = cam;
  _t1 = X;
  _t2 = proj;
  project(cam, X, proj);
  _t4 = w[0];
  _t3 = (proj[0] - feat[0]);
  err[0] = _t4 * _t3;
  _t6 = w[0];
  _t5 = (proj[1] - feat[1]);
  err[1] = _t6 * _t5;
  {
    double _r_d1 = _d_err[1];
    double _r5 = _r_d1 * _t5;
    _d_w[0] += _r5;
    double _r6 = _t6 * _r_d1;
    _d_proj[1] += _r6;
    _d_feat[1] += -_r6;
    _d_err[1] -= _r_d1;
    _d_err[1];
  }
  {
    double _r_d0 = _d_err[0];
    double _r3 = _r_d0 * _t3;
    _d_w[0] += _r3;
    double _r4 = _t4 * _r_d0;
    _d_proj[0] += _r4;
    _d_feat[0] += -_r4;
    _d_err[0] -= _r_d0;
    _d_err[0];
  }
  {
    project_pullback(_t0, _t1, _t2, _d_cam, _d_X, _d_proj);
    clad::array<double> _r0(_d_cam);
    clad::array<double> _r1(_d_X);
    clad::array<double> _r2(_d_proj);
  }
}

void computeZachWeightError_pullback(const double *const w, double *err,
                                     clad::array_ref<double> _d_w,
                                     clad::array_ref<double> _d_err) {
  double _t0;
  double _t1;
  _t1 = w[0];
  _t0 = w[0];
  err[0] = 1 - _t1 * _t0;
  {
    double _r_d0 = _d_err[0];
    double _r0 = -_r_d0 * _t0;
    _d_w[0] += _r0;
    double _r1 = _t1 * -_r_d0;
    _d_w[0] += _r1;
    _d_err[0] -= _r_d0;
    _d_err[0];
  }
}