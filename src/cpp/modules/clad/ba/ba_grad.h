#include "clad/Differentiator/ArrayRef.h"
void sqsum_pullback(int n, const double *const x, double _d_y,
                    clad::array_ref<int> _d_n, clad::array_ref<double> _d_x);

void rodrigues_rotate_point_pullback(const double rot[], const double pt[],
                                     double rotatedPt[],
                                     clad::array_ref<double> _d_rot,
                                     clad::array_ref<double> _d_pt,
                                     clad::array_ref<double> _d_rotatedPt);

void radial_distort_pullback(const double rad_params[], double proj[],
                             clad::array_ref<double> _d_rad_params,
                             clad::array_ref<double> _d_proj);

void project_pullback(const double cam[], const double X[], double proj[],
                      clad::array_ref<double> _d_cam,
                      clad::array_ref<double> _d_X,
                      clad::array_ref<double> _d_proj);

void computeReprojError_pullback(const double cam[], const double X[],
                                 const double w[], const double feat[],
                                 double err[], clad::array_ref<double> _d_cam,
                                 clad::array_ref<double> _d_X,
                                 clad::array_ref<double> _d_w,
                                 clad::array_ref<double> _d_feat,
                                 clad::array_ref<double> _d_err);

void computeZachWeightError_pullback(const double *const w, double *err,
                                     clad::array_ref<double> _d_w,
                                     clad::array_ref<double> _d_err);