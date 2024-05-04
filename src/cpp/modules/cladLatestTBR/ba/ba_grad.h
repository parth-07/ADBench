#include "clad/Differentiator/ArrayRef.h"

void computeReprojError_pullback(const double *cam, const double *X,
                                 const double *w, const double *feat,
                                 double *err, clad::array_ref<double> _d_cam,
                                 clad::array_ref<double> _d_X,
                                 clad::array_ref<double> _d_w,
                                 clad::array_ref<double> _d_feat,
                                 clad::array_ref<double> _d_err);

void computeZachWeightError_pullback(const double *const w, double *err,
                                     clad::array_ref<double> _d_w,
                                     clad::array_ref<double> _d_err);