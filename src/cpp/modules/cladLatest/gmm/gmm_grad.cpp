#include "gmm.h"
// CLADFIX-Add-Start
#include "clad/Differentiator/Differentiator.h"
// CLADFIX-Add-End

void gmm_objective_pullback(int d, int k, int n, const double *const alphas, const double *const means, const double *const icf, const double *const x, Wishart wishart, double *err, int *_d_d, int *_d_k, int *_d_n, double *_d_alphas, double *_d_means, double *_d_icf, double *_d_x, Wishart *_d_wishart, double *_d_err);
void gmm_objective_wrapper_grad(int d, int k, int n, const double *const alphas, const double *const means, const double *const icf, const double *const x, Wishart wishart, double *err, int *_d_d, int *_d_k, int *_d_n, double *_d_alphas, double *_d_means, double *_d_icf, double *_d_x, Wishart *_d_wishart, double *_d_err) {
    gmm_objective(d, k, n, alphas, means, icf, x, wishart, err);
    goto _label0;
  _label0:
    _d_err[0] += 1;
    {
        int _r0 = 0;
        int _r1 = 0;
        int _r2 = 0;
        Wishart _r3 = {};
        gmm_objective_pullback(d, k, n, alphas, means, icf, x, wishart, err, &_r0, &_r1, &_r2, _d_alphas, _d_means, _d_icf, _d_x, &_r3, _d_err);
        *_d_d += _r0;
        *_d_k += _r1;
        *_d_n += _r2;
    }
}
void preprocess_qs_pullback(int d, int k, const double *const icf, double *sum_qs, double *Qdiags, int *_d_d, int *_d_k, double *_d_icf, double *_d_sum_qs, double *_d_Qdiags);
void subtract_pullback(int d, const double *const x, const double *const y, double *out, int *_d_d, double *_d_x, double *_d_y, double *_d_out);
void Qtimesx_pullback(int d, const double *const Qdiag, const double *const ltri, const double *const x, double *out, int *_d_d, double *_d_Qdiag, double *_d_ltri, double *_d_x, double *_d_out);
void sqnorm_pullback(int n, const double *const x, double _d_y, int *_d_n, double *_d_x);
void logsumexp_pullback(int n, const double *const x, double _d_y, int *_d_n, double *_d_x);
void log_wishart_prior_pullback(int p, int k, Wishart wishart, const double *const sum_qs, const double *const Qdiags, const double *const icf, double _d_y, int *_d_p, int *_d_k, Wishart *_d_wishart, double *_d_sum_qs, double *_d_Qdiags, double *_d_icf);
void gmm_objective_pullback(int d, int k, int n, const double *const alphas, const double *const means, const double *const icf, const double *const x, Wishart wishart, double *err, int *_d_d, int *_d_k, int *_d_n, double *_d_alphas, double *_d_means, double *_d_icf, double *_d_x, Wishart *_d_wishart, double *_d_err) {
    double _t0;
    double _d_CONSTANT = 0;
    int _d_icf_sz = 0;
    double _d_slse = 0;
    unsigned long _t1;
    int _d_ix = 0;
    int ix = 0;
    clad::tape<unsigned long> _t2 = {};
    clad::tape<int> _t3 = {};
    int _d_ik = 0;
    int ik = 0;
    clad::tape<double> _t4 = {};
    clad::tape<double> _t5 = {};
    clad::tape<double> _t6 = {};
    double _d_lse_alphas = 0;
    double _t7;
    double _t8;
    _t0 = log(2 * 3.1415926535900001);
    const double CONSTANT = -n * d * 0.5 * _t0;
    int icf_sz = d * (d + 1) / 2;
    double _d_Qdiags[d * k];
    clad::zero_init(_d_Qdiags, d * k);
    double Qdiags[d * k];
    double _d_sum_qs[k];
    clad::zero_init(_d_sum_qs, k);
    double sum_qs[k];
    double _d_xcentered[d];
    clad::zero_init(_d_xcentered, d);
    double xcentered[d];
    double _d_Qxcentered[d];
    clad::zero_init(_d_Qxcentered, d);
    double Qxcentered[d];
    double _d_main_term[k];
    clad::zero_init(_d_main_term, k);
    double main_term[k];
    preprocess_qs(d, k, icf, &sum_qs[0], &Qdiags[0]);
    double slse = 0.;
    _t1 = 0;
    for (ix = 0; ix < n; ix++) {
        _t1++;
        clad::push(_t2, 0UL);
        for (clad::push(_t3, ik) , ik = 0; ik < k; ik++) {
            clad::back(_t2)++;
            subtract(d, &x[ix * d], &means[ik * d], &xcentered[0]);
            Qtimesx(d, &Qdiags[ik * d], &icf[ik * icf_sz + d], &xcentered[0], &Qxcentered[0]);
            clad::push(_t4, main_term[ik]);
            main_term[ik] = alphas[ik] + sum_qs[ik] - 0.5 * clad::push(_t5, sqnorm(d, &Qxcentered[0]));
        }
        clad::push(_t6, slse);
        slse = slse + logsumexp(k, &main_term[0]);
    }
    double lse_alphas = logsumexp(k, alphas);
    _t7 = *err;
    *err = CONSTANT + slse - n * lse_alphas;
    _t8 = *err;
    *err = *err + log_wishart_prior(d, k, wishart, &sum_qs[0], &Qdiags[0], icf);
    {
        *err = _t8;
        double _r_d3 = *_d_err;
        *_d_err -= _r_d3;
        *_d_err += _r_d3;
        int _r7 = 0;
        int _r8 = 0;
        Wishart _r9 = {};
        log_wishart_prior_pullback(d, k, wishart, &sum_qs[0], &Qdiags[0], icf, _r_d3, &_r7, &_r8, &_r9, &_d_sum_qs[0], &_d_Qdiags[0], _d_icf);
        *_d_d += _r7;
        *_d_k += _r8;
    }
    {
        *err = _t7;
        double _r_d2 = *_d_err;
        *_d_err -= _r_d2;
        _d_CONSTANT += _r_d2;
        _d_slse += _r_d2;
        *_d_n += -_r_d2 * lse_alphas;
        _d_lse_alphas += n * -_r_d2;
    }
    {
        int _r6 = 0;
        logsumexp_pullback(k, alphas, _d_lse_alphas, &_r6, _d_alphas);
        *_d_k += _r6;
    }
    for (; _t1; _t1--) {
        ix--;
        {
            slse = clad::pop(_t6);
            double _r_d1 = _d_slse;
            _d_slse -= _r_d1;
            _d_slse += _r_d1;
            int _r5 = 0;
            logsumexp_pullback(k, &main_term[0], _r_d1, &_r5, &_d_main_term[0]);
            *_d_k += _r5;
        }
        {
            for (; clad::back(_t2); clad::back(_t2)--) {
                ik--;
                {
                    main_term[ik] = clad::pop(_t4);
                    double _r_d0 = _d_main_term[ik];
                    _d_main_term[ik] -= _r_d0;
                    _d_alphas[ik] += _r_d0;
                    _d_sum_qs[ik] += _r_d0;
                    int _r4 = 0;
                    sqnorm_pullback(d, &Qxcentered[0], 0.5 * -_r_d0, &_r4, &_d_Qxcentered[0]);
                    *_d_d += _r4;
                }
                {
                    int _r3 = 0;
                    Qtimesx_pullback(d, &Qdiags[ik * d], &icf[ik * icf_sz + d], &xcentered[0], &Qxcentered[0], &_r3, &_d_Qdiags[ik * d], &_d_icf[ik * icf_sz + d], &_d_xcentered[0], &_d_Qxcentered[0]);
                    *_d_d += _r3;
                }
                {
                    int _r2 = 0;
                    subtract_pullback(d, &x[ix * d], &means[ik * d], &xcentered[0], &_r2, &_d_x[ix * d], &_d_means[ik * d], &_d_xcentered[0]);
                    *_d_d += _r2;
                }
            }
            {
                _d_ik = 0;
                ik = clad::pop(_t3);
            }
            clad::pop(_t2);
        }
    }
    {
        int _r0 = 0;
        int _r1 = 0;
        preprocess_qs_pullback(d, k, icf, &sum_qs[0], &Qdiags[0], &_r0, &_r1, _d_icf, &_d_sum_qs[0], &_d_Qdiags[0]);
        *_d_d += _r0;
        *_d_k += _r1;
    }
    {
        *_d_d += _d_icf_sz / 2 * (d + 1);
        *_d_d += d * _d_icf_sz / 2;
    }
    {
        *_d_n += -_d_CONSTANT * _t0 * 0.5 * d;
        *_d_d += -n * _d_CONSTANT * _t0 * 0.5;
    }
}
void preprocess_qs_pullback(int d, int k, const double *const icf, double *sum_qs, double *Qdiags, int *_d_d, int *_d_k, double *_d_icf, double *_d_sum_qs, double *_d_Qdiags) {
    int _d_icf_sz = 0;
    unsigned long _t0;
    int _d_ik = 0;
    int ik = 0;
    clad::tape<double> _t1 = {};
    clad::tape<unsigned long> _t2 = {};
    clad::tape<int> _t3 = {};
    int _d_id = 0;
    int id = 0;
    clad::tape<double> _t4 = {};
    double _d_q = 0;
    double q = 0;
    clad::tape<double> _t5 = {};
    clad::tape<double> _t6 = {};
    int icf_sz = d * (d + 1) / 2;
    _t0 = 0;
    for (ik = 0; ik < k; ik++) {
        _t0++;
        clad::push(_t1, sum_qs[ik]);
        sum_qs[ik] = 0.;
        clad::push(_t2, 0UL);
        for (clad::push(_t3, id) , id = 0; id < d; id++) {
            clad::back(_t2)++;
            clad::push(_t4, q) , q = icf[ik * icf_sz + id];
            clad::push(_t5, sum_qs[ik]);
            sum_qs[ik] = sum_qs[ik] + q;
            clad::push(_t6, Qdiags[ik * d + id]);
            Qdiags[ik * d + id] = exp(q);
        }
    }
    for (; _t0; _t0--) {
        ik--;
        {
            for (; clad::back(_t2); clad::back(_t2)--) {
                id--;
                {
                    Qdiags[ik * d + id] = clad::pop(_t6);
                    double _r_d2 = _d_Qdiags[ik * d + id];
                    _d_Qdiags[ik * d + id] -= _r_d2;
                    double _r0 = 0;
                    _r0 += _r_d2 * clad::custom_derivatives::exp_pushforward(q, 1.).pushforward;
                    _d_q += _r0;
                }
                {
                    sum_qs[ik] = clad::pop(_t5);
                    double _r_d1 = _d_sum_qs[ik];
                    _d_sum_qs[ik] -= _r_d1;
                    _d_sum_qs[ik] += _r_d1;
                    _d_q += _r_d1;
                }
                {
                    _d_icf[ik * icf_sz + id] += _d_q;
                    _d_q = 0;
                    q = clad::pop(_t4);
                }
            }
            {
                _d_id = 0;
                id = clad::pop(_t3);
            }
            clad::pop(_t2);
        }
        {
            sum_qs[ik] = clad::pop(_t1);
            double _r_d0 = _d_sum_qs[ik];
            _d_sum_qs[ik] -= _r_d0;
        }
    }
    {
        *_d_d += _d_icf_sz / 2 * (d + 1);
        *_d_d += d * _d_icf_sz / 2;
    }
}
void subtract_pullback(int d, const double *const x, const double *const y, double *out, int *_d_d, double *_d_x, double *_d_y, double *_d_out) {
    unsigned long _t0;
    int _d_id = 0;
    int id = 0;
    clad::tape<double> _t1 = {};
    _t0 = 0;
    for (id = 0; id < d; id++) {
        _t0++;
        clad::push(_t1, out[id]);
        out[id] = x[id] - y[id];
    }
    for (; _t0; _t0--) {
        id--;
        {
            out[id] = clad::pop(_t1);
            double _r_d0 = _d_out[id];
            _d_out[id] -= _r_d0;
            _d_x[id] += _r_d0;
            _d_y[id] += -_r_d0;
        }
    }
}
void Qtimesx_pullback(int d, const double *const Qdiag, const double *const ltri, const double *const x, double *out, int *_d_d, double *_d_Qdiag, double *_d_ltri, double *_d_x, double *_d_out) {
    unsigned long _t0;
    int _d_id = 0;
    int id = 0;
    clad::tape<double> _t1 = {};
    int _d_Lparamsidx = 0;
    unsigned long _t2;
    int _d_i = 0;
    int i = 0;
    clad::tape<unsigned long> _t3 = {};
    clad::tape<int> _t4 = {};
    int _d_j = 0;
    int j = 0;
    clad::tape<double> _t5 = {};
    _t0 = 0;
    for (id = 0; id < d; id++) {
        _t0++;
        clad::push(_t1, out[id]);
        out[id] = Qdiag[id] * x[id];
    }
    int Lparamsidx = 0;
    _t2 = 0;
    for (i = 0; i < d; i++) {
        _t2++;
        clad::push(_t3, 0UL);
        for (clad::push(_t4, j) , j = i + 1; j < d; j++) {
            clad::back(_t3)++;
            clad::push(_t5, out[j]);
            out[j] = out[j] + ltri[Lparamsidx] * x[i];
            Lparamsidx++;
        }
    }
    for (; _t2; _t2--) {
        i--;
        {
            for (; clad::back(_t3); clad::back(_t3)--) {
                j--;
                Lparamsidx--;
                {
                    out[j] = clad::pop(_t5);
                    double _r_d1 = _d_out[j];
                    _d_out[j] -= _r_d1;
                    _d_out[j] += _r_d1;
                    _d_ltri[Lparamsidx] += _r_d1 * x[i];
                    _d_x[i] += ltri[Lparamsidx] * _r_d1;
                }
            }
            {
                _d_i += _d_j;
                _d_j = 0;
                j = clad::pop(_t4);
            }
            clad::pop(_t3);
        }
    }
    for (; _t0; _t0--) {
        id--;
        out[id] = clad::pop(_t1);
        double _r_d0 = _d_out[id];
        _d_out[id] -= _r_d0;
        _d_Qdiag[id] += _r_d0 * x[id];
        _d_x[id] += Qdiag[id] * _r_d0;
    }
}
void sqnorm_pullback(int n, const double *const x, double _d_y, int *_d_n, double *_d_x) {
    double _d_res = 0;
    unsigned long _t0;
    int _d_i = 0;
    int i = 0;
    clad::tape<double> _t1 = {};
    double res = x[0] * x[0];
    _t0 = 0;
    for (i = 1; i < n; i++) {
        _t0++;
        clad::push(_t1, res);
        res = res + x[i] * x[i];
    }
    goto _label0;
  _label0:
    _d_res += _d_y;
    for (; _t0; _t0--) {
        i--;
        res = clad::pop(_t1);
        double _r_d0 = _d_res;
        _d_res -= _r_d0;
        _d_res += _r_d0;
        _d_x[i] += _r_d0 * x[i];
        _d_x[i] += x[i] * _r_d0;
    }
    {
        _d_x[0] += _d_res * x[0];
        _d_x[0] += x[0] * _d_res;
    }
}
void arr_max_pullback(int n, const double *const x, double _d_y, int *_d_n, double *_d_x);
void logsumexp_pullback(int n, const double *const x, double _d_y, int *_d_n, double *_d_x) {
    double _d_mx = 0;
    double _d_semx = 0;
    unsigned long _t0;
    int _d_i = 0;
    int i = 0;
    clad::tape<double> _t1 = {};
    double mx = arr_max(n, x);
    double semx = 0.;
    _t0 = 0;
    for (i = 0; i < n; i++) {
        _t0++;
        clad::push(_t1, semx);
        semx = semx + exp(x[i] - mx);
    }
    goto _label0;
  _label0:
    {
        double _r2 = 0;
        _r2 += _d_y * clad::custom_derivatives::log_pushforward(semx, 1.).pushforward;
        _d_semx += _r2;
        _d_mx += _d_y;
    }
    for (; _t0; _t0--) {
        i--;
        {
            semx = clad::pop(_t1);
            double _r_d0 = _d_semx;
            _d_semx -= _r_d0;
            _d_semx += _r_d0;
            double _r1 = 0;
            _r1 += _r_d0 * clad::custom_derivatives::exp_pushforward(x[i] - mx, 1.).pushforward;
            _d_x[i] += _r1;
            _d_mx += -_r1;
        }
    }
    {
        int _r0 = 0;
        arr_max_pullback(n, x, _d_mx, &_r0, _d_x);
        *_d_n += _r0;
    }
}
inline void log_gamma_distrib_pullback(double a, double p, double _d_y, double *_d_a, double *_d_p);
void log_wishart_prior_pullback(int p, int k, Wishart wishart, const double *const sum_qs, const double *const Qdiags, const double *const icf, double _d_y, int *_d_p, int *_d_k, Wishart *_d_wishart, double *_d_sum_qs, double *_d_Qdiags, double *_d_icf) {
    int _d_n = 0;
    int _d_icf_sz = 0;
    double _t0;
    double _t1;
    double _d_C = 0;
    double _d_out = 0;
    unsigned long _t2;
    int _d_ik = 0;
    int ik = 0;
    clad::tape<double> _t3 = {};
    double _d_frobenius = 0;
    double frobenius = 0;
    clad::tape<double> _t4 = {};
    int n = p + wishart.m + 1;
    int icf_sz = p * (p + 1) / 2;
    _t1 = log(2);
    _t0 = (log(wishart.gamma) - 0.5 * _t1);
    double C = n * p * _t0 - log_gamma_distrib(0.5 * n, p);
    double out = 0;
    _t2 = 0;
    for (ik = 0; ik < k; ik++) {
        _t2++;
        clad::push(_t3, frobenius) , frobenius = sqnorm(p, &Qdiags[ik * p]) + sqnorm(icf_sz - p, &icf[ik * icf_sz + p]);
        clad::push(_t4, out);
        out = out + 0.5 * wishart.gamma * wishart.gamma * (frobenius) - wishart.m * sum_qs[ik];
    }
    goto _label0;
  _label0:
    {
        _d_out += _d_y;
        *_d_k += -_d_y * C;
        _d_C += k * -_d_y;
    }
    for (; _t2; _t2--) {
        ik--;
        {
            out = clad::pop(_t4);
            double _r_d0 = _d_out;
            _d_out -= _r_d0;
            _d_out += _r_d0;
            (*_d_wishart).gamma += 0.5 * _r_d0 * (frobenius) * wishart.gamma;
            (*_d_wishart).gamma += 0.5 * wishart.gamma * _r_d0 * (frobenius);
            _d_frobenius += 0.5 * wishart.gamma * wishart.gamma * _r_d0;
            (*_d_wishart).m += -_r_d0 * sum_qs[ik];
            _d_sum_qs[ik] += wishart.m * -_r_d0;
        }
        {
            int _r3 = 0;
            sqnorm_pullback(p, &Qdiags[ik * p], _d_frobenius, &_r3, &_d_Qdiags[ik * p]);
            *_d_p += _r3;
            int _r4 = 0;
            sqnorm_pullback(icf_sz - p, &icf[ik * icf_sz + p], _d_frobenius, &_r4, &_d_icf[ik * icf_sz + p]);
            _d_icf_sz += _r4;
            *_d_p += -_r4;
            _d_frobenius = 0;
            frobenius = clad::pop(_t3);
        }
    }
    {
        _d_n += _d_C * _t0 * p;
        *_d_p += n * _d_C * _t0;
        double _r0 = 0;
        _r0 += n * p * _d_C * clad::custom_derivatives::log_pushforward(wishart.gamma, 1.).pushforward;
        (*_d_wishart).gamma += _r0;
        double _r1 = 0;
        double _r2 = 0;
        log_gamma_distrib_pullback(0.5 * n, p, -_d_C, &_r1, &_r2);
        _d_n += 0.5 * _r1;
        *_d_p += _r2;
    }
    {
        *_d_p += _d_icf_sz / 2 * (p + 1);
        *_d_p += p * _d_icf_sz / 2;
    }
    {
        *_d_p += _d_n;
        (*_d_wishart).m += _d_n;
    }
}
void arr_max_pullback(int n, const double *const x, double _d_y, int *_d_n, double *_d_x) {
    double _d_m = 0;
    unsigned long _t0;
    int _d_i = 0;
    int i = 0;
    clad::tape<bool> _t2 = {};
    clad::tape<double> _t3 = {};
    double m = x[0];
    _t0 = 0;
    for (i = 1; i < n; i++) {
        _t0++;
        bool _t1 = m < x[i];
        {
            if (_t1) {
                clad::push(_t3, m);
                m = x[i];
            }
            clad::push(_t2, _t1);
        }
    }
    goto _label0;
  _label0:
    _d_m += _d_y;
    for (; _t0; _t0--) {
        i--;
        if (clad::pop(_t2)) {
            m = clad::pop(_t3);
            double _r_d0 = _d_m;
            _d_m -= _r_d0;
            _d_x[i] += _r_d0;
        }
    }
    _d_x[0] += _d_m;
}
inline void log_gamma_distrib_pullback(double a, double p, double _d_y, double *_d_a, double *_d_p) {
    double _t0;
    double _d_out = 0;
    unsigned long _t1;
    int _d_j = 0;
    int j = 0;
    clad::tape<double> _t2 = {};
    _t0 = log(3.1415926535900001);
    double out = 0.25 * p * (p - 1) * _t0;
    _t1 = 0;
    for (j = 1; j <= p; j++) {
        _t1++;
        clad::push(_t2, out);
        out = out + lgamma(a + 0.5 * (1 - j));
    }
    goto _label0;
  _label0:
    _d_out += _d_y;
    for (; _t1; _t1--) {
        j--;
        {
            out = clad::pop(_t2);
            double _r_d0 = _d_out;
            _d_out -= _r_d0;
            _d_out += _r_d0;
            double _r0 = 0;
            _r0 += _r_d0 * numerical_diff::forward_central_difference(lgamma, a + 0.5 * (1 - j), 0, 0, a + 0.5 * (1 - j));
            *_d_a += _r0;
            _d_j += -0.5 * _r0;
        }
    }
    {
        *_d_p += 0.25 * _d_out * _t0 * (p - 1);
        *_d_p += 0.25 * p * _d_out * _t0;
    }
}