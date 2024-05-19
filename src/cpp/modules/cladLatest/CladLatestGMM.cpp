// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

#include "CladLatestGMM.h"

#include "gmm/gmm.h"
#include "gmm/gmm_grad.h"

// This function must be called before any other function.
void CladLatestGMM::prepare(GMMInput&& input)
{
    this->input = input;
    int Jcols = (this->input.k * (this->input.d + 1) * (this->input.d + 2)) / 2;
    result = { 0, std::vector<double>(Jcols) };
}



GMMOutput CladLatestGMM::output()
{
    return result;
}



void CladLatestGMM::calculate_objective(int times)
{
    for (int i = 0; i < times; i++)
    {
        gmm_objective(
            input.d,
            input.k,
            input.n,
            input.alphas.data(),
            input.means.data(),
            input.icf.data(),
            input.x.data(),
            input.wishart,
            &result.objective
        );
    }
}



void CladLatestGMM::calculate_jacobian(int times)
{
    double* alphas_gradient_part = result.gradient.data();
    double* means_gradient_part = result.gradient.data() + input.alphas.size();
    double* icf_gradient_part =
        result.gradient.data() +
        input.alphas.size() +
        input.means.size();
    std::vector<double> d_x(input.x.size(), 0);

    for (int i = 0; i < times; i++)
    {
        double tmp = 0.0;       // stores fictive result
                                // (Tapenade doesn't calculate an original function in reverse mode)

        double errb = 1.0;      // stores dY
                                // (equals to 1.0 for gradient calculation)
        
        int d_d = 0, d_k = 0, d_n = 0;
        Wishart d_wishart{0, 0};
        std::fill(result.gradient.begin(), result.gradient.end(), 0);
        std::fill(d_x.begin(), d_x.end(), 0);

        gmm_objective_pullback(
            input.d,
            input.k,
            input.n,
            input.alphas.data(),
            input.means.data(),
            input.icf.data(),
            input.x.data(),
            input.wishart,
            &tmp,
            &d_d,
            &d_k,
            &d_n,
            alphas_gradient_part,
            means_gradient_part,
            icf_gradient_part,
            d_x.data(),
            &d_wishart,
            &errb
        );
    }
}



extern "C" DLL_PUBLIC ITest<GMMInput, GMMOutput>* get_gmm_test()
{
    return new CladLatestGMM();
}