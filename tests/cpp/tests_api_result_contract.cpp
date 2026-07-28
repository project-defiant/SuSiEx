#include <cmath>
#include <iostream>

#include "susiex/api.hpp"

int main()
{
    ::npop = 1;
    softpar par;
    par.n_gwas.push_back(10000);
    par.out_dir = ".";
    par.out_name = "test_api_result_contract";
    par.chr = "1";
    par.start = 1;
    par.end = 1000000;
    par.n_sig = 2;
    par.max_iter = 10;
    par.tol = 1e-6;
    par.nthreads = 1;

    int nsnp = 5;
    double beta_in[1 * 5] = {0.2, 0.0, 0.0, 0.0, 0.0};
    double pval_in[1 * 5] = {1e-6, 1.0, 1.0, 1.0, 1.0};
    char ind_in[1 * 5] = {1, 1, 1, 1, 1};
    int mkIdx_in[5] = {0, 1, 2, 3, 4};

    LDTYPE*** ld_in = new LDTYPE**[1];
    ld_in[0] = new LDTYPE*[nsnp];
    for(int j = 0; j < nsnp; ++j)
    {
        ld_in[0][j] = new LDTYPE[nsnp];
        for(int k = 0; k < nsnp; ++k)
            ld_in[0][j][k] = (j == k) ? 1.0f : 0.05f;
    }

    ms_result out{};
    int rc = susiex_multisusie_fit(
        1,
        nsnp,
        beta_in,
        pval_in,
        ind_in,
        ld_in,
        mkIdx_in,
        &par,
        &out
    );

    for(int j = 0; j < nsnp; ++j)
        delete[] ld_in[0][j];
    delete[] ld_in[0];
    delete[] ld_in;

    if(rc != SSEX_OK)
    {
        std::cerr << "expected SSEX_OK but got " << rc << std::endl;
        return 1;
    }
    if(out.converged != 1)
    {
        std::cerr << "expected converged result" << std::endl;
        return 1;
    }
    if(!out.cs_purity || !out.cs_min_p || !out.cs_min_p_by_population || !out.logbf_by_population || !out.population_causal_prob || !out.cs_filtered)
    {
        std::cerr << "expected diagnostics arrays to be populated" << std::endl;
        ms_result_free(&out);
        return 1;
    }
    if(std::isnan(out.population_causal_prob[0]) || out.population_causal_prob[0] < 0.0 || out.population_causal_prob[0] > 1.0)
    {
        std::cerr << "expected valid population causal probability" << std::endl;
        ms_result_free(&out);
        return 1;
    }
    if(out.cs_purity[0] <= 0.0)
    {
        std::cerr << "expected positive component purity" << std::endl;
        ms_result_free(&out);
        return 1;
    }

    ms_result_free(&out);
    std::cout << "API result-contract test: OK" << std::endl;
    return 0;
}
