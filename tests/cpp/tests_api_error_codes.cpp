#include <iostream>
#include "susiex/api.hpp"

int main()
{
    ::npop = 1;
    softpar par;
    par.n_gwas.push_back(10000);
    par.out_dir = ".";
    par.out_name = "test_api_error_codes";
    par.chr = "1";
    par.start = 1;
    par.end = 1000000;
    par.n_sig = 2;
    par.max_iter = 5;
    par.tol = 1e-6;
    par.nthreads = 1;

    int nsnp = 3;
    double beta_in[1 * 3] = {0.1, 0.0, 0.0};
    double pval_in[1 * 3] = {1e-6, 1.0, 1.0};
    char ind_in[1 * 3] = {1,1,1};
    int mkIdx_in[3] = {0,1,2};

    // allocate ld_in with asymmetry: ld[0][1] != ld[1][0]
    LDTYPE*** ld_in = new LDTYPE**[1];
    ld_in[0] = new LDTYPE*[nsnp];
    for(int j = 0; j < nsnp; ++j)
    {
        ld_in[0][j] = new LDTYPE[nsnp];
        for(int k = 0; k < nsnp; ++k)
            ld_in[0][j][k] = (j == k) ? 1.0f : 0.05f;
    }
    // introduce asymmetry
    ld_in[0][0][1] = 0.1f;
    ld_in[0][1][0] = 0.2f;

    ms_result out;
    int rc = susiex_multisusie_fit(1, nsnp, beta_in, pval_in, ind_in, ld_in, mkIdx_in, &par, &out);

    // expect validation error
    if(rc != SSEX_VALIDATION_ERROR)
    {
        std::cerr << "API error-code test: expected SSEX_VALIDATION_ERROR (" << SSEX_VALIDATION_ERROR << ") but got " << rc << std::endl;
        return 1;
    }

    std::cout << "API error-code test: OK (caught validation error via return code)" << std::endl;
    return 0;
}
