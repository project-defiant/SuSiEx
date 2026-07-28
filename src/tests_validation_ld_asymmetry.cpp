#include <iostream>
#include <cstring>
#include <cstdlib>
#include "data.hpp"
#include "model.hpp"

int main()
{
    ::npop = 1;
    softpar par;
    par.n_gwas.push_back(10000);
    par.out_dir = ".";
    par.out_name = "test_susiex_invalid_ld_asym";
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

    dataset dat;
    try {
        dat.load_from_memory(par, nsnp, beta_in, pval_in, ind_in, ld_in, mkIdx_in);
    } catch(const std::exception &e) {
        std::cerr << "Validation asymmetry test: caught expected validation error: " << e.what() << std::endl;
        return 0; // success: expected failure
    }

    std::cerr << "Validation asymmetry test: expected failure due to asymmetric LD, but load_from_memory succeeded" << std::endl;
    return 1; // indicate failing test
}
