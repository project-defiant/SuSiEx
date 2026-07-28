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
    par.out_name = "test_susiex_invalid_ld";
    par.chr = "1";
    par.start = 1;
    par.end = 1000000;
    par.n_sig = 2;
    par.max_iter = 5;
    par.tol = 1e-6;
    par.nthreads = 1;

    int nsnp = 4;
    double beta_in[1 * 4] = {0.1, 0.0, 0.0, 0.0};
    double pval_in[1 * 4] = {1e-6, 1.0, 1.0, 1.0};
    char ind_in[1 * 4] = {1,1,1,1};
    int mkIdx_in[4] = {0,1,2,3};

    // allocate ld_in with invalid diagonal (not 1.0)
    LDTYPE*** ld_in = new LDTYPE**[1];
    ld_in[0] = new LDTYPE*[nsnp];
    for(int j = 0; j < nsnp; ++j)
    {
        ld_in[0][j] = new LDTYPE[nsnp];
        for(int k = 0; k < nsnp; ++k)
        {
            if(j == k)
                ld_in[0][j][k] = 0.9f; // invalid diagonal
            else
                ld_in[0][j][k] = 0.05f;
        }
    }

    dataset dat;
    try {
        dat.load_from_memory(par, nsnp, beta_in, pval_in, ind_in, ld_in, mkIdx_in);
    } catch(const std::exception &e) {
        std::cerr << "Validation diag test: caught expected validation error: " << e.what() << std::endl;
        return 0; // success: expected failure
    }

    std::cerr << "Validation diag test: expected failure due to invalid LD diagonal, but load_from_memory succeeded" << std::endl;
    return 1; // indicate failing test
}
