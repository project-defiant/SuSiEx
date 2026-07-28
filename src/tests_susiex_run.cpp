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
    par.out_name = "test_susiex";
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
    char ind_in[1 * 5] = {1,1,1,1,1};
    int mkIdx_in[5] = {0,1,2,3,4};

    // allocate ld_in
    LDTYPE*** ld_in = new LDTYPE**[1];
    ld_in[0] = new LDTYPE*[nsnp];
    for(int j = 0; j < nsnp; ++j)
    {
        ld_in[0][j] = new LDTYPE[nsnp];
        for(int k = 0; k < nsnp; ++k)
            ld_in[0][j][k] = (j == k) ? 1.0f : 0.05f;
    }

    dataset dat;
    dat.load_from_memory(par, nsnp, beta_in, pval_in, ind_in, ld_in, mkIdx_in);

    susiex model(::npop, dat.nsnp, dat, par);

    model.susie_sst_xethn();

    std::cout << "susiex run: OK" << std::endl;

    // cleanup
    for(int j = 0; j < nsnp; ++j)
        delete [] ld_in[0][j];
    delete [] ld_in[0];
    delete [] ld_in;

    return 0;
}
