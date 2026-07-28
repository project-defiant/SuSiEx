#include <iostream>
#include <cstring>
#include <cstdlib>
#include "data.hpp"

int main()
{
    // Set global npop before constructing dataset
    ::npop = 1;

    softpar par;

    int nsnp = 2;
    double beta_in[1 * 2] = {0.1, 0.2};
    double pval_in[1 * 2] = {0.01, 0.02};
    char ind_in[1 * 2] = {1, 1};
    int mkIdx_in[2] = {0, 1};

    // allocate ld_in
    LDTYPE*** ld_in = new LDTYPE**[1];
    ld_in[0] = new LDTYPE*[nsnp];
    for(int j = 0; j < nsnp; ++j)
    {
        ld_in[0][j] = new LDTYPE[nsnp];
        for(int k = 0; k < nsnp; ++k)
            ld_in[0][j][k] = (j == k) ? 1.0f : 0.1f;
    }

    dataset dat;
    dat.load_from_memory(par, nsnp, beta_in, pval_in, ind_in, ld_in, mkIdx_in);

    if(dat.nsnp != nsnp)
    {
        std::cerr << "nsnp mismatch" << std::endl;
        return 2;
    }
    if(dat.beta[0] != beta_in[0] || dat.beta[1] != beta_in[1])
    {
        std::cerr << "beta mismatch" << std::endl;
        return 3;
    }
    if(dat.pval[0] != pval_in[0] || dat.pval[1] != pval_in[1])
    {
        std::cerr << "pval mismatch" << std::endl;
        return 4;
    }
    if(dat.ind[0] != ind_in[0] || dat.ind[1] != ind_in[1])
    {
        std::cerr << "ind mismatch" << std::endl;
        return 5;
    }
    if(dat.mkIdx[0] != 0 || dat.mkIdx[1] != 1)
    {
        std::cerr << "mkIdx mismatch" << std::endl;
        return 6;
    }
    // check ld
    if(dat.ld[0][0][0] != 1.0f || dat.ld[0][0][1] != 0.1f || dat.ld[0][1][0] != 0.1f || dat.ld[0][1][1] != 1.0f)
    {
        std::cerr << "ld mismatch" << std::endl;
        return 7;
    }

    std::cout << "memory_loader_test: OK" << std::endl;

    // cleanup ld_in
    for(int j = 0; j < nsnp; ++j)
        delete [] ld_in[0][j];
    delete [] ld_in[0];
    delete [] ld_in;

    return 0;
}
