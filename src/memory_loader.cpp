#include "data.hpp"
#include "validation.hpp"
#include <cstring>
#include <iostream>
#include <stdexcept>

// Load dataset from in-memory arrays.
// Caller must set global `npop` before constructing dataset instance.
void dataset::load_from_memory(const softpar & par,
    int _nsnp,
    const double* beta_in,
    const double* pval_in,
    const char* ind_in,
    LDTYPE*** ld_in,
    const int* mkIdx_in)
{
    if(npop <= 0)
    {
        std::cerr << "Error: global npop must be set before calling load_from_memory" << std::endl;
        exit(1);
    }

    nsnp = _nsnp;

    // allocate arrays
    beta = new double[npop * nsnp];
    pval = new double[npop * nsnp];
    mkIdx = new int[nsnp];
    ind = new char[npop * nsnp];

    // copy
    memcpy(beta, beta_in, sizeof(double) * npop * nsnp);
    memcpy(pval, pval_in, sizeof(double) * npop * nsnp);
    memcpy(ind, ind_in, sizeof(char) * npop * nsnp);
    memcpy(mkIdx, mkIdx_in, sizeof(int) * nsnp);

    // tau_sq: set to max beta^2 per population (fallback to small positive if all zero)
    if(tau_sq)
        delete []tau_sq;
    tau_sq = new double[npop];
    for(int i = 0; i < npop; ++i)
    {
        double maxb = 0.0;
        for(int j = 0; j < nsnp; ++j)
        {
            double b = beta_in[i * nsnp + j];
            double b2 = b * b;
            if(b2 > maxb) maxb = b2;
        }
        if(maxb <= 0)
            tau_sq[i] = 1e-6;
        else
            tau_sq[i] = maxb;
    }

    // allocate and copy LD matrices
    ld = new LDTYPE**[npop];

    // validate using reusable validator
    std::string verr;
    if(!validate_ld_matrices(npop, nsnp, ld_in, verr))
    {
        // throw structured error instead of exiting; callers can catch
        throw std::runtime_error(verr);
    }

    for(int i = 0; i < npop; ++i)
    {
        ld[i] = new LDTYPE*[nsnp];
        for(int j = 0; j < nsnp; ++j)
        {
            ld[i][j] = new LDTYPE[nsnp];
            for(int k = 0; k < nsnp; ++k)
                ld[i][j][k] = ld_in[i][j][k];
        }
    }

    // populate mks and mkIdx to be compatible with existing write paths
    mks.clear();
    mks.resize(nsnp);
    for(int j = 0; j < nsnp; ++j)
    {
        snp &s = mks[j];
        s.coord.id = std::to_string(j);
        s.coord.pos = j + 1;
        s.coord.a1 = "A";
        s.coord.a2 = "T";
        // set per-population idxs and stats
        for(int i = 0; i < npop; ++i)
        {
            s.idxs[i * 3] = 5; // present in LD ref and GWAS
            s.idxs[i * 3 + 1] = j; // index
            s.idxs[i * 3 + 2] = j; // gwas idx
            s.stats[i * 4] = 0.5; // freq
            s.stats[i * 4 + 1] = beta_in[i * nsnp + j];
            s.stats[i * 4 + 2] = pval_in[i * nsnp + j];
            double lp = (pval_in[i * nsnp + j] > 0) ? -log10(pval_in[i * nsnp + j]) : 0.0;
            s.stats[i * 4 + 3] = lp;
        }
    }
    for(int j = 0; j < nsnp; ++j)
        mkIdx[j] = j;
}
