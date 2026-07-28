#pragma once

#include "susiex/data.hpp"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    SSEX_OK = 0,
    SSEX_VALIDATION_ERROR = 1,
    SSEX_INVALID_ARGS = 2,
    SSEX_INTERNAL_ERROR = -1
} ssex_code;

// Minimal result struct for now; expand later
typedef struct {
    int converged;
    int L;      // number of components (nsig)
    int ncs;    // number of credible sets
    int K;      // npop
    int P;      // nsnp

    double *alpha; // length L * P (alpha[l*P + p])
    double *pip;   // length P
    double *mu;    // length K * L * P, layout: l major then k then p as mu[(l*K + k)*P + p]
    double *lbf;   // length L

    int *cs_counts; // length L
    int *cs_indices; // flattened indices (length sum cs_counts)

    char *err_msg; // optional null-terminated message (caller freed by ms_result_free)
} ms_result;

// C API entrypoint: runs susiex using in-memory arrays. Returns ssex_code.
// - npop_in: number of populations (sets global ::npop)
// - nsnp: number of SNPs
// - beta, pval: length npop_in * nsnp, row-major (pop-major)
// - ind: length npop_in * nsnp
// - ld: LDTYPE*** shaped [npop][nsnp][nsnp]
// - mkIdx: length nsnp
// - par_in: pointer to a filled softpar instance (caller-owned)
// - out: optional pointer to ms_result (caller-allocated)
int susiex_multisusie_fit(int npop_in,
                          int nsnp,
                          const double* beta,
                          const double* pval,
                          const char* ind,
                          LDTYPE*** ld,
                          const int* mkIdx,
                          const softpar* par_in,
                          ms_result* out);

void ms_result_free(ms_result* r);

#ifdef __cplusplus
}
#endif
