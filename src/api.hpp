#pragma once

#include "data.hpp"

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
