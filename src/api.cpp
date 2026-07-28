#include "api.hpp"
#include "validation.hpp"
#include "model.hpp"
#include "data.hpp"
#include <exception>
#include <new>

extern int npop; // declared in data.hpp

int susiex_multisusie_fit(int npop_in,
                          int nsnp,
                          const double* beta,
                          const double* pval,
                          const char* ind,
                          LDTYPE*** ld,
                          const int* mkIdx,
                          const softpar* par_in,
                          ms_result* out)
{
    if(npop_in <= 0 || nsnp <= 0 || !beta || !pval || !ind || !ld || !mkIdx || !par_in)
        return SSEX_INVALID_ARGS;

    try
    {
        // set global
        ::npop = npop_in;

        // copy softpar locally
        softpar par = *par_in;

        dataset dat;
        try {
            dat.load_from_memory(par, nsnp, beta, pval, ind, ld, mkIdx);
        } catch(const std::exception &e) {
            // validation errors are surfaced as runtime_error in loader
            return SSEX_VALIDATION_ERROR;
        }

        susiex model(::npop, dat.nsnp, dat, par);
        model.susie_sst_xethn();

        if(out)
            out->converged = 1; // currently no detailed diagnostics

        return SSEX_OK;
    }
    catch(const std::bad_alloc &)
    {
        return SSEX_INTERNAL_ERROR;
    }
    catch(const std::exception &)
    {
        return SSEX_INTERNAL_ERROR;
    }
    catch(...)
    {
        return SSEX_INTERNAL_ERROR;
    }
}

void ms_result_free(ms_result* r)
{
    // no heap allocations currently; placeholder for future fields
    (void)r;
}
