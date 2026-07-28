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
        {
            // populate result
            out->converged = 1; // assume success if we reached here
            out->L = model.nsig;
            out->ncs = model.ncs;
            out->K = model.npop;
            out->P = model.nsnp;
            int L = out->L;
            int K = out->K;
            int P = out->P;

            // alpha: model.alpha is size L * P (indexed [l*P + p])
            out->alpha = new double[(size_t)L * P];
            memcpy(out->alpha, model.alpha, sizeof(double) * (size_t)L * P);

            // pip
            out->pip = new double[P];
            memcpy(out->pip, model.pip, sizeof(double) * P);

            // mu: use model.b (size L * K * P) with layout b[l * (K*P) + k*P + p]
            out->mu = new double[(size_t)K * L * P];
            for(int l = 0; l < L; ++l)
                for(int k = 0; k < K; ++k)
                    for(int p = 0; p < P; ++p)
                    {
                        size_t src_idx = (size_t)l * (K * P) + k * P + p;
                        size_t dst_idx = (size_t)l * (K * P) + k * P + p;
                        out->mu[dst_idx] = model.b[src_idx];
                    }

            // lbf: use model.loglik (length L)
            out->lbf = new double[(size_t)L];
            memcpy(out->lbf, model.loglik, sizeof(double) * L);

            // cs counts & indices
            out->cs_counts = new int[L];
            int total = 0;
            for(int l = 0; l < L; ++l)
            {
                if(!model.csset[l].fltOut)
                    out->cs_counts[l] = (int)model.csset[l].idx.size();
                else
                    out->cs_counts[l] = 0;
                total += out->cs_counts[l];
            }
            out->cs_indices = nullptr;
            if(total > 0)
            {
                out->cs_indices = new int[total];
                int pos = 0;
                for(int l = 0; l < L; ++l)
                {
                    for(size_t j = 0; j < model.csset[l].idx.size(); ++j)
                    {
                        out->cs_indices[pos++] = model.csset[l].idx[j];
                    }
                }
            }

            out->err_msg = nullptr;
        }

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
    if(!r) return;
    if(r->alpha) delete [] r->alpha;
    if(r->pip) delete [] r->pip;
    if(r->mu) delete [] r->mu;
    if(r->lbf) delete [] r->lbf;
    if(r->cs_counts) delete [] r->cs_counts;
    if(r->cs_indices) delete [] r->cs_indices;
    if(r->err_msg) delete [] r->err_msg;
    // zero out to be safe
    r->alpha = nullptr;
    r->pip = nullptr;
    r->mu = nullptr;
    r->lbf = nullptr;
    r->cs_counts = nullptr;
    r->cs_indices = nullptr;
    r->err_msg = nullptr;
}
