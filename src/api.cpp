#include "susiex/api.hpp"
#include "susiex/validation.hpp"
#include "susiex/model.hpp"
#include "susiex/data.hpp"
#include <exception>
#include <new>
#include <cmath>
#include <algorithm>

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
    if(out)
        *out = ms_result{};

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
            out->converged = model.fit_status == 1;
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

            out->logbf_by_population = new double[(size_t)L * K * P];
            memcpy(out->logbf_by_population, model.logBF,
                   sizeof(double) * (size_t)L * K * P);
            out->cs_purity = new double[L];
            out->cs_min_p = new double[L];
            out->cs_min_p_by_population = new double[(size_t)L * K];
            out->population_causal_prob = new double[(size_t)L * K];
            out->cs_filtered = new int[L];
            for(int l = 0; l < L; ++l)
            {
                const bool filtered = l >= static_cast<int>(model.csset.size()) || model.csset[l].fltOut;
                out->cs_filtered[l] = filtered ? 1 : 0;
                out->cs_purity[l] = filtered ? 0.0 : model.csset[l].purity;
                out->cs_min_p[l] = filtered ? 1.0 : model.csset[l].minP;
                for(int k = 0; k < K; ++k)
                {
                    const size_t offset = (size_t)l * K + k;
                    out->cs_min_p_by_population[offset] =
                        filtered ? 1.0 : model.csset[l].minPalPop[k];
                    double max_logbf = -INFINITY;
                    for(int p = 0; p < P; ++p)
                        max_logbf = std::max(max_logbf, model.logBF[(size_t)l * K * P + (size_t)k * P + p]);
                    double total = 0.0;
                    for(int p = 0; p < P; ++p)
                        total += std::exp(model.logBF[(size_t)l * K * P + (size_t)k * P + p] - max_logbf);
                    out->population_causal_prob[offset] = total > 0.0 ? 1.0 : 0.0;
                }
            }

            // cs counts & indices
            out->cs_counts = new int[L];
            int total = 0;
            for(int l = 0; l < L; ++l)
            {
                if(l < static_cast<int>(model.csset.size()) && !model.csset[l].fltOut)
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
                for(int l = 0; l < L && l < static_cast<int>(model.csset.size()); ++l)
                {
                    if(model.csset[l].fltOut)
                        continue;
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
    if(r->logbf_by_population) delete [] r->logbf_by_population;
    if(r->cs_purity) delete [] r->cs_purity;
    if(r->cs_min_p) delete [] r->cs_min_p;
    if(r->cs_min_p_by_population) delete [] r->cs_min_p_by_population;
    if(r->population_causal_prob) delete [] r->population_causal_prob;
    if(r->cs_filtered) delete [] r->cs_filtered;
    if(r->cs_counts) delete [] r->cs_counts;
    if(r->cs_indices) delete [] r->cs_indices;
    if(r->err_msg) delete [] r->err_msg;
    // zero out to be safe
    r->alpha = nullptr;
    r->pip = nullptr;
    r->mu = nullptr;
    r->lbf = nullptr;
    r->logbf_by_population = nullptr;
    r->cs_purity = nullptr;
    r->cs_min_p = nullptr;
    r->cs_min_p_by_population = nullptr;
    r->population_causal_prob = nullptr;
    r->cs_filtered = nullptr;
    r->cs_counts = nullptr;
    r->cs_indices = nullptr;
    r->err_msg = nullptr;
}
