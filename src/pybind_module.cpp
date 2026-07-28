#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "api.hpp"

namespace py = pybind11;

// Minimal convenience wrapper: uses default softpar with provided n_gwas per population or fallback 10000
int susiex_pyfit(py::array_t<double> beta_in,
                 py::array_t<double> pval_in,
                 py::array_t<uint8_t> ind_in,
                 py::array_t<float> ld_in,
                 py::array_t<int> mkIdx_in,
                 std::vector<int> n_gwas = std::vector<int>())
{
    // Expect shapes: beta (npop, nsnp), pval (npop, nsnp), ind (npop, nsnp), ld (npop, nsnp, nsnp), mkIdx (nsnp,)
    if(beta_in.ndim() != 2 || pval_in.ndim() != 2 || ind_in.ndim() != 2 || ld_in.ndim() != 3 || mkIdx_in.ndim() != 1)
        throw std::invalid_argument("Unexpected array dimensions");

    ssize_t npop_s = beta_in.shape(0);
    ssize_t nsnp_s = beta_in.shape(1);
    int npop = (int)npop_s;
    int nsnp = (int)nsnp_s;

    if(pval_in.shape(0) != npop_s || pval_in.shape(1) != nsnp_s) throw std::invalid_argument("pval shape mismatch");
    if(ind_in.shape(0) != npop_s || ind_in.shape(1) != nsnp_s) throw std::invalid_argument("ind shape mismatch");
    if(ld_in.shape(0) != npop_s || ld_in.shape(1) != nsnp_s || ld_in.shape(2) != nsnp_s) throw std::invalid_argument("ld shape mismatch");
    if(mkIdx_in.shape(0) != nsnp_s) throw std::invalid_argument("mkIdx length mismatch");

    // copy beta and pval into contiguous vectors (pop-major)
    auto beta_buf = beta_in.unchecked<2>();
    auto pval_buf = pval_in.unchecked<2>();
    std::vector<double> beta((size_t)npop * nsnp);
    std::vector<double> pval((size_t)npop * nsnp);
    for(int i = 0; i < npop; ++i)
        for(int j = 0; j < nsnp; ++j)
        {
            beta[(size_t)i * nsnp + j] = beta_buf(i, j);
            pval[(size_t)i * nsnp + j] = pval_buf(i, j);
        }

    // copy ind (uint8 -> char)
    auto ind_buf = ind_in.unchecked<2>();
    std::vector<char> ind((size_t)npop * nsnp);
    for(int i = 0; i < npop; ++i)
        for(int j = 0; j < nsnp; ++j)
            ind[(size_t)i * nsnp + j] = static_cast<char>(ind_buf(i, j));

    // copy mkIdx
    auto mkBuf = mkIdx_in.unchecked<1>();
    std::vector<int> mkIdx((size_t)nsnp);
    for(int j = 0; j < nsnp; ++j) mkIdx[j] = mkBuf(j);

    // build ld triple pointer
    auto ld_buf = ld_in.unchecked<3>();
    LDTYPE*** ld = new LDTYPE**[npop];
    for(int i = 0; i < npop; ++i)
    {
        ld[i] = new LDTYPE*[nsnp];
        for(int j = 0; j < nsnp; ++j)
        {
            ld[i][j] = new LDTYPE[nsnp];
            for(int k = 0; k < nsnp; ++k)
                ld[i][j][k] = ld_buf(i, j, k);
        }
    }

    // prepare softpar
    softpar par;
    if((int)n_gwas.size() == npop)
    {
        par.n_gwas = n_gwas;
    }
    else
    {
        par.n_gwas.clear();
        for(int i = 0; i < npop; ++i) par.n_gwas.push_back(10000);
    }

    ms_result out;
    int rc = susiex_multisusie_fit(npop, nsnp, beta.data(), pval.data(), ind.data(), ld, mkIdx.data(), &par, &out);

    // free ld
    for(int i = 0; i < npop; ++i)
    {
        for(int j = 0; j < nsnp; ++j)
            delete [] ld[i][j];
        delete [] ld[i];
    }
    delete [] ld;

    if(rc != SSEX_OK)
    {
        ms_result_free(&out);
        throw std::runtime_error(std::string("susiex_multisusie_fit failed with code: ") + std::to_string(rc));
    }

    // map ms_result to Python objects
    py::dict res;
    int L = out.L;
    int K = out.K;
    int P = out.P;

    // alpha (L, P)
    py::array_t<double> alpha({(size_t)L, (size_t)P});
    auto alpha_buf = alpha.mutable_unchecked<2>();
    for(int l = 0; l < L; ++l)
        for(int p = 0; p < P; ++p)
            alpha_buf(l, p) = out.alpha[(size_t)l * P + p];

    // pip (P,)
    py::array_t<double> pip({(size_t)P});
    auto pip_buf = pip.mutable_unchecked<1>();
    for(int p = 0; p < P; ++p) pip_buf(p) = out.pip[p];

    // mu (K, L, P) - out->mu stored as l*(K*P) + k*P + p
    py::array_t<double> mu({(size_t)K, (size_t)L, (size_t)P});
    auto mu_buf = mu.mutable_unchecked<3>();
    for(int l = 0; l < L; ++l)
        for(int k = 0; k < K; ++k)
            for(int p = 0; p < P; ++p)
            {
                size_t idx = (size_t)l * (K * P) + k * P + p;
                mu_buf(k, l, p) = out.mu[idx];
            }

    // lbf
    py::array_t<double> lbf({(size_t)L});
    auto lbf_buf = lbf.mutable_unchecked<1>();
    for(int l = 0; l < L; ++l) lbf_buf(l) = out.lbf[l];

    // cs: list of arrays
    py::list cs_list;
    int pos = 0;
    for(int l = 0; l < L; ++l)
    {
        int cnt = out.cs_counts ? out.cs_counts[l] : 0;
        py::array_t<int> arr({(size_t)cnt});
        auto arr_buf = arr.mutable_unchecked<1>();
        for(int j = 0; j < cnt; ++j)
            arr_buf(j) = out.cs_indices[pos++];
        cs_list.append(arr);
    }

    res["alpha"] = alpha;
    res["pip"] = pip;
    res["mu"] = mu;
    res["lbf"] = lbf;
    res["cs"] = cs_list;
    res["converged"] = true;

    ms_result_free(&out);
    return res;
}

PYBIND11_MODULE(susiex_python, m) {
    m.doc() = "Pybind11 wrapper for SuSiEx minimal C API";
    m.def("susiex_pyfit", &susiex_pyfit, py::arg("beta"), py::arg("pval"), py::arg("ind"), py::arg("ld"), py::arg("mkIdx"), py::arg("n_gwas") = std::vector<int>());
}

PYBIND11_MODULE(susiex_python, m) {
    m.doc() = "Pybind11 wrapper for SuSiEx minimal C API";
    m.def("susiex_pyfit", &susiex_pyfit, py::arg("beta"), py::arg("pval"), py::arg("ind"), py::arg("ld"), py::arg("mkIdx"), py::arg("n_gwas") = std::vector<int>());
}
