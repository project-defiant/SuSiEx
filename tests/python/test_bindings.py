import numpy as np
import pytest
import susiex_python

# validation error should raise a runtime error from the binding

def test_susiex_pyfit_validation_error():
    npop = 1
    nsnp = 3
    beta = np.zeros((npop, nsnp), dtype=np.float64)
    beta[0,0] = 0.1
    pval = np.ones((npop, nsnp), dtype=np.float64)
    pval[0,0] = 1e-6
    ind = np.ones((npop, nsnp), dtype=np.uint8)
    ld = np.zeros((npop, nsnp, nsnp), dtype=np.float32)
    for i in range(nsnp):
        for j in range(nsnp):
            ld[0,i,j] = 1.0 if i==j else 0.05
    # introduce asymmetry
    ld[0,0,1] = 0.1
    ld[0,1,0] = 0.2
    mkIdx = np.arange(nsnp, dtype=np.int32)

    with pytest.raises(RuntimeError):
        susiex_python.susiex_pyfit(beta, pval, ind, ld, mkIdx)


def test_susiex_pyfit_success_smoke():
    npop = 1
    nsnp = 5
    beta = np.zeros((npop, nsnp), dtype=np.float64)
    # causal at index 0
    beta[0,0] = 0.2
    pval = np.ones((npop, nsnp), dtype=np.float64)
    pval[0,0] = 1e-8
    ind = np.ones((npop, nsnp), dtype=np.uint8)
    ld = np.zeros((npop, nsnp, nsnp), dtype=np.float32)
    for i in range(nsnp):
        for j in range(nsnp):
            ld[0,i,j] = 1.0 if i==j else (0.1 if abs(i-j)==1 else 0.01)
    mkIdx = np.arange(nsnp, dtype=np.int32)

    res = susiex_python.susiex_pyfit(beta, pval, ind, ld, mkIdx)

    assert isinstance(res, dict)
    assert "alpha" in res and "pip" in res and "mu" in res and "cs" in res
    pip = res["pip"]
    # expect the causal SNP to have the highest pip
    assert pip[0] >= max(pip[1:])
