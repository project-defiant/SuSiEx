import numpy as np
import susiex_python

# minimal smoke test (does not run here; for CI when pybind env available)

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

    rc = susiex_python.susiex_pyfit(beta, pval, ind, ld, mkIdx)
    assert rc == 1  # SSEX_VALIDATION_ERROR expected
