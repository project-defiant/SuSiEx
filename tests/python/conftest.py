import numpy as np
import pytest


@pytest.fixture
def valid_fit_arrays() -> dict[str, np.ndarray]:
    n_populations = 1
    n_variants = 5
    beta = np.zeros((n_populations, n_variants), dtype=np.float64)
    beta[0, 0] = 0.2
    pval = np.ones((n_populations, n_variants), dtype=np.float64)
    pval[0, 0] = 1e-8
    ind = np.ones((n_populations, n_variants), dtype=np.uint8)
    ld = np.zeros((n_populations, n_variants, n_variants), dtype=np.float32)
    for i in range(n_variants):
        for j in range(n_variants):
            ld[0, i, j] = 1.0 if i == j else (0.1 if abs(i - j) == 1 else 0.01)
    return {
        "beta": beta,
        "pval": pval,
        "ind": ind,
        "ld": ld,
        "mk_idx": np.arange(n_variants, dtype=np.int32),
    }
