import numpy as np
import pytest

from susiex_cli.models import FitParameters
from susiex_cli.runner import fit


@pytest.mark.parametrize(
    ("name", "value", "message"),
    [
        ("pval", np.zeros((1, 2)), "pval must have shape"),
        ("ind", np.zeros((1, 2), dtype=np.uint8), "ind must have shape"),
        ("ld", np.zeros((1, 3, 2), dtype=np.float32), "ld must have shape"),
        ("mk_idx", np.arange(2), "mk_idx must have shape"),
    ],
)
def test_fit_rejects_inconsistent_array_shapes(name, value, message) -> None:
    arrays = {
        "beta": np.zeros((1, 3)),
        "pval": np.zeros((1, 3)),
        "ind": np.ones((1, 3), dtype=np.uint8),
        "ld": np.eye(3, dtype=np.float32)[None, :, :],
        "mk_idx": np.arange(3),
    }
    arrays[name] = value
    with pytest.raises(ValueError, match=message):
        fit(**arrays, n_gwas=[100])


def test_fit_rejects_sample_size_population_mismatch() -> None:
    with pytest.raises(ValueError, match="n_gwas length"):
        fit(
            np.zeros((2, 3)),
            np.ones((2, 3)),
            np.ones((2, 3), dtype=np.uint8),
            np.broadcast_to(np.eye(3, dtype=np.float32), (2, 3, 3)),
            np.arange(3),
            n_gwas=[100],
        )


@pytest.mark.parametrize("sample_sizes", [[], [0], [-1]])
def test_fit_parameters_reject_invalid_sample_sizes(sample_sizes) -> None:
    with pytest.raises(ValueError):
        FitParameters(n_gwas=sample_sizes)
