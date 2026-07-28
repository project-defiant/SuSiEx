"""Pure Python adapter around the compiled SuSiEx extension."""

from __future__ import annotations

from importlib import import_module
from typing import Any

import numpy as np

from .models import FitParameters


def fit(
    beta: np.ndarray,
    pval: np.ndarray,
    ind: np.ndarray,
    ld: np.ndarray,
    mk_idx: np.ndarray,
    *,
    n_gwas: list[int],
) -> dict[str, Any]:
    """Run SuSiEx through the native extension using validated NumPy arrays.

    The array contract mirrors ``susiex_python.susiex_pyfit`` while keeping
    extension imports and shape validation inside the Python package. File and
    dataset adapters will be added separately once the pipeline schema is
    finalized.
    """
    parameters = FitParameters(n_gwas=n_gwas)
    arrays = {
        "beta": np.asarray(beta, dtype=np.float64),
        "pval": np.asarray(pval, dtype=np.float64),
        "ind": np.asarray(ind, dtype=np.uint8),
        "ld": np.asarray(ld, dtype=np.float32),
        "mk_idx": np.asarray(mk_idx, dtype=np.int32),
    }
    _validate_shapes(arrays)
    if len(parameters.n_gwas) != arrays["beta"].shape[0]:
        raise ValueError("n_gwas length must equal the number of populations")

    native_module = import_module("susiex_python")
    susiex_pyfit = native_module.susiex_pyfit

    return susiex_pyfit(
        arrays["beta"],
        arrays["pval"],
        arrays["ind"],
        arrays["ld"],
        arrays["mk_idx"],
        parameters.n_gwas,
    )


def _validate_shapes(arrays: dict[str, np.ndarray]) -> None:
    beta = arrays["beta"]
    if beta.ndim != 2:
        raise ValueError("beta must have shape (n_populations, n_variants)")
    n_populations, n_variants = beta.shape
    expected = {
        "pval": (n_populations, n_variants),
        "ind": (n_populations, n_variants),
        "ld": (n_populations, n_variants, n_variants),
        "mk_idx": (n_variants,),
    }
    for name, shape in expected.items():
        if arrays[name].shape != shape:
            raise ValueError(
                f"{name} must have shape {shape}, got {arrays[name].shape}"
            )
