"""Typer command line interface for the native SuSiEx adapter."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

import numpy as np
import typer
from loguru import logger

from .runner import fit

app = typer.Typer(no_args_is_help=True)


@app.command()
def run(
    beta: Path = typer.Option(..., exists=True, readable=True),
    pval: Path = typer.Option(..., exists=True, readable=True),
    ind: Path = typer.Option(..., exists=True, readable=True),
    ld: Path = typer.Option(..., exists=True, readable=True),
    mk_idx: Path = typer.Option(..., exists=True, readable=True),
    n_gwas: str = typer.Option(..., help="Comma-separated sample sizes by population."),
    output: Path = typer.Option(
        ..., help="JSON output for the native result metadata."
    ),
) -> None:
    """Run SuSiEx from NumPy arrays; dataset adapters are a separate layer."""
    try:
        result = fit(
            np.load(beta),
            np.load(pval),
            np.load(ind),
            np.load(ld),
            np.load(mk_idx),
            n_gwas=[int(value) for value in n_gwas.split(",")],
        )
    except (OSError, ValueError, RuntimeError) as error:
        logger.error("SuSiEx run failed: {}", error)
        raise typer.Exit(code=1) from error

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(_summary(result), sort_keys=True) + "\n")


def _summary(result: dict[str, object]) -> dict[str, object]:
    """Reduce native arrays to a small machine-readable run summary."""
    alpha = np.asarray(result["alpha"])
    pip = np.asarray(result["pip"])
    cs = cast(list[Any], result["cs"])
    return {
        "converged": bool(result["converged"]),
        "n_components": int(alpha.shape[0]),
        "n_populations": int(np.asarray(result["mu"]).shape[0]),
        "n_variants": int(pip.shape[0]),
        "credible_set_sizes": [int(np.asarray(values).size) for values in cs],
    }
