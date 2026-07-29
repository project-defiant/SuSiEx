"""Typer command line interface for the native SuSiEx adapter."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

import numpy as np
import typer
from loguru import logger

from .runner import fit
from .models import ApplicationInput
from .anndata_output import write_anndata
from .outputs import study_locus_records, write_study_locus_parquet
from .preparation import prepare_arrays
from .stats import SuSiExStats, stats_from_result, write_stats

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
    run_id: str = typer.Option("array-run"),
    fine_mapping_locus_set_id: str = typer.Option("array-locus"),
    stats_output: Path | None = typer.Option(None),
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
        if stats_output is not None:
            write_stats(
                stats_output,
                SuSiExStats(
                    runId=run_id,
                    fineMappingLocusSetId=fine_mapping_locus_set_id,
                    status="FAILED",
                    reason=str(error),
                ),
            )
        logger.error("SuSiEx run failed: {}", error)
        raise typer.Exit(code=1) from error

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(_summary(result), sort_keys=True) + "\n")
    if stats_output is not None:
        write_stats(
            stats_output,
            stats_from_result(
                result,
                run_id=run_id,
                fine_mapping_locus_set_id=fine_mapping_locus_set_id,
            ),
        )


@app.command("pipeline")
def pipeline_run(
    fine_mapping_locus_set: Path = typer.Option(..., exists=True, readable=True),
    multi_ancestry_pairwise_ld: Path = typer.Option(..., exists=True, readable=True),
    study_metadata: Path = typer.Option(..., exists=True, readable=True),
    run_id: str = typer.Option(...),
    fine_mapping_locus_set_id: str = typer.Option(...),
    study_locus_output: Path = typer.Option(...),
    extended_results_output: Path = typer.Option(...),
    stats_output: Path = typer.Option(...),
    n_sig: int = typer.Option(5, min=1),
    max_iter: int = typer.Option(100, min=1),
    level: float = typer.Option(0.95, min=0.000001, max=0.999999),
    min_purity: float = typer.Option(0.5, min=0, max=1),
    pth: float = typer.Option(1e-5, min=0.000000000001, max=1),
    tol: float = typer.Option(1e-4, min=0.000000000001),
    nthreads: int = typer.Option(1, min=1),
    mult_step: bool = typer.Option(False),
) -> None:
    """Run SuSiEx from pipeline parquet and JSONL inputs."""

    inputs = ApplicationInput(
        run_id=run_id,
        fine_mapping_locus_set_id=fine_mapping_locus_set_id,
        fine_mapping_locus_set_path=fine_mapping_locus_set,
        multi_ancestry_pairwise_ld_path=multi_ancestry_pairwise_ld,
        study_metadata_path=study_metadata,
        study_locus_output_path=study_locus_output,
        extended_results_output_path=extended_results_output,
        stats_output_path=stats_output,
    )
    try:
        prepared = prepare_arrays(inputs)
        result = fit(
            prepared.beta,
            prepared.pval,
            prepared.ind,
            prepared.ld,
            prepared.mk_idx,
            n_gwas=prepared.sample_sizes,
            n_sig=n_sig,
            max_iter=max_iter,
            level=level,
            min_purity=min_purity,
            pth=pth,
            tol=tol,
            nthreads=nthreads,
            mult_step=mult_step,
        )
    except (OSError, ValueError, RuntimeError) as error:
        write_stats(
            stats_output,
            SuSiExStats(
                runId=run_id,
                fineMappingLocusSetId=fine_mapping_locus_set_id,
                status="FAILED",
                reason=str(error),
            ),
        )
        logger.error("SuSiEx pipeline run failed: {}", error)
        raise typer.Exit(code=1) from error

    stats = stats_from_result(
        result, run_id=run_id, fine_mapping_locus_set_id=fine_mapping_locus_set_id
    )
    write_stats(stats_output, stats)
    if not bool(result["converged"]):
        logger.warning("SuSiEx pipeline fit did not converge")
        return

    variants = [
        {"variantId": variant_id, "chromosome": chromosome, "position": position}
        for variant_id, chromosome, position in zip(
            prepared.variant_ids,
            prepared.chromosomes,
            prepared.positions,
            strict=True,
        )
    ]
    records = study_locus_records(
        result,
        variants,
        run_id=run_id,
        fine_mapping_locus_set_id=fine_mapping_locus_set_id,
        study_id="|".join(prepared.study_ids),
        chromosome=prepared.chromosomes[0],
        locus_start=min(prepared.positions),
        locus_end=max(prepared.positions),
    )
    write_study_locus_parquet(records, study_locus_output)
    write_anndata(
        result,
        prepared,
        run_id=run_id,
        fine_mapping_locus_set_id=fine_mapping_locus_set_id,
        output=extended_results_output,
    )


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
