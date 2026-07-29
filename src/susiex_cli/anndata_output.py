"""Write native SuSiEx results as an extended AnnData object."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .preparation import PreparedArrays


def write_anndata(
    result: dict[str, Any],
    prepared: PreparedArrays,
    *,
    run_id: str,
    fine_mapping_locus_set_id: str,
    output: Path,
) -> None:
    """Persist all native components, posteriors, and provenance in H5AD."""

    alpha = np.asarray(result["alpha"], dtype=np.float32)
    n_components, n_variants = alpha.shape
    cs = result["cs"]
    obs = pd.DataFrame(
        {
            "componentIndex": np.arange(n_components),
            "logBF": np.asarray(result["lbf"], dtype=np.float32),
            "componentPurity": np.asarray(result["component_purity"], dtype=np.float32),
            "componentMinPValue": np.asarray(
                result["component_min_p_values"], dtype=np.float32
            ),
            "componentFiltered": np.asarray(
                result["component_is_filtered"], dtype=bool
            ),
            "credibleSetSize": [int(np.asarray(values).size) for values in cs],
        },
        index=[f"component_{index}" for index in range(n_components)],
    )
    var = pd.DataFrame(
        {
            "variantId": prepared.variant_ids,
            "chromosome": prepared.chromosomes,
            "position": prepared.positions,
            "pip": np.asarray(result["pip"], dtype=np.float32),
        },
        index=prepared.variant_ids,
    )
    mu = np.asarray(result["mu"], dtype=np.float32)
    layers = {
        f"mu__{ancestry}": mu[index]
        for index, ancestry in enumerate(prepared.ancestries)
        if index < mu.shape[0]
    }
    adata = ad.AnnData(X=alpha, obs=obs, var=var, layers=layers)
    adata.obsm["populationCausalProb"] = np.asarray(
        result["population_causal_prob"], dtype=np.float32
    )
    adata.obsm["logBFByPopulation"] = np.asarray(
        result["logbf_by_population"], dtype=np.float32
    )
    adata.varm["variantPresent"] = np.asarray(prepared.ind, dtype=bool).T
    adata.uns.update(
        {
            "runId": run_id,
            "fineMappingLocusSetId": fine_mapping_locus_set_id,
            "studyIds": prepared.study_ids,
            "ancestries": prepared.ancestries,
            "sampleSizes": prepared.sample_sizes,
            "converged": bool(result["converged"]),
            "inputMode": "summary-statistics",
        }
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
