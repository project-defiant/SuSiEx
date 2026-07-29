"""Machine-readable SuSiEx run statistics."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .models import SuSiExStats


def stats_from_result(
    result: dict[str, Any], *, run_id: str, fine_mapping_locus_set_id: str
) -> SuSiExStats:
    """Summarize a native result without serializing numerical arrays."""

    alpha = np.asarray(result["alpha"])
    pip = np.asarray(result["pip"])
    credible_sets = result["cs"]
    converged = bool(result["converged"])
    return SuSiExStats(
        runId=run_id,
        fineMappingLocusSetId=fine_mapping_locus_set_id,
        status="SUCCESS" if converged else "NON_CONVERGED",
        converged=converged,
        nComponents=int(alpha.shape[0]),
        nPopulations=int(np.asarray(result["mu"]).shape[0]),
        nVariants=int(pip.shape[0]),
        nCredibleSets=len(credible_sets),
        reason=None if converged else "SuSiEx fit did not converge",
    )


def write_stats(output: Path, stats: SuSiExStats) -> None:
    """Write one deterministic JSON status record."""

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(stats.model_dump(exclude_none=True), sort_keys=True) + "\n")
