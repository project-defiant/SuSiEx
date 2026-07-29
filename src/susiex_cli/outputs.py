"""Conversion of native SuSiEx results to StudyLocus records."""

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

STUDY_LOCUS_COLUMNS = (
    "studyLocusId",
    "studyId",
    "variantId",
    "chromosome",
    "position",
    "beta",
    "sampleSize",
    "pValueMantissa",
    "pValueExponent",
    "effectAlleleFrequencyFromSource",
    "standardError",
    "qualityControls",
    "locusStart",
    "locusEnd",
    "locus",
)


def study_locus_records(
    result: Mapping[str, Any],
    variants: Sequence[Mapping[str, Any]],
    *,
    run_id: str,
    fine_mapping_locus_set_id: str,
    study_id: str,
    chromosome: str,
    locus_start: int,
    locus_end: int,
) -> list[dict[str, Any]]:
    """Build one shared StudyLocus record per native credible set."""

    alpha = np.asarray(result["alpha"])
    log_bf = np.asarray(result["lbf"])
    records: list[dict[str, Any]] = []
    for component, indices in enumerate(result["cs"]):
        ordered_indices = [int(index) for index in np.asarray(indices)]
        if not ordered_indices:
            continue
        if any(index < 0 or index >= len(variants) for index in ordered_indices):
            raise ValueError("Credible-set variant index is outside the variant table")

        probabilities = [float(alpha[component, index]) for index in ordered_indices]
        total = sum(probabilities)
        if total > 0:
            probabilities = [probability / total for probability in probabilities]
        cumulative = 0.0
        nested: list[dict[str, Any]] = []
        for index, probability in zip(ordered_indices, probabilities, strict=True):
            cumulative += probability
            variant = dict(variants[index])
            nested.append(
                {
                    "is95CredibleSet": cumulative >= 0.95,
                    "is99CredibleSet": cumulative >= 0.99,
                    "logBF": float(log_bf[component]),
                    "posteriorProbability": probability,
                    "variantId": variant["variantId"],
                    "pValueMantissa": variant.get("pValueMantissa"),
                    "pValueExponent": variant.get("pValueExponent"),
                    "beta": variant.get("beta"),
                    "standardError": variant.get("standardError"),
                    "r2Overall": None,
                }
            )

        lead_offset = int(np.argmax(probabilities))
        lead = dict(variants[ordered_indices[lead_offset]])
        stable_key = f"{run_id}:{fine_mapping_locus_set_id}:susiex:{component}:{lead['variantId']}"
        records.append(
            {
                "studyLocusId": hashlib.md5(stable_key.encode()).hexdigest(),
                "studyId": study_id,
                "variantId": lead["variantId"],
                "chromosome": chromosome,
                "position": lead.get("position"),
                "beta": lead.get("beta"),
                "sampleSize": lead.get("sampleSize"),
                "pValueMantissa": lead.get("pValueMantissa"),
                "pValueExponent": lead.get("pValueExponent"),
                "effectAlleleFrequencyFromSource": lead.get(
                    "effectAlleleFrequencyFromSource"
                ),
                "standardError": lead.get("standardError"),
                "qualityControls": [],
                "locusStart": locus_start,
                "locusEnd": locus_end,
                "locus": nested,
            }
        )
    return records


def write_study_locus_parquet(records: list[dict[str, Any]], output: Path) -> None:
    """Write shared StudyLocus records as one flat parquet file."""

    if not records:
        raise ValueError("Cannot write StudyLocus output without records")
    output.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(records).write_parquet(output)
