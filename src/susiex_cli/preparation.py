"""Prepare pipeline parquet inputs for the native SuSiEx API."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict

from .models import ApplicationInput


class PreparedArrays(BaseModel):
    """Deterministically ordered arrays accepted by the native binding."""

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    variant_ids: list[str]
    chromosomes: list[str]
    positions: list[int]
    beta: np.ndarray
    pval: np.ndarray
    ind: np.ndarray
    ld: np.ndarray
    mk_idx: np.ndarray


def prepare_arrays(inputs: ApplicationInput) -> PreparedArrays:
    """Read and align locus, LD, and metadata files for one application run."""

    locus = pl.read_parquet(_parquet_path(inputs.fine_mapping_locus_set_path))
    metadata = _read_metadata(inputs.study_metadata_path)
    _validate_locus(locus, inputs.fine_mapping_locus_set_id, metadata)
    rows = _variant_rows(locus)
    ordered = sorted(
        (variant for study in rows.values() for variant in study.values()),
        key=lambda row: (row["chromosome"], row["position"], row["variant_id"]),
    )
    by_id = {row["variant_id"]: row for row in ordered}
    variant_ids = list(by_id)
    index = {variant_id: i for i, variant_id in enumerate(variant_ids)}
    ld = pl.read_parquet(_parquet_path(inputs.multi_ancestry_pairwise_ld_path))
    _validate_ld(ld)

    beta = np.zeros((len(metadata), len(variant_ids)), dtype=np.float64)
    pval = np.ones_like(beta)
    ind = np.zeros((len(metadata), len(variant_ids)), dtype=np.uint8)
    ld_matrices = np.zeros(
        (len(metadata), len(variant_ids), len(variant_ids)), dtype=np.float32
    )
    for population, metadata_row in enumerate(metadata):
        study = rows[metadata_row["studyId"]]
        for variant_id, variant in study.items():
            i = index[variant_id]
            beta[population, i] = variant["beta"]
            pval[population, i] = variant["p_value"]
            ind[population, i] = 1
        np.fill_diagonal(ld_matrices[population], ind[population])
        _fill_ld(
            ld_matrices[population],
            ld.filter(pl.col("ancestry") == metadata_row["ancestry"]),
            index,
        )

    return PreparedArrays(
        variant_ids=variant_ids,
        chromosomes=[str(row["chromosome"]) for row in ordered],
        positions=[int(row["position"]) for row in ordered],
        beta=beta,
        pval=pval,
        ind=ind,
        ld=ld_matrices,
        mk_idx=np.arange(len(variant_ids), dtype=np.int32),
    )


def _parquet_path(path: Path) -> Path | str:
    return str(path / "**" / "*.parquet") if path.is_dir() else path


def _read_metadata(path: Path) -> list[dict[str, Any]]:
    records = [json.loads(line) for line in path.read_text().splitlines() if line]
    if not records:
        raise ValueError("Study metadata must contain at least one row")
    if len({row.get("studyId") for row in records}) != len(records):
        raise ValueError("Study metadata contains duplicate studyId values")
    if len({row.get("ancestry") for row in records}) != len(records):
        raise ValueError("Study metadata must contain one study per ancestry")
    if any(int(row.get("sampleSize", 0)) <= 0 for row in records):
        raise ValueError("Study metadata sampleSize values must be positive")
    return sorted(records, key=lambda row: str(row["studyId"]))


def _validate_locus(
    locus: pl.DataFrame, locus_set_id: str, metadata: list[dict[str, Any]]
) -> None:
    required = {"fineMappingLocusSetId", "studyId", "locus"}
    missing = required - set(locus.columns)
    if missing:
        raise ValueError(f"FineMappingLocusSet is missing columns: {sorted(missing)}")
    observed_ids = set(locus["fineMappingLocusSetId"].drop_nulls().to_list())
    if observed_ids != {locus_set_id}:
        raise ValueError("FineMappingLocusSet contains an unexpected locus-set ID")
    observed_studies = set(locus["studyId"].drop_nulls().to_list())
    expected_studies = {row["studyId"] for row in metadata}
    if observed_studies != expected_studies:
        raise ValueError("Study metadata and locus studyId values differ")


def _variant_rows(locus: pl.DataFrame) -> dict[str, dict[str, dict[str, Any]]]:
    rows: dict[str, dict[str, dict[str, Any]]] = {}
    for row in locus.select("studyId", "locus").iter_rows(named=True):
        study_id = str(row["studyId"])
        if row["locus"] is None:
            raise ValueError(f"Study locus has no variants: {study_id}")
        study = rows.setdefault(study_id, {})
        for variant in row["locus"]:
            variant_id = str(variant["variantId"])
            if variant_id in study:
                raise ValueError(f"Duplicate variant in study locus: {variant_id}")
            beta = variant.get("beta")
            se = variant.get("standardError")
            if (
                beta is None
                or se is None
                or not math.isfinite(float(beta))
                or not math.isfinite(float(se))
                or float(se) <= 0
            ):
                raise ValueError(
                    f"Invalid beta or standardError for {study_id}/{variant_id}"
                )
            chromosome, position = _parse_variant_id(variant_id)
            mantissa = variant.get("pValueMantissa")
            exponent = variant.get("pValueExponent")
            p_value = (
                1.0
                if mantissa is None or exponent is None
                else float(mantissa) * 10.0 ** int(exponent)
            )
            if not math.isfinite(p_value) or p_value <= 0:
                raise ValueError(f"Invalid p-value for {study_id}/{variant_id}")
            study[variant_id] = {
                "variant_id": variant_id,
                "chromosome": chromosome,
                "position": position,
                "beta": float(beta),
                "p_value": min(p_value, 1.0),
            }
    return rows


def _validate_ld(ld: pl.DataFrame) -> None:
    required = {"ancestry", "variantIdI", "variantIdJ", "r"}
    missing = required - set(ld.columns)
    if missing:
        raise ValueError(
            f"MultiAncestryPairwiseLD is missing columns: {sorted(missing)}"
        )


def _fill_ld(matrix: np.ndarray, rows: pl.DataFrame, index: dict[str, int]) -> None:
    for row in rows.iter_rows(named=True):
        first = str(row["variantIdI"])
        second = str(row["variantIdJ"])
        if first not in index or second not in index:
            continue
        value = float(row["r"])
        if not math.isfinite(value) or not -1 <= value <= 1:
            raise ValueError(f"Invalid LD value for {first}/{second}: {value}")
        i, j = index[first], index[second]
        matrix[i, j] = value
        matrix[j, i] = value


def _parse_variant_id(variant_id: str) -> tuple[str, int]:
    fields = variant_id.split("_")
    if len(fields) < 2:
        raise ValueError(
            f"Cannot parse chromosome and position from variantId: {variant_id}"
        )
    try:
        return fields[0], int(fields[1])
    except ValueError as error:
        raise ValueError(
            f"Cannot parse position from variantId: {variant_id}"
        ) from error
