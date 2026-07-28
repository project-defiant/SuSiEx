"""Validated models for the first Python/native SuSiEx boundary."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, field_validator


class FitParameters(BaseModel):
    """Parameters forwarded to the native SuSiEx implementation."""

    model_config = ConfigDict(extra="forbid")

    n_gwas: list[int] = Field(min_length=1)

    @field_validator("n_gwas")
    @classmethod
    def sample_sizes_are_positive(cls, value: list[int]) -> list[int]:
        if any(sample_size <= 0 for sample_size in value):
            raise ValueError("n_gwas values must be positive")
        return value


class FitResult(BaseModel):
    """Serializable scalar and array-shaped result metadata."""

    model_config = ConfigDict(extra="forbid")

    converged: bool
    n_components: int = Field(ge=0)
    n_populations: int = Field(gt=0)
    n_variants: int = Field(gt=0)
    credible_set_sizes: list[int]


class ApplicationInput(BaseModel):
    """File and identity contract for one fine-mapping locus-set run."""

    model_config = ConfigDict(extra="forbid")

    run_id: str = Field(min_length=1)
    fine_mapping_locus_set_id: str = Field(min_length=1)
    fine_mapping_locus_set_path: Path
    multi_ancestry_pairwise_ld_path: Path
    study_metadata_path: Path
    study_locus_output_path: Path
    extended_results_output_path: Path
    stats_output_path: Path

    @field_validator("run_id", "fine_mapping_locus_set_id")
    @classmethod
    def identifiers_are_not_blank(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("identifiers must not be blank")
        return value
