"""Validated models for the first Python/native SuSiEx boundary."""

from __future__ import annotations

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
