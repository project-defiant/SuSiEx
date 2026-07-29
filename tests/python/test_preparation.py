import json

import polars as pl

from susiex_cli.models import ApplicationInput
from susiex_cli.preparation import prepare_arrays


def test_prepare_arrays_aligns_union_and_missing_population(tmp_path) -> None:
    locus_path = tmp_path / "locus.parquet"
    pl.DataFrame(
        {
            "fineMappingLocusSetId": ["LOCUS_A", "LOCUS_A"],
            "studyId": ["STUDY_A", "STUDY_B"],
            "locus": [
                [
                    {
                        "variantId": "1_110_A_G",
                        "beta": 2.0,
                        "standardError": 1.0,
                        "pValueMantissa": 2.0,
                        "pValueExponent": -3,
                    },
                    {
                        "variantId": "1_100_C_T",
                        "beta": 1.0,
                        "standardError": 1.0,
                        "pValueMantissa": 1.0,
                        "pValueExponent": -2,
                    },
                ],
                [
                    {
                        "variantId": "1_110_A_G",
                        "beta": 3.0,
                        "standardError": 1.5,
                        "pValueMantissa": 3.0,
                        "pValueExponent": -4,
                    }
                ],
            ],
        }
    ).write_parquet(locus_path)
    ld_path = tmp_path / "ld.parquet"
    pl.DataFrame(
        {
            "ancestry": ["eur"],
            "variantIdI": ["1_100_C_T"],
            "variantIdJ": ["1_110_A_G"],
            "r": [0.25],
        }
    ).write_parquet(ld_path)
    metadata = tmp_path / "metadata.jsonl"
    metadata.write_text(
        "\n".join(
            json.dumps(row)
            for row in [
                {"studyId": "STUDY_A", "ancestry": "eur", "sampleSize": 1000},
                {"studyId": "STUDY_B", "ancestry": "eur2", "sampleSize": 2000},
            ]
        )
        + "\n"
    )
    prepared = prepare_arrays(
        ApplicationInput(
            run_id="RUN_A",
            fine_mapping_locus_set_id="LOCUS_A",
            fine_mapping_locus_set_path=locus_path,
            multi_ancestry_pairwise_ld_path=ld_path,
            study_metadata_path=metadata,
            study_locus_output_path=tmp_path / "study.parquet",
            extended_results_output_path=tmp_path / "fit.h5ad",
            stats_output_path=tmp_path / "stats.json",
        )
    )

    assert prepared.variant_ids == ["1_100_C_T", "1_110_A_G"]
    assert prepared.beta.shape == (2, 2)
    assert prepared.ind.tolist() == [[1, 1], [0, 1]]
    assert prepared.ld[0].tolist() == [[1.0, 0.25], [0.25, 1.0]]
    assert prepared.ld[1].tolist() == [[0.0, 0.0], [0.0, 1.0]]
