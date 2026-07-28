import numpy as np

from susiex_cli.outputs import study_locus_records


def test_study_locus_records_preserve_native_component_alpha(valid_fit_arrays) -> None:
    result = {
        "alpha": np.array([[0.8, 0.2]]),
        "lbf": np.array([2.0]),
        "cs": [np.array([0, 1], dtype=np.int32)],
    }
    variants = [
        {"variantId": "1_100_A_G", "chromosome": "1", "position": 100},
        {"variantId": "1_110_C_T", "chromosome": "1", "position": 110},
    ]

    output = study_locus_records(
        result,
        variants,
        run_id="RUN_A",
        fine_mapping_locus_set_id="LOCUS_A",
        study_id="STUDY_A",
        chromosome="1",
        locus_start=90,
        locus_end=120,
    )

    assert len(output) == 1
    assert output[0]["variantId"] == "1_100_A_G"
    assert output[0]["locus"][0]["posteriorProbability"] == 0.8
    assert output[0]["locus"][1]["posteriorProbability"] == 0.2
    assert output[0]["locus"][0]["is95CredibleSet"] is False
