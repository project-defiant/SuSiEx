import numpy as np

from susiex_cli.stats import stats_from_result, write_stats


def test_stats_report_native_shape_and_non_convergence(tmp_path) -> None:
    result = {
        "alpha": np.zeros((3, 4)),
        "pip": np.zeros(4),
        "mu": np.zeros((2, 4)),
        "cs": [np.array([0, 1]), np.array([2]), np.array([], dtype=np.int32)],
        "converged": False,
    }

    stats = stats_from_result(
        result, run_id="RUN_A", fine_mapping_locus_set_id="LOCUS_A"
    )
    output = tmp_path / "nested" / "stats.json"
    write_stats(output, stats)

    assert stats.status == "NON_CONVERGED"
    assert stats.nComponents == 3
    assert stats.nPopulations == 2
    assert stats.nVariants == 4
    assert stats.nCredibleSets == 3
    assert '"reason": "SuSiEx fit did not converge"' in output.read_text()
