import anndata as ad
import numpy as np

from susiex_cli.anndata_output import write_anndata
from susiex_cli.preparation import PreparedArrays


def test_write_anndata_preserves_native_axes_and_provenance(tmp_path) -> None:
    prepared = PreparedArrays(
        variant_ids=["1_100_A_G", "1_110_C_T", "1_120_G_A"],
        chromosomes=["1", "1", "1"],
        positions=[100, 110, 120],
        study_ids=["STUDY_A", "STUDY_B"],
        ancestries=["eur", "afr"],
        sample_sizes=[1000, 2000],
        beta=np.zeros((2, 3)),
        pval=np.ones((2, 3)),
        ind=np.ones((2, 3), dtype=np.uint8),
        ld=np.stack([np.eye(3), np.eye(3)]),
        mk_idx=np.arange(3),
    )
    result = {
        "alpha": np.ones((2, 3), dtype=np.float64) / 3,
        "pip": np.ones(3),
        "lbf": np.array([1.0, 2.0]),
        "component_purity": np.array([0.8, 0.7]),
        "component_min_p_values": np.array([1e-8, 1e-6]),
        "component_is_filtered": np.array([0, 1]),
        "population_causal_prob": np.zeros((2, 2)),
        "logbf_by_population": np.zeros((2, 2, 3)),
        "mu": np.zeros((2, 2, 3)),
        "cs": [np.array([0, 1]), np.array([2])],
        "converged": True,
    }
    output = tmp_path / "fit.h5ad"

    write_anndata(
        result,
        prepared,
        run_id="RUN_A",
        fine_mapping_locus_set_id="LOCUS_A",
        output=output,
    )
    observed = ad.read_h5ad(output)

    assert observed.shape == (2, 3)
    assert list(observed.var["variantId"]) == prepared.variant_ids
    assert {"mu__eur", "mu__afr"}.issubset(set(observed.layers.keys()))
    assert observed.uns["runId"] == "RUN_A"
    assert observed.varm["variantPresent"].shape == (3, 2)
