import pytest
import susiex_python


def test_susiex_pyfit_validation_error(valid_fit_arrays):
    valid_fit_arrays["ld"][0, 0, 1] = 0.1
    valid_fit_arrays["ld"][0, 1, 0] = 0.2
    with pytest.raises(RuntimeError):
        susiex_python.susiex_pyfit(*valid_fit_arrays.values())


def test_susiex_pyfit_result_contract(valid_fit_arrays):
    res = susiex_python.susiex_pyfit(*valid_fit_arrays.values())

    assert isinstance(res, dict)
    assert "alpha" in res and "pip" in res and "mu" in res and "cs" in res
    assert "component_purity" in res
    assert "component_min_p_values" in res
    assert "component_min_p_values_by_population" in res
    assert "component_is_filtered" in res
    assert "population_causal_prob" in res
    assert "logbf_by_population" in res
    assert res["converged"] is True
    assert res["alpha"].shape == (5, 5)
    assert res["pip"].shape == (5,)
    assert res["mu"].shape == (1, 5, 5)
    assert res["lbf"].shape == (5,)
    assert res["component_purity"].shape == (5,)
    assert res["component_min_p_values"].shape == (5,)
    assert res["component_min_p_values_by_population"].shape == (5, 1)
    assert res["component_is_filtered"].shape == (5,)
    assert res["population_causal_prob"].shape == (5, 1)
    assert res["logbf_by_population"].shape == (5, 1, 5)
    assert len(res["cs"]) == 5
    assert all(0 <= index < 5 for values in res["cs"] for index in values)
    assert res["pip"][0] >= max(res["pip"][1:])
    active = res["component_is_filtered"] == 0
    assert (res["component_purity"][active] > 0).all()
    assert (res["component_min_p_values"] >= 0).all()
    assert (res["population_causal_prob"] >= 0).all()
    assert (res["population_causal_prob"] <= 1).all()


def test_susiex_pyfit_reports_non_convergence(valid_fit_arrays):
    res = susiex_python.susiex_pyfit(
        *valid_fit_arrays.values(),
        n_gwas=[10000],
        max_iter=1,
    )

    assert res["converged"] is False
