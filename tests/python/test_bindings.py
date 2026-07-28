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
    assert res["converged"] is True
    assert res["alpha"].shape == (5, 5)
    assert res["pip"].shape == (5,)
    assert res["mu"].shape == (1, 5, 5)
    assert res["lbf"].shape == (5,)
    assert len(res["cs"]) == 5
    assert all(0 <= index < 5 for values in res["cs"] for index in values)
    assert res["pip"][0] >= max(res["pip"][1:])
