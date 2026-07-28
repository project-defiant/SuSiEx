import json

import numpy as np
from typer.testing import CliRunner

from susiex_cli.cli import app


def test_cli_writes_summary_to_nested_output(valid_fit_arrays, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_dir = tmp_path / "inputs"
    output = tmp_path / "results" / "fit.json"
    input_dir.mkdir()
    paths = {}
    for name, values in valid_fit_arrays.items():
        path = input_dir / f"{name}.npy"
        np.save(path, values)
        paths[name] = path

    result = CliRunner().invoke(
        app,
        [
            "--beta",
            str(paths["beta"]),
            "--pval",
            str(paths["pval"]),
            "--ind",
            str(paths["ind"]),
            "--ld",
            str(paths["ld"]),
            "--mk-idx",
            str(paths["mk_idx"]),
            "--n-gwas",
            "10000",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.stdout
    summary = json.loads(output.read_text())
    assert summary["converged"] is True
    assert summary["n_components"] == 5
    assert summary["n_populations"] == 1
    assert summary["n_variants"] == 5
    assert len(summary["credible_set_sizes"]) == 5
    assert summary["credible_set_sizes"][0] > 0


def test_cli_reports_invalid_sample_size(valid_fit_arrays, tmp_path):
    paths = {}
    for name, values in valid_fit_arrays.items():
        path = tmp_path / f"{name}.npy"
        np.save(path, values)
        paths[name] = path

    result = CliRunner().invoke(
        app,
        [
            "--beta",
            str(paths["beta"]),
            "--pval",
            str(paths["pval"]),
            "--ind",
            str(paths["ind"]),
            "--ld",
            str(paths["ld"]),
            "--mk-idx",
            str(paths["mk_idx"]),
            "--n-gwas",
            "not-an-integer",
            "--output",
            str(tmp_path / "result.json"),
        ],
    )

    assert result.exit_code == 1
    assert not (tmp_path / "result.json").exists()
