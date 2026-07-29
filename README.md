# SuSiEx

> [!IMPORTANT]
> **Open Targets fork.** This repository is maintained by
> [Open Targets](https://www.opentargets.org/) for integration with the
> [`opentargets/nf-fine-mapping`](https://github.com/opentargets/nf-fine-mapping)
> pipeline. It restructures the native implementation and adds Python
> bindings, the pipeline-facing CLI, validated Gentropy-compatible schemas,
> status and extended-result outputs, reproducible packaging, containers, and
> CI. The scientific method and original C++ implementation originate from the
> canonical [`getian107/SuSiEx`](https://github.com/getian107/SuSiEx)
> project. This fork is not the canonical upstream distribution; please retain
> the original attribution and citation when using it.

SuSiEx is a C++ cross-ancestry fine-mapping implementation exposed through a
Python package. The numerical core consumes validated in-memory arrays; the
Python layer provides array validation, logging, and the command-line boundary.

## Pipeline interface

```bash
uv sync --extra test
uv run susiex \
  pipeline \
  --fine-mapping-locus-set fine_mapping_locus_set.parquet \
  --multi-ancestry-pairwise-ld multi_ancestry_pairwise_ld \
  --study-metadata metadata.jsonl \
  --run-id RUN_001 \
  --fine-mapping-locus-set-id LOCUS_SET_001 \
  --study-locus-output results/study_locus.parquet \
  --extended-results-output results/susiex.h5ad \
  --stats-output results/stats.json
```

The Open Targets interface consumes one Gentropy-compatible fine-mapping locus
set, one multi-ancestry pairwise LD dataset, and JSONL study metadata. A
successful converged fit writes a flat StudyLocus Parquet file, an AnnData H5AD
file containing the full numerical result, and machine-readable run
statistics. Use `uv run susiex pipeline --help` for all method parameters.

The lower-level `susiex run` command remains available for direct NumPy-array
inputs and native-core development.

## Development

```bash
make dev
make native-test
make test
make lint
```

Project layout:

```text
include/susiex/   Public C++ headers
src/              Native implementation
src/python/       pybind11 extension binding
tests/cpp/        Native regression tests
tests/python/     Python package tests
```

## Citation

Please cite the SuSiEx method:

Yuan, K., Longchamps, R. J., Pardiñas, A. F., et al. Fine-mapping across
diverse ancestries drives the discovery of putative causal variants underlying
human complex traits and diseases. *Nature Genetics* **56**, 1841–1850 (2024).
https://doi.org/10.1038/s41588-024-01870-z
