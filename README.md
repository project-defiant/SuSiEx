# SuSiEx

SuSiEx is a C++ cross-ancestry fine-mapping implementation exposed through a
Python package. The numerical core consumes validated in-memory arrays; the
Python layer provides array validation, logging, and the command-line boundary.

## Getting started

```bash
uv sync --extra test
uv run susiex \
  --beta beta.npy \
  --pval pval.npy \
  --ind ind.npy \
  --ld ld.npy \
  --mk-idx mk_idx.npy \
  --n-gwas 50000,50000 \
  --output result.json
```

The current Python boundary intentionally uses NumPy arrays. Parquet inputs and
StudyLocus-compatible outputs will be added once the SuSiEx-specific schema and
allele/population semantics are finalized.

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
