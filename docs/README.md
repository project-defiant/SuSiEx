# SuSiEx documentation

SuSiEx provides a Python interface over the C++ cross-ancestry fine-mapping
implementation. The native core operates on in-memory arrays; the Python layer
validates inputs, calls the extension, and exposes the command-line interface.

## Install locally

Install the locked development environment with [uv](https://docs.astral.sh/uv/):

```bash
uv sync --extra test
uv run susiex --help
```

The package requires Python 3.11–3.13, a C++ compiler, and OpenMP. The compiler
and OpenMP linker flags are used while building the native extension.

## Run the CLI

The current boundary accepts NumPy arrays:

```bash
uv run susiex \
  --beta beta.npy \
  --pval pval.npy \
  --ind ind.npy \
  --ld ld.npy \
  --mk-idx mk_idx.npy \
  --n-gwas 50000,50000 \
  --output result.json
```

The arrays have the following shapes:

| Argument | Shape | Meaning |
| --- | --- | --- |
| `beta` | `(populations, variants)` | Effect-size statistics |
| `pval` | `(populations, variants)` | P-values |
| `ind` | `(populations, variants)` | Variant-presence indicators |
| `ld` | `(populations, variants, variants)` | LD correlation matrices |
| `mk-idx` | `(variants,)` | Native marker indices |

`n-gwas` supplies one positive GWAS sample size per population. The LD
matrices must have unit diagonal and be symmetric; invalid matrices fail before
model fitting.

## Run in Docker

Build the image from the repository root:

```bash
docker build --tag susiex:local .
docker run --rm susiex:local --help
```

To run a fit, mount a directory containing the NumPy inputs and output:

```bash
docker run --rm \
  --volume "$PWD/data:/data" \
  susiex:local \
  --beta /data/beta.npy \
  --pval /data/pval.npy \
  --ind /data/ind.npy \
  --ld /data/ld.npy \
  --mk-idx /data/mk_idx.npy \
  --n-gwas 50000,50000 \
  --output /data/result.json
```

The image entrypoint is `susiex`, so CLI options can be passed directly after
the image name.

## Development and tests

```bash
make dev
make native-test
make test
make lint
```

The native test runner compiles into a temporary directory. It covers memory
loading, model execution, API error codes, LD diagonal validation, LD symmetry,
and NaN handling. Python tests cover the pybind smoke path and Python-side
shape/sample-size validation.

## Project layout

```text
include/susiex/   Public C++ headers
src/              Native implementation
src/python/       pybind11 extension
src/susiex_cli/   Python adapter and CLI
tests/cpp/        Native regression tests
tests/python/     Python tests
```

Parquet and StudyLocus-compatible adapters are intentionally not documented as
available yet; they require a finalized SuSiEx-specific schema contract.
