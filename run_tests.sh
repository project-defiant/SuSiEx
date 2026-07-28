#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${TMPDIR:-/tmp}/susiex-tests"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$BUILD_DIR"

compile_and_run() {
    local name="$1"
    shift
    echo "Compiling ${name}..."
    g++ "$@" -o "$BUILD_DIR/$name"
    echo "Running ${name}..."
    "$BUILD_DIR/$name"
}

COMMON=("-std=c++11" "-I" "$ROOT_DIR/include")
OPENMP=("-fopenmp")
NATIVE=("$ROOT_DIR/src/data.cpp" "$ROOT_DIR/src/memory_loader.cpp" "$ROOT_DIR/src/validation.cpp")
MODEL=("$ROOT_DIR/src/model.cpp")

compile_and_run tests_memory_loader \
    "${COMMON[@]}" \
    "$ROOT_DIR/tests/cpp/tests_memory_loader.cpp" "${NATIVE[@]}"

(
    cd "$BUILD_DIR"
    compile_and_run tests_susiex_run \
        "${COMMON[@]}" "${OPENMP[@]}" \
        "$ROOT_DIR/tests/cpp/tests_susiex_run.cpp" \
        "${NATIVE[@]}" "${MODEL[@]}"
)

compile_and_run tests_validation_unit \
    "${COMMON[@]}" \
    "$ROOT_DIR/tests/cpp/tests_validation_unit.cpp" \
    "$ROOT_DIR/src/data.cpp" "$ROOT_DIR/src/validation.cpp"

compile_and_run tests_api_error_codes \
    "${COMMON[@]}" "${OPENMP[@]}" \
    "$ROOT_DIR/tests/cpp/tests_api_error_codes.cpp" "$ROOT_DIR/src/api.cpp" \
    "${NATIVE[@]}" "${MODEL[@]}"

compile_and_run tests_validation_ld_diag \
    "${COMMON[@]}" "${OPENMP[@]}" \
    "$ROOT_DIR/tests/cpp/tests_validation_ld_diag.cpp" "${NATIVE[@]}" "${MODEL[@]}"

compile_and_run tests_validation_ld_asymmetry \
    "${COMMON[@]}" "${OPENMP[@]}" \
    "$ROOT_DIR/tests/cpp/tests_validation_ld_asymmetry.cpp" "${NATIVE[@]}" "${MODEL[@]}"

echo "All native tests passed"
