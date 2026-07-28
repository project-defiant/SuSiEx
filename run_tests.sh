#!/bin/bash
set -euo pipefail

echo "Compiling memory loader test..."
g++ -std=c++11 -I src src/tests_memory_loader.cpp src/memory_loader.cpp src/validation.cpp src/data.cpp -o tests_memory_loader

if [ $? -ne 0 ]; then
    echo "Compile failed"; exit 1
fi

echo "Running memory loader test..."
./tests_memory_loader

if [ $? -ne 0 ]; then
    echo "memory loader test failed"; exit 1
fi

# compile and run susiex run test

echo "Compiling susiex run test..."
g++ -std=c++11 -fopenmp -I src src/tests_susiex_run.cpp src/memory_loader.cpp src/validation.cpp src/data.cpp src/model.cpp -o tests_susiex_run

if [ $? -ne 0 ]; then
    echo "Compile susiex run failed"; exit 1
fi

echo "Running susiex run test..."
./tests_susiex_run

# compile and run validation unit test (direct validator tests)

echo "Compiling validation unit test..."
g++ -std=c++11 -I src src/tests_validation_unit.cpp src/validation.cpp src/data.cpp -o tests_validation_unit

if [ $? -ne 0 ]; then
    echo "Compile validation unit test failed"; exit 1
fi

echo "Running validation unit test..."
./tests_validation_unit

# compile and run validation tests

echo "Compiling validation diag test..."
g++ -std=c++11 -fopenmp -I src src/tests_validation_ld_diag.cpp src/memory_loader.cpp src/validation.cpp src/data.cpp src/model.cpp -o tests_validation_ld_diag

if [ $? -ne 0 ]; then
    echo "Compile validation diag test failed"; exit 1
fi

echo "Running validation diag test (expected to fail)..."
./tests_validation_ld_diag || true

echo "Validation diag test finished (expected failure)"

# compile and run validation asymmetry test

echo "Compiling validation asymmetry test..."
g++ -std=c++11 -fopenmp -I src src/tests_validation_ld_asymmetry.cpp src/memory_loader.cpp src/validation.cpp src/data.cpp src/model.cpp -o tests_validation_ld_asymmetry

if [ $? -ne 0 ]; then
    echo "Compile validation asymmetry test failed"; exit 1
fi

echo "Running validation asymmetry test (expected to fail)..."
./tests_validation_ld_asymmetry || true

echo "Validation asymmetry test finished (expected failure)"

echo "All tests finished"
