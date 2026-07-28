from setuptools import setup, Extension
import sys
import setuptools

try:
    import pybind11
    include_dirs = [pybind11.get_include()]
except Exception:
    include_dirs = []

ext_modules = [
    Extension(
        'susiex_python',
        ['src/pybind_module.cpp', 'src/api.cpp', 'src/memory_loader.cpp', 'src/validation.cpp', 'src/data.cpp', 'src/model.cpp'],
        include_dirs=include_dirs + ['src'],
        extra_compile_args=['-std=c++11', '-fopenmp'],
        extra_link_args=['-fopenmp'],
        language='c++'
    )
]

setup(
    name='susiex_python',
    version='0.0.1',
    author='project-defiant',
    description='Python bindings for SuSiEx (minimal)',
    ext_modules=ext_modules,
    zip_safe=False,
)
