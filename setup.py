from setuptools import Extension, setup

import pybind11

ext_modules = [
    Extension(
        "susiex_python",
        [
        "src/python/pybind_module.cpp",
            "src/api.cpp",
            "src/memory_loader.cpp",
            "src/validation.cpp",
            "src/data.cpp",
            "src/model.cpp",
        ],
        include_dirs=[pybind11.get_include(), "include"],
        extra_compile_args=["-std=c++11", "-fopenmp"],
        extra_link_args=["-fopenmp"],
        language="c++",
    )
]

setup(
    ext_modules=ext_modules,
    package_dir={"": "src"},
    packages=["susiex_cli"],
    zip_safe=False,
)
