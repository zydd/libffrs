import pathlib
import setuptools
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent / "pybind11"))

import pybind11.setup_helpers

__version__ = "0.1.0"

ext_modules = [
    pybind11.setup_helpers.Pybind11Extension(
        "libffrs",
        [
            "libffrs/pyffrs.cpp",
        ],
        define_macros=[("VERSION_INFO", f'"{__version__}"')],
        include_dirs=["libffrs"],
        cxx_std=17,
        extra_compile_args=[
            "-Wextra",
            "-Wpedantic",
        ],
        # extra_link_args=[],
    ),
]

setuptools.setup(
    name="ffrs",
    version=__version__,
    author="Gabriel Machado",
    author_email="gabriel_machado@live.com",
    url="https://github.com/zydd",
    description="Fairly Fast & Flexible ECC",
    long_description="",
    ext_modules=ext_modules,
    packages=["ffrs"],
    cmdclass={"build_ext": pybind11.setup_helpers.build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)
