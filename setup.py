#!/usr/bin/env python
"""
setup.py file for SWIG example
"""
from setuptools import setup, Extension

distlink_dir = "./distlink-2.0"  # path to cpp library
distlink_module = Extension(
    name="_distlink",
    sources=["distlink_wrap.cxx", f"./{distlink_dir}/distlink.cpp"],
    include_dirs=[distlink_dir, ],
    extra_compile_args=["-O3", "-march=native", "-mfpmath=sse"]  # taken from README.txt
)

setup(
    name="distlink",
    version="2.0",
    author="Maximilian Pierzyna",
    description="""Python wrapper for distlink library by R.V. Baluev and D.V. Mikryukov""",
    ext_modules=[distlink_module],
    py_modules=["distlink"],
)
