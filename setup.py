#!/usr/bin/python

try:
  from setuptools import setup
  from setuptools.extension import Extension
except Exception:
  from distutils.core import setup
  from distutils.extension import Extension

from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

_extra = ['-O3', '-ffast-math']


extensions = [
  Extension(
    'ccvt',
    sources = ['./src/ccvt.pyx'],
    extra_compile_args = _extra
  )
]

setup(
  name = "ccvt",
  version = '0.0.1',
  author = '@inconvergent',
  install_requires = ['numpy>=1.8.2', 'cython>=0.20.0'],
  license = 'MIT',
  cmdclass={'build_ext' : build_ext},
  include_dirs = [numpy.get_include()],
  ext_modules = cythonize(extensions,include_path = [numpy.get_include()])
)

