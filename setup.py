#!/usr/bin/python

try:
  from setuptools import setup, find_packages
except Exception:
  from distutils.core import setup, find_packages

packages = find_packages()

setup(
    name = 'ccvt',
    version = '0.0.1',
    author = '@inconvergent',
    install_requires = ['numpy>=1.8.2'],
    license = 'MIT',
    packages = packages
)

