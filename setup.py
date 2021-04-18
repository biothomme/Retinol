#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
from setuptools import setup

CURRENT = pathlib.Path(__file__).parent
README = (CURRENT / "README.md").read_text()

# main setup
setup(name="bumbleview",
      version="0.2.3",
      description=
      "Convert physical spectra to excitation potential in insect eyes",
      long_description=README,
      long_description_content_type="text/markdown",
      url="https://github.com/biothomme/Retinol.git",
      author="thomas m huber",
      author_email="thomas.huber@evobio.eu",
      license="MIT",
      keywords=["pollination", "ecology", "insect", "vision", "wavelength",
                "color triangle", "color hexagon", "trichromatic", "erg",
                "electroretinogram"],
      packages=["bumbleview"],
      include_package_data=True,
      # python_requires="3.7",
      install_requires=[
          "pandas", "matplotlib", "seaborn", "scipy", "statsmodels"]
      )
