#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   This file is part of the magpurify package, available at:
#   https://github.com/snayfach/MAGpurify
#
#   Magpurify is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <https://www.gnu.org/licenses/>.

from setuptools import find_packages, setup

setup(
    name="magpurify",
    version="2.0",
    packages=find_packages(),
    license="GNU General Public License v3.0",
    description="Identify and remove incorrectly binned contigs from metagenome-assembled genomes.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "biopython",
        "numpy",
        "pandas",
        "sklearn",
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["magpurify=magpurify.cli:cli"]},
    url="https://github.com/snayfach/MAGpurify",
    keywords=[
        "bioinformatics",
        "metagenomics",
        "metagenome-assembled genomes"
    ],
    author="Stephen Nayfach, Antonio Pedro Camargo",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Software Development :: Libraries",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python :: 3",
    ],
)
