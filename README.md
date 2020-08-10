<div align="center">
    <h1>ARTIC tools</h1>
    <h3>a bioinformatics pipeline for working with virus sequencing data sequenced with nanopore</h3>
    <hr>
    <a href="https://travis-ci.org/will-rowe/artic_tools"><img src="https://travis-ci.org/will-rowe/artic_tools.svg?branch=master" alt="travis"></a>
    <a href='https://artic-tools.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/artic-tools/badge/?version=latest' alt='Documentation Status'/></a>
    <a href="https://bioconda.github.io/recipes/artic_tools/README.html"><img src="https://anaconda.org/bioconda/artic/badges/downloads.svg" alt="bioconda"></a>
    <a href="https://github.com/will-rowe/artic_tools/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
</div>

---

## Overview

Work in progress...

This is/will be a set of tools to use in the ARTIC pipeline.

At the moment it can serve as a drop in replacement for the `align_trim` functionality of the pipeline, and can also record a few extra stats and allow some extra alignment filtering.

## Install

`artic_tools` requires CMake and a C++17 compiler.

To download and install `artic_tools` run:

```
git clone --recursive https://github.com/will-rowe/artic_tools.git
mkdir artic_tools/build
cd artic_tools/build
cmake ..
make -j4
make test
../bin/artic_tools -h
```
