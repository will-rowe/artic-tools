<div align="center">
    <h1>ARTIC tools</h1>
    <h3>a set of tools for working with the ARTIC pipeline</h3>
    <hr>
    <a href="https://travis-ci.org/will-rowe/artic-tools"><img src="https://travis-ci.org/will-rowe/artic-tools.svg?branch=master" alt="travis"></a>
    <a href='https://artic-tools.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/artic-tools/badge/?version=latest' alt='Documentation Status'/></a>
    <a href="https://bioconda.github.io/recipes/artic-tools/README.html"><img src="https://anaconda.org/bioconda/artic-tools/badges/downloads.svg" alt="bioconda"></a>
    <a href="https://github.com/will-rowe/artic-tools/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
</div>

---

## Overview

Work in progress...

At the moment it can serve as a drop in replacement for the `align_trim` functionality of the pipeline, and can also record a few extra stats and allow some extra alignment filtering. See the [docs](https://artic-tools.readthedocs.io/en/latest/softmasker/) for more info on softmasking or checkout the [ARTIC pipeline](https://github.com/artic-network/fieldbioinformatics).

## Install

`artic-tools` requires CMake and a C++17 compiler.

To download and install `artic-tools` run:

```
git clone --recursive https://github.com/will-rowe/artic-tools.git
mkdir artic-tools/build
cd artic-tools/build
cmake ..
make -j4
make test
../bin/artic-tools -h
```
