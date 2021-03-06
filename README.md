<div align="center">
    <h1>ARTIC tools</h1>
    <h3>a set of tools for working with the ARTIC pipeline</h3>
    <hr>
    <a href="https://github.com/will-rowe/artic-tools/actions"><img src="https://github.com/will-rowe/artic-tools/workflows/CI/badge.svg" alt="CI Status" /></a>
    <a href="https://artic-tools.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/artic-tools/badge/?version=latest" alt="Documentation Status"/></a>
    <a href="https://bioconda.github.io/recipes/artic-tools/README.html"><img src="https://anaconda.org/bioconda/artic-tools/badges/downloads.svg" alt="bioconda"></a>
    <a href="https://github.com/will-rowe/artic-tools/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
    <a href="https://zenodo.org/badge/latestdoi/280260660"><img src="https://zenodo.org/badge/280260660.svg" alt="DOI"></a>
</div>

---

## Overview

This collection of tools is being worked on and more utility is due to be added. For now, we have:

- download primer schemes and reference sequences
- primer scheme validation
- alignment softmasking
- vcf filtering

Read the [docs](https://artic-tools.readthedocs.io/en/latest/) for more info and checkout the [ARTIC pipeline](https://github.com/artic-network/fieldbioinformatics).

## Install

### conda

```
conda install -c bioconda artic-tools
```

### source

You need CMake, Boost, HTSlib and a C++17 compiler. To download and install `artic-tools` run:

```
git clone --recursive https://github.com/will-rowe/artic-tools.git
mkdir artic-tools/build
cd artic-tools/build
cmake ..
make -j4
make test
../bin/artic-tools -h
```
