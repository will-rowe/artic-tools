# Installation

***

## Bioconda

```
conda install -c bioconda artic-tools
```

## Source

`artic-tools` requires CMake, HTSlib and a C++17 compiler.

```
git clone --recursive https://github.com/will-rowe/artic-tools.git
mkdir artic-tools/build
cd artic-tools/build
cmake ..
make -j4
make test
../bin/artic-tools -h
```
