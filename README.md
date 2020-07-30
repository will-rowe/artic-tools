# ARTIC tools

[![Build Status](https://travis-ci.org/will-rowe/artic_tools.svg?branch=master)](https://travis-ci.org/will-rowe/artic_tools)

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
