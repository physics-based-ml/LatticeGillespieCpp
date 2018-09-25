# Gillespie on a Lattice in the single occupancy limit

## Build

Use the makefile
```
make
make install
```
The default location is `/usr/local/lib` and `/usr/local/include`.

## Including

Use the convenient `include <lattGillespie>` (note capital `G`).

## Linking

Example:
```
g++ -std=c++14 -O3 -llattgillespie -o main.o main.cpp
```
Be sure to put it into your `DYLD_LIBRARY_PATH`:
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/absolute/path/to/lib/folder
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/oernst/github_public_repos/LatticeGillespieCpp/lib
```
else use the `-L` flag when linking. (`-L` needed at link time; `DYLD_LIBRARY_PATH` needed at runtime).

## Namespace

The namespace is `latg`.

## Features

Unimolecular and bimolecular reactions on 1D, 2D, 3D lattices.

![example](/figures/example.png)