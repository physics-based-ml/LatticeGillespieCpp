# Gillespie on a Lattice in the single occupancy limit

## Build

Use the makefile
```
make
```
Then link your code against the library using
```
g++ -std=c++14 -O3 -llatticegillespie -o main.o main.cpp
```
Be sure to put it into your `DYLD_LIBRARY_PATH`:
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/absolute/path/to/lib/folder
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/oernst/github_public_repos/LatticeGillespieCpp/lib
```
else use the `-L` flag when linking. (`-L` needed at link time; `DYLD_LIBRARY_PATH` needed at runtime).

## Features

Unimolecular and bimolecular reactions on 1D, 2D, 3D lattices.

![example](/figures/example.png)