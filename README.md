# Gillespie on a Lattice in the single occupancy limit

## Build

Use the makefile
```
make
```
Then link against the library using
```
g++ -std=c++14 -O3 -L../lib -llatticegillespie -o main.o main.cpp
```
Be sure to put it into your `DYLD_LIBRARY_PATH`:
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/absolute/path/to/lib/folder
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/oernst/github_public_repos/LatticeGillespieCpp/lib
```
(`-L` needed at link time; `DYLD_LIBRARY_PATH` needed at runtime).