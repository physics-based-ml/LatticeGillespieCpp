# Gillespie on a Lattice in the single occupancy limit

Both 1D and 3D versions

## 3D

Make a shared library and link main against it as:
```
g++ -std=c++14 -O3 -c -fpic gillespie3d.cpp reactions.cpp species.cpp lattice3d.cpp
g++ -std=c++14 -O3 -shared -o libgillespie3d.so gillespie3d.o reactions.o species.o lattice3d.o
g++ -std=c++14 -O3 -L./ -o main.o main.cpp -lgillespie3d
```

## 1D

Make a shared library and link main against it as:
```
g++ -std=c++14 -O3 -c -fpic gillespie1d.cpp reactions.cpp species.cpp lattice1d.cpp
g++ -std=c++14 -O3 -shared -o libgillespie1d.so gillespie1d.o reactions.o species.o lattice1d.o
g++ -std=c++14 -O3 -L./ -o main.o main.cpp -lgillespie1d
```