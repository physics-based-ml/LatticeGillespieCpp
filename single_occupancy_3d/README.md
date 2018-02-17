# 3D lattice in single site occupancy limit

Make a shared library and link main against it as:

```
g++ -std=c++14 -O3 -c -fpic gillespie3d.cpp reactions.cpp species.cpp lattice.cpp
g++ -std=c++14 -O3 -shared -o libgillespie3d.so gillespie3d.o reactions.o species.o lattice.o
g++ -std=c++14 -O3 -L./ -o main_pre.o main_pre.cpp -lgillespie3d
```