// iostream
#ifndef IOSTREAM_h
#define IOSTREAM_h
#include <iostream>
#endif

// Math
#ifndef MATH_h
#define MATH_h
#include <math.h>
#endif

// Get the header
#include "gillespie3d.hpp"

using namespace Gillespie3D;

int main() {

	// Make a simulation
	Simulation sim(0.1,5);

	sim.add_species("A1",true);

	std::map<std::string,int> counts0;
	counts0["A1"] = 10;
	sim.populate_lattice(counts0);

	// sim.add_bi_rxn("rxn", 0.5, "A1","A1");

	// Write the initial lattice
	sim.write_lattice(0);
	sim.run(1,true,false);
	sim.write_lattice(1);
	sim.run(1,true,false);
	sim.write_lattice(2);

	return 0;
}