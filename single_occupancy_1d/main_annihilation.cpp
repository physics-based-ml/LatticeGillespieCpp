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
#include "gillespie1d.hpp"

using namespace Gillespie1D;

// Function to get the mean of a vector
double get_mean(std::vector<int> v) {
	double ctr = 0.0;
	for (auto i: v) { ctr += i; };
	return ctr / v.size();
};

// Function to get the std dev of a vector
double get_std(std::vector<int> v) {
	double mean = get_mean(v);
	double ctr = 0.0;
	for (auto i: v) { ctr += pow(i - mean,2); };
	return sqrt(ctr / (v.size()-1));
};

int main() {

	// Seed random no
	srand (time(NULL));
	// srand (2);

	// Initital counts of species
	std::map<std::string,int> counts0;

	counts0["A"] = 100;

	// Box length
	int box_length = 1000;

	// Timestep
	double dt = 0.01;

	// Number of steps to run
	int n_steps = 101;

	// Two modes
	for (int imode=0; imode<2; imode++)
	{
		int min,max;
		if (imode==0) {min=1; max=11;}
		else if (imode==1) {min=11; max=21;};

		/****************************************
		Make a simulation!
		****************************************/

		for (int write_version=min; write_version<max; write_version++)
		{

			// Make a simulation
			Simulation sim(dt,box_length);

			/********************
			Add species
			********************/

			sim.add_species("A");

			/********************
			Add bimol rxns
			********************/

			// M1
			sim.add_bi_rxn("ann", 0.01, "A","A");

			/********************
			Populate the lattice
			********************/

			if (imode==0) {
				sim.populate_lattice_mode_1(100);
			} else if (imode==1) {
				sim.populate_lattice_mode_2(100);
			};

			/********************
			Run
			********************/

			// Run
			bool verbose = true;
			bool write_counts = true;
			bool write_nns = true;
			bool write_latt = true;
			int write_step = 1;
			sim.run(n_steps,verbose,write_counts,write_nns,write_latt,write_step,write_version);

		};
	};

	return 0;
}