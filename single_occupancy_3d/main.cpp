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
	counts0["B"] = 100;
	counts0["C"] = 100;

	/****************************************
	Semi-attractor dynamics
	****************************************/
	
	/********************
	Version 1
	********************/

	/*
	// Version	
	int write_version_no = 1;

	double kp1=30.0, km1=0.25, kp2=1.0, kp3=10.0, kp4=0.6, kp5=16.5, km5=1.0;
	// Box length
	int box_length = 20;

	// Timestep
	double dt = 0.001;

	// Number of steps to run
	int n_steps = 10001;

	// Reaction probs
	double m = 0.15;
	double probm1 = m*km1;
	double probp2 = m*kp2;
	double probp4 = m*kp4;
	double probm5 = m*km5;
	*/

	/********************
	Version 2
	********************/

	/*
	// Version	
	int write_version_no = 2;

	double kp1=30.0, km1=0.25, kp2=1.0, kp3=10.0, kp4=0.6, kp5=16.5, km5=2.0;
	// Box length
	int box_length = 20;

	// Timestep
	double dt = 0.0001;

	// Number of steps to run
	int n_steps = 10001;

	// Reaction probs
	double m = 0.02;
	double probm1 = m*km1;
	double probp2 = m*kp2;
	double probp4 = m*kp4;
	double probm5 = m*km5;
	*/

	/********************
	Version 3
	********************/

	/*
	// Version	
	int write_version_no = 3;

 	double kp1=30.0, km1=0.25, kp2=1.0, kp3=10.0, kp4=0.6, kp5=16.5, km5=0.6;
	// Box length
	int box_length = 20;

	// Timestep
	double dt = 0.01;

	// Number of steps to run
	int n_steps = 1001;

	// Reaction probs
	double m = 0.5;
	double probm1 = m*km1;
	double probp2 = m*kp2;
	double probp4 = m*kp4;
	double probm5 = m*km5;
	*/

	/********************
	Version 4
	********************/

	/*
	// Version	
	int write_version_no = 4;

 	double kp1=30.0, km1=0.25, kp2=1.0, kp3=10.0, kp4=0.6, kp5=16.5, km5=0.9;
	// Box length
	int box_length = 10;

	// Timestep
	double dt = 0.01;

	// Number of steps to run
	int n_steps = 1001;

	// Reaction probs
	double m = 0.4;
	double probm1 = m*km1;
	double probp2 = m*kp2;
	double probp4 = m*kp4;
	double probm5 = m*km5;
	*/

	/****************************************
	New
	****************************************/

	// Version	
	int write_version_no = 5;

	counts0["A"] = 800;
	counts0["B"] = 800;
	counts0["C"] = 800;

 	double kp1=30.0, km1=0.25, kp2=1.0, kp3=10.0, kp4=0.6, kp5=16.5, km5=0.9;
	// Box length
	int box_length = 20;

	// Timestep
	double dt = 0.01;

	// Number of steps to run
	int n_steps = 1001;

	// Reaction probs
	double m = 0.4;
	double probm1 = m*km1;
	double probp2 = m*kp2;
	double probp4 = m*kp4;
	double probm5 = m*km5;

	/****************************************
	Make a simulation!
	****************************************/

	// Make a simulation
	Simulation sim(dt,box_length);

	/********************
	Add species
	********************/

	sim.add_species("A");
	sim.add_species("B");
	sim.add_species("C");

	/********************
	Add unimol rxns
	********************/

	// P1
	sim.add_uni_rxn("rxnP1", kp1, "A","A","A");
	// P3
	sim.add_uni_rxn("rxnP3", kp3, "B");
	// P5
	sim.add_uni_rxn("rxnP5", kp5, "C","C","C");

	/********************
	Add bimol rxns
	********************/

	// M1
	sim.add_bi_rxn("rxnM1", probm1, "A","A","A");
	// P2
	sim.add_bi_rxn("rxnP2", probp2, "A","B","B","B");
	// P4
	sim.add_bi_rxn("rxnP4", probp4, "A","C");
	// M5
	sim.add_bi_rxn("rxnM5", probm5, "C","C","C");

	/********************
	Populate the lattice
	********************/

	sim.populate_lattice(counts0);

	/********************
	Run
	********************/

	// Run
	bool verbose = true;
	bool write_counts = true;
	bool write_nns = true;
	bool write_latt = true;
	int write_step = 1;
	sim.run(n_steps,verbose,write_counts,write_nns,write_latt,write_step,write_version_no);

	/********************
	Write
	********************/

	// Write the lattice to a file
	// sim.write_lattice(0);

	// Write the lattice to a file
	// sim.write_lattice(n_steps);

	return 0;
}