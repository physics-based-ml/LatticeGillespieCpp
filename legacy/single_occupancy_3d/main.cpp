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

	// Concentrations of species to keep constant
	double a1 = 30.0, a2 = 0.01, a3 = 0.01, a4 = 16.5, a5 = 10.0;

	// Desired quantity of A2,A3
	counts0["A2"] = 1;
	counts0["A3"] = 1;

	// Effective box volume, length
	int box_vol_eff = counts0["A2"] / a2;
	double box_len_eff = cbrt(box_vol_eff);

	// Quantities of species to keep constant
	counts0["A1"] = a1 * box_vol_eff;
	counts0["A4"] = a4 * box_vol_eff;
	counts0["A5"] = a5 * box_vol_eff;

	// X,Y,Z concentrations
	double x0 = 10.0, y0 = 80.0, z0 = 0.1;

	// X,Y,Z quantities
	counts0["X"] = x0*box_vol_eff;
	counts0["Y"] = y0*box_vol_eff;
	counts0["Z"] = z0*box_vol_eff;

	// Box length
	// Total of 14662 particles - min size is 25
	int box_length = 100;

	// Times
	double dt = 0.0000001;

	// Effective diffusion step, in terms of box length
	double box_step_eff = box_len_eff/box_length;
	double deff = box_step_eff / (dt * dt);

	// Rates
	double km1 = 0.25, km2 = 0.001, km5 = 0.5, kp1 = 1.0, kp2 = 1.0, kp3 = 1.0, kp4 = 1.0, kp5 = 1.0, km3 = 1.0, km4 = 1.0;

	// Reaction probabilities
	// Bimolecular reactions need to be scaled to treat effective rates
	// Units here are len^3 / time
	// Scale by box_step_effective^3
	std::map<std::string,double> rxn_probs;
	double aint = 1.2E-8; // Interaction area
	rxn_probs["rxnP1"] = sqrt(M_PI) * pow(box_step_eff,3) * kp1 / (8.0 * sqrt(deff) * aint);
	rxn_probs["rxnM1"] = sqrt(M_PI) * pow(box_step_eff,3) * km1 / (8.0 * sqrt(deff) * aint);
	rxn_probs["rxnP2"] = sqrt(M_PI) * pow(box_step_eff,3) * kp2 / (8.0 * sqrt(deff) * aint);
	rxn_probs["rxnM2"] = sqrt(M_PI) * pow(box_step_eff,3) * km2 / (8.0 * sqrt(deff) * aint);
	rxn_probs["rxnP3"] = sqrt(M_PI) * pow(box_step_eff,3) * kp3 / (8.0 * sqrt(deff) * aint);
	rxn_probs["rxnP4"] = sqrt(M_PI) * pow(box_step_eff,3) * kp4 / (8.0 * sqrt(deff) * aint);
	rxn_probs["rxnP5"] = sqrt(M_PI) * pow(box_step_eff,3) * kp5 / (8.0 * sqrt(deff) * aint);
	rxn_probs["rxnM5"] = sqrt(M_PI) * pow(box_step_eff,3) * km5 / (8.0 * sqrt(deff) * aint);
	std::cout << rxn_probs["rxnP3"] << std::endl;

	// Reaction probabilities
	/*
	std::map<std::string,double> rxn_probs;
	rxn_probs["rxnP1"] = 0.139785;
	rxn_probs["rxnM1"] = 0.100797;
	rxn_probs["rxnP2"] = 0.146503;
	rxn_probs["rxnM2"] = 0.000129401;
	rxn_probs["rxnP3"] = 0.123583;
	rxn_probs["rxnP4"] = 0.139535;
	rxn_probs["rxnP5"] = 0.056512;
	rxn_probs["rxnM5"] = 0.25;
	*/

	// Unimol Reactions
	// M3
	UniReaction rxnM3("rxnM3", km3, "A2","A5","Y");
	// M4
	UniReaction rxnM4("rxnM4", km4, "A3","X","Z");

	// Bimol Reactions
	// P1
	BiReaction rxnP1("rxnP1", rxn_probs["rxnP1"], "A1","X","X","X");
	// M1
	BiReaction rxnM1("rxnM1", rxn_probs["rxnM1"], "X","X","A1","X");
	// P2
	BiReaction rxnP2("rxnP2", rxn_probs["rxnP2"], "X","Y","Y","Y");
	// M2
	BiReaction rxnM2("rxnM2", rxn_probs["rxnM2"], "Y","Y","X","Y");
	// P3
	BiReaction rxnP3("rxnP3", rxn_probs["rxnP3"], "A5","Y","A2");
	// P4
	BiReaction rxnP4("rxnP4", rxn_probs["rxnP4"], "X","Z","A3");
	// P5
	BiReaction rxnP5("rxnP5", rxn_probs["rxnP5"], "A4","Z","Z","Z");
	// M5
	BiReaction rxnM5("rxnM5", rxn_probs["rxnM5"], "Z","Z","A4","Z");

	// Make a Lattice
	Lattice latt(box_length);

	// Populate the lattice randomly according to the counts
	latt.populate_lattice(counts0);

	// Make a simulation
	Simulation sim(latt,dt);

	// Populate the lattice according to some hdict, jdict
	/*
	std::map<std::string,double> h_dict;
	h_dict["A"] = 1.0;
	std::map<SpeciesPair,double> j_dict;
	j_dict[SpeciesPair("A","A")] = 0.0;
	latt.populate_lattice(h_dict,j_dict,200,false,true);
	*/
			
	// Add all the reactions
	sim.add_reaction(rxnM1);
	sim.add_reaction(rxnM2);
	sim.add_reaction(rxnM3);
	sim.add_reaction(rxnM4);
	sim.add_reaction(rxnM5);
	sim.add_reaction(rxnP1);
	sim.add_reaction(rxnP2);
	sim.add_reaction(rxnP3);
	sim.add_reaction(rxnP4);
	sim.add_reaction(rxnP5);

	// Add all the conservation laws
	sim.set_conserved("A1");
	sim.set_conserved("A2");
	sim.set_conserved("A3");
	sim.set_conserved("A4");
	sim.set_conserved("A5");

	// Number of steps to run
	int n_steps = 10001;

	// Run
	sim.run(n_steps,true,false);

	// Write the lattice to a file
	// sim.write_lattice(0);

	// Write the lattice to a file
	// sim.write_lattice(n_steps);

	return 0;
}