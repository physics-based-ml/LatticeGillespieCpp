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

// Writing to file
#ifndef FSTREAM_h
#define FSTREAM_h
#include <fstream>
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
	// srand (time(NULL));
	srand (2);

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

	// Bimolecular reaction probabilities
	std::map<std::string,double> rxn_probs;
	rxn_probs["rxnP1"] = 0.139785;
	rxn_probs["rxnM1"] = 0.100797;
	rxn_probs["rxnP2"] = 0.146503;
	rxn_probs["rxnM2"] = 0.000129401;
	rxn_probs["rxnP3"] = 0.123583;
	rxn_probs["rxnP4"] = 0.139535;
	rxn_probs["rxnP5"] = 0.056512;
	rxn_probs["rxnM5"] = 0.25;
	
	// Desired reaction counts (from Mathematica)
	std::map<std::string,double> desired_rxn_counts;
	desired_rxn_counts["rxnP1"] = 234.24;
	desired_rxn_counts["rxnM1"] = 15.42;
	desired_rxn_counts["rxnP2"] = 618.71;
	desired_rxn_counts["rxnM2"] = 6.36;
	desired_rxn_counts["rxnP3"] = 792.92;
	desired_rxn_counts["rxnP4"] = 0.72;
	desired_rxn_counts["rxnP5"] = 1.51;
	desired_rxn_counts["rxnM5"] = 0.01;

	// Number of optimization steps
	int n_opt_steps = 3;

	// Go over all optimization steps
	for (auto i_opt_step = 0; i_opt_step < n_opt_steps; i_opt_step++) 
	{
		std::cout << "--------------------" << std::endl;
		std::cout << "--- " << "Opt step: " << i_opt_step << " ---" << std::endl;
		std::cout << "--------------------" << std::endl;

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

		// Make a simulation
		Simulation sim(latt,dt);

		// Number of trials
		int n_trials = 50;

		// Number of steps
		int n_steps = 1001;

		// Store reaction counts
		std::map<std::string,std::vector<int>> rxn_count;
		std::map<std::string,std::vector<int>> rxn_count_st;
		std::map<std::string,std::vector<int>> coll_count;
		std::map<std::string,std::vector<int>> coll_count_st;

		for (auto i_trial = 0; i_trial < n_trials; i_trial++)
		{
			std::cout << "   > Trial: " << i_trial+1 << " / " << n_trials << std::endl;

			// Make a Lattice
			latt = Lattice(box_length);

			// Populate the lattice randomly according to the counts
			latt.populate_lattice(counts0);

			// Make a simulation
			sim = Simulation(latt,dt);
			
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

			// Run
			sim.run(n_steps,false,false);

			// Store the final reaction, coll counts
			rxn_count = sim.get_rxn_count();
			for (auto rc: rxn_count) {
				rxn_count_st[rc.first].push_back(rc.second.back());
			};
			coll_count = sim.get_coll_count();
			for (auto cc: coll_count) {
				coll_count_st[cc.first].push_back(cc.second.back());
			};
		};	

		// Update the probabilities 
		for (auto rc: rxn_count_st) {
			// Find this reaction (if bimol)
			auto it = desired_rxn_counts.find(rc.first);
			if (it != desired_rxn_counts.end()) {

				std::cout << "Rxn " << rc.first << " mean: " << get_mean(rc.second) << " +- " << get_std(rc.second) << " vs desired: " << it->second << " for n coll: " << get_mean(coll_count_st[rc.first]) << " +- " << get_std(coll_count_st[rc.first]) << std::endl;

				// Update to the desired / no coll
				std::cout << "   Prob: " << rxn_probs[rc.first] << " -> ";
				if (get_mean(coll_count_st[rc.first]) > 0.0001) {
					rxn_probs[rc.first] = it->second / get_mean(coll_count_st[rc.first]);
				};
				std::cout << rxn_probs[rc.first] << std::endl;
			};
		};
	};

	// Print the new reaction probabilities
	std::cout << "-----------------" << std::endl;
	std::cout << "-----------------" << std::endl;
	std::cout << "New reaction probabilities:" << std::endl;
	std::map<std::string,double> rxn_probs_new;
	for (auto rp: rxn_probs) {
		std::cout << rp.first << " " << rp.second << std::endl;
		rxn_probs_new[rp.first] = rp.second;
	};

	// Write to file
	std::ofstream f;
	f.open ("rxn_probs.txt");
	for (auto rpn : rxn_probs_new) {
		f << rpn.first << " " << rpn.second << "\n";
	};
	f.close();

	return 0;
}