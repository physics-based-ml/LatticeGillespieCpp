// Include header
#include "gillespie3d.hpp"

// iostream
#ifndef IOSTREAM_h
#define IOSTREAM_h
#include <iostream>
#include <iomanip>
#include <sstream>
#endif

// algorithm for sort
#ifndef ALGORITHM_h
#define ALGORITHM_h
#include <algorithm>
#endif

// Random numbers
#ifndef RAND_h
#define RAND_h
#include <stdlib.h>     /* srand, rand */
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

/************************************
* Namespace for Gillespie3D
************************************/

namespace Gillespie3D {

	/****************************************
	General functions
	****************************************/

	// Write vector to file
	void write_vector_to_file(std::string fname, std::vector<int> v)
	{
		std::ofstream f;
		f.open (fname);
		for (auto i : v) {
			f << i << "\n";
		};
		f.close();
	};

	// Print mvec
	std::ostream& operator<<(std::ostream& os, const std::vector<Species>& vs)
	{
		for (auto v : vs) { os << v.name << ": " << v.count << " "; };
		return os;
	};

	/****************************************
	Main simulation
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	// Constructors
	Simulation::Simulation(double dt, int box_length) : _lattice(box_length)
	{
		this->_uni_next = nullptr;
		this->_t = 0.0;
		this->_dt = dt;
	};

	// Destructor
	Simulation::~Simulation() {};

	/********************
	Add a species
	********************/

	void Simulation::add_species(std::string name, bool conserved) {
		// Add
		this->_species.push_back(Species(name,conserved));
		// Go through all reactions, add these
		for (auto i=0; i<this->_uni_rxns.size(); i++) {
			this->_species.back().add_rxn(&(this->_uni_rxns[i]));
		};
		for (auto i=0; i<this->_bi_rxns.size(); i++) {
			this->_species.back().add_rxn(&(this->_bi_rxns[i]));
		};
	};

	/********************
	 Add a reaction
	********************/

	// Unimolecular rxn
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r) {
		// Find the species
		Species *sr = _find_species(r);
		// Make the reaction
		this->_uni_rxns.push_back(UniReaction(name,kr,sr));
		// Add reaction to all species
		for (auto species = this->_species.begin(); species != this->_species.end(); species++) {
			species->add_rxn(&(this->_uni_rxns.back()));
		};
	};
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r, std::string p) {
		// Find the species
		Species *sr = _find_species(r);
		Species *sp = _find_species(p);
		// Make the reaction
		this->_uni_rxns.push_back(UniReaction(name,kr,sr,sp));
		// Add reaction to all species
		for (auto species = this->_species.begin(); species != this->_species.end(); species++) {
			species->add_rxn(&(this->_uni_rxns.back()));
		};
	};
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r, std::string p1, std::string p2) {
		// Find the species
		Species *sr = _find_species(r);
		Species *sp1 = _find_species(p1);
		Species *sp2 = _find_species(p2);
		// Make the reaction
		this->_uni_rxns.push_back(UniReaction(name,kr,sr,sp1,sp2));
		// Add reaction to all species
		for (auto species = this->_species.begin(); species != this->_species.end(); species++) {
			species->add_rxn(&(this->_uni_rxns.back()));
		};
	};

	// Bimolecular rxn
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2) {
		// Find the species
		Species *sr1 = _find_species(r1);
		Species *sr2 = _find_species(r2);
		// Make the reaction
		this->_bi_rxns.push_back(BiReaction(name,prob,sr1,sr2));
		// Add reaction to all species
		for (auto species = this->_species.begin(); species != this->_species.end(); species++) {
			species->add_rxn(&(this->_bi_rxns.back()));
		};
	};
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p) {
		// Find the species
		Species *sr1 = _find_species(r1);
		Species *sr2 = _find_species(r2);
		Species *sp = _find_species(p);
		// Make the reaction
		this->_bi_rxns.push_back(BiReaction(name,prob,sr1,sr2,sp));
		// Add reaction to all species
		for (auto species = this->_species.begin(); species != this->_species.end(); species++) {
			species->add_rxn(&(this->_bi_rxns.back()));
		};
	};
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2) {
		// Find the species
		Species *sr1 = _find_species(r1);
		Species *sr2 = _find_species(r2);
		Species *sp1 = _find_species(p1);
		Species *sp2 = _find_species(p2);
		// Make the reaction
		this->_bi_rxns.push_back(BiReaction(name,prob,sr1,sr2,sp1,sp2));
		// Add reaction to all species
		for (auto species = this->_species.begin(); species != this->_species.end(); species++) {
			species->add_rxn(&(this->_bi_rxns.back()));
		};
	};

	/********************
	Populate lattice
	********************/

	void Simulation::populate_lattice(std::map<std::string,int> counts) {
		std::map<Species*,int> mcounts;
		for (auto c: counts) {
			mcounts[_find_species(c.first)] = c.second;
		};
		this->_lattice.populate_lattice(mcounts);
	};

	/********************
	Run simulation
	********************/

	void Simulation::run(int n_timesteps, bool verbose, bool write_statistics)
	{
		// Go over all timesteps
		double t_next;
		for (int i_step=0; i_step < n_timesteps; i_step++) 
		{
			// The next time
			t_next = this->_t + this->_dt;

			// Print if needed
			if (verbose) {
				if (int(this->_t) % 20 == 0) {
					std::cout << "Timestep: " << this->_t << " " << this->_species << std::endl;
				};
			};

			// Do we need to schedule unimolecular reactions?
			// (Either start of sim, or there were none available before)
			if (this->_uni_next == nullptr) {
				// Schedule unimolecular reactions
				_schedule_uni();
			};

			// Do uni reactions
			while (this->_uni_next != nullptr && this->_t_uni_next < t_next) {
				// Do it
				this->_lattice.do_uni_rxn(this->_uni_next);
				// Advance time
				this->_t = this->_t_uni_next;
				// Schedule
				_schedule_uni();
			};

			// Diffuse and do bimol reactions
			this->_lattice.diffuse_mols();

			// It is now the next time
			this->_t = t_next;

		};
	};

	/********************
	Write lattice
	********************/

	void Simulation::write_lattice(int index)
	{
		std::stringstream fname;
		fname << "latt_" << std::setfill('0') << std::setw(4) << index << ".txt";
		this->_lattice.write_to_file(fname.str());
		fname.str("");
	};

	/****************************************
	Main simulation - PRIVATE
	****************************************/

	/********************
	Find a species by name
	********************/

	Species* Simulation::_find_species(std::string name) {
		// Go through the species
		for (auto i=0; i< this->_species.size(); i++) 
		{
			if (this->_species[i].name == name) 
			{ 
				return &(this->_species[i]); 
			};
		};

		return nullptr;
	};

	/********************
	Schedule the next uni reaction
	********************/

	void Simulation::_schedule_uni() {
		double props_cum = 0.0;
		std::vector<double> props;
		props.push_back(0.0);

		// Go through all possible reagants, calculate propensities
		for (auto u: this->_uni_rxns) 
		{
			props_cum += u.r->count * u.kr;
			props.push_back(props_cum);
		};

		// Check that at least one reaction is possible
		if (!(props_cum > 0)) {
			this->_t_uni_next = this->_t;
			this->_uni_next = nullptr; // Based on nullptr, will check again later
			return;
		};

		// Choose a reaction
		double r = randD(0.0,props_cum);
		int n = 0;
		while (!(props[n] < r && r < props[n+1])) {
			n++;
		};
		this->_uni_next = &(this->_uni_rxns[n]);

		// Time of next reaction
		this->_t_uni_next = this->_t + log(1.0/randD(0.0,1.0))/props_cum;
	};
};