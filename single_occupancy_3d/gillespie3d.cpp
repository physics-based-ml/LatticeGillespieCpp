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
	std::ostream& operator<<(std::ostream& os, const std::list<Species>& vs)
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
		this->_t_step = 0;
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
		auto itu = this->_uni_rxns.begin();
		while (itu != this->_uni_rxns.end()) {
			this->_species.back().add_rxn(&(*itu));
			itu++;
		};
		auto itb = this->_bi_rxns.begin();
		while (itb != this->_bi_rxns.end()) {
			this->_species.back().add_rxn(&(*itb));
			itb++;
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
	    // Go through all species
		Species *s;
		for (auto c: counts) 
		{
			s = _find_species(c.first);
	    	for (auto i=0; i<c.second; i++) 
	    	{
	    		// Make
		    	this->_lattice.make_mol_random(s);	
	    	};
		};
	};
	
	/********************
	Do a uni reaction
	********************/

	void Simulation::do_uni_rxn(UniReaction *rxn) {

		// Declarations
		std::pair<bool,SiteIt> get_it;
		SiteIt sit;
		Site s,snbr;
		std::pair<bool,Site> free_pair;

		// Try to do the reaction at several sites
		// Failure can arise if there is not enough room for the products
		int ctr = 0;
		while (ctr < 20) {

			// Grab a random site
			get_it = this->_lattice.get_mol_random_it(rxn->r);
			if (!(get_it.first)) {
				// No sites with this species exist; stop
				return;
			};
			sit = get_it.second;
			s = Site(sit);

			// Check if there is room for the products
			if (rxn->p.size() == 2) {
				free_pair = this->_lattice.get_free_neighbor_random(sit);
				if (!(free_pair.first)) {
					// Not enough room for the two products; try again with a different random site
					ctr += 1;
					continue;
				} else {
					snbr = free_pair.second;
				};
			};

			// Remove the reactant
			this->_lattice.erase_mol_it(sit);

			// Conserve reactants
			if (rxn->r->conserved) {
				this->_lattice.make_mol_random(rxn->r);
			};

			// Place products, if needed at the neighbor site
			if (rxn->p.size() == 1) {
				this->_lattice.make_mol(s, rxn->p[0]);
			} else if (rxn->p.size() == 2) {
				this->_lattice.make_mol(s, rxn->p[0]);
				this->_lattice.make_mol(snbr, rxn->p[1]);
			};

			// Conserve products
			for (auto p: rxn->p)
			{
				if (p->conserved) {
					this->_lattice.erase_mol_random(p);
				};
			};

			// Sucess
			return;
		};
	};

	/********************
	Diffuse all the mols and do bimol reactions
	********************/

	void Simulation::diffuse_mols() 
	{
		// Copy the old map
		Lattice todo = this->_lattice;

		// Clear the current lattice
		this->_lattice.clear();

		// Declarations needed
		std::pair<bool,SiteIt> old_pair;
		SiteIt sOldIt,sNewIt;
		Site sNew,sOld;
		Mol *mOld,*mColl;
		std::pair<Site,std::pair<bool,SiteIt>> nbr_pair;
		bool occ_todo, occ_done;
		std::pair<bool,SiteIt> coll_done_pair;
		std::pair<bool,BiReaction*> rxn_pair;
		BiReaction *rxn; 

		// Go over all mols to move
		while (todo.size() > 0) {

			if (DIAG_DIFFUSE) { std::cout << std::endl; };

			// Reset
			occ_todo = false;
			occ_done = false;

			// Grab some element
			old_pair = todo.get_mol_random_it();

			sOldIt = old_pair.second;
			mOld = &(sOldIt.it_2->second);

			if (DIAG_DIFFUSE) { std::cout << "diffuse_mols: got element..." << std::flush; };

			// Move
			nbr_pair = todo.get_neighbor_random(sOldIt);
			sNew = nbr_pair.first;

			if (DIAG_DIFFUSE) { std::cout << "got neighbor..." << std::flush; };

			// Check if occupied in todo lattice
			occ_todo = nbr_pair.second.first;
			if (occ_todo) {
				// Yes; its collided with something in the todo pile
				sNewIt = nbr_pair.second.second;
				mColl = &(sNewIt.it_2->second);
			} else {
				// No; check if it's collided with something in the done pile
				coll_done_pair = this->_lattice.get_mol_it(sNew);
				occ_done = coll_done_pair.first;
				if (occ_done) {
					// Yes; its collided with something in the done pile
					sNewIt = coll_done_pair.second;
					mColl = &(sNewIt.it_2->second);
				};
			};

			if (DIAG_DIFFUSE) { std::cout << "check colls..." << std::flush; };

			// If unoccupied, just commit the move (diffuse)
			if (!occ_todo && !occ_done) {
				// Move
				this->_lattice.make_mol(sNew, mOld->sp);
				// Continue
				todo.erase_mol_it(sOldIt);
				continue;
			};

			// Collision; check if reaction occurs
			rxn_pair = mOld->check_bi_rxns_mol(mColl);

			if (DIAG_DIFFUSE) { std::cout << "checked rxn..." << std::flush; };

			// No reaction occurred?
			if (!(rxn_pair.first)) {
				// No reaction & don't move
				this->_lattice.make_mol(Site(sOldIt), mOld->sp);
				// Continue
				todo.erase_mol_it(sOldIt);
				continue;
			};

			// Reaction occurred
			rxn = rxn_pair.second;
			rxn->count++;

			if (DIAG_DIFFUSE) { std::cout << "rxn occurred..." << std::flush; };

			// Remove the reactants
			sOld = Site(sOldIt); // grab it before it's erased
			todo.erase_mol_it(sOldIt);
			if (occ_todo) { 
				// Note: the old iterator sNewIt has been invalidated :(
				todo.erase_mol(sNew);
			} else if (occ_done) { 
				this->_lattice.erase_mol(sNewIt); 
			};

			if (DIAG_DIFFUSE) { std::cout << "removed r..." << std::flush; };

			// Conserve reactants
			if (rxn->r1->conserved) {
				this->_lattice.make_mol_random(rxn->r1);
			};
			if (rxn->r2->conserved) {
				this->_lattice.make_mol_random(rxn->r2);
			};

			if (DIAG_DIFFUSE) { std::cout << "conserved r..." << std::flush; };

			// Place products, if needed at the old site
			if (rxn->p.size() == 1) {
				this->_lattice.make_mol(sNew, rxn->p[0]);
			} else if (rxn->p.size() == 2) {
				this->_lattice.make_mol(sNew, rxn->p[0]);
				this->_lattice.make_mol(sOld, rxn->p[1]);
			};

			if (DIAG_DIFFUSE) { std::cout << "added p..." << std::flush; };

			// Conserve products
			for (auto p: rxn->p)
			{
				if (p->conserved) {
					this->_lattice.erase_mol_random(p);
				};
			};

			if (DIAG_DIFFUSE) { std::cout << "conserved p..." << std::flush; };

			// Finished
		};
	};


	/********************
	Run simulation
	********************/

	void Simulation::run(int n_timesteps, bool verbose, bool write)
	{
		// Clear data in files if writing
		std::ofstream ofs;
		std::stringstream fname;
		if (write) {
			for (auto s: this->_species) {
				fname << "species_" << s.name << ".txt";
				ofs.open(fname.str(), std::ofstream::out | std::ofstream::trunc);
				ofs.close();
				fname.str("");
			};
		};

		// Go over all timesteps
		double t_next;
		for (int i_step=0; i_step < n_timesteps; i_step++) 
		{
			// The next time
			t_next = this->_t + this->_dt;

			// Print/Write if needed
			if (this->_t_step % 20 == 0) {
				if (verbose) {
					std::cout << "Time: " << this->_t << " " << this->_species << std::endl;
				};
				if (write) {
					// Write species counts to file
					for (auto s: this->_species) {
						fname << "species_" << s.name << ".txt";
						ofs.open(fname.str(), std::ofstream::app);
						ofs << this->_t << " " << s.count << "\n";
						ofs.close();
						fname.str("");
					};
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
				do_uni_rxn(this->_uni_next);
				// Advance time
				this->_t = this->_t_uni_next;
				// Schedule
				_schedule_uni();
			};

			// Diffuse and do bimol reactions
			diffuse_mols();

			// It is now the next time
			this->_t = t_next;
			this->_t_step++;

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
		auto it = this->_species.begin();
		while (it != this->_species.end()) {
			if (it->name == name) {
				return &(*it);
			};
			it++;
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
		auto it = this->_uni_rxns.begin();
		int n = 0;
		while (!(props[n] < r && r < props[n+1])) {
			n++;
			it++;
		};
		this->_uni_next = &(*it);

		// Time of next reaction
		this->_t_uni_next = this->_t + log(1.0/randD(0.0,1.0))/props_cum;
	};
};