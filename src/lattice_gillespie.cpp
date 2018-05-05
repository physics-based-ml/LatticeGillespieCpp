#include "../include/lattice_gillespie.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <math.h>
#include <fstream>

#include "lattice.hpp"

/************************************
* Namespace for LatticeGillespie
************************************/

namespace LatticeGillespie {

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
	Main simulation IMPLEMENTATION HEADER
	****************************************/

	class Simulation::Impl
	{
	private:

		// The dimensionality of the lattice
		int _dim;

		// The lattice
		Lattice *_lattice;

		// List of species
		std::list<Species> _species;

		// List of bimol rxns
		std::list<BiReaction> _bi_rxns;

		// List of unimol rxns
		std::list<UniReaction> _uni_rxns;

		// Box length
		int _box_length;

		// Current time
		double _t;
		int _t_step;

		// Timestep
		double _dt;

		// Time and value of the next unimolecular reaction
		double _t_uni_next;
		UniReaction *_uni_next; // Null indicates there is none scheduled

		/********************
		Find a species by name
		********************/

		Species* _find_species(std::string name);

		/********************
		Schedule the next uni reaction
		********************/

		void _schedule_uni();

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _reset();

	public:

		/********************
		Constructor/Destructor
		********************/

		Impl(double dt, int box_length, int dim);
		Impl(Impl&& other);
	    Impl& operator=(Impl&& other);
		~Impl();

		/********************
		Add species
		********************/

		void add_species(std::string name, bool conserved = false);

		/********************
		 Add a reaction
		********************/

		// Unimolecular rxn
		void add_uni_rxn(std::string name, double kr, std::string r);
		void add_uni_rxn(std::string name, double kr, std::string r, std::string p);
		void add_uni_rxn(std::string name, double kr, std::string r, std::string p1, std::string p2);

		// Bimolecular rxn
		void add_bi_rxn(std::string name, double prob, std::string r1, std::string r2);
		void add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p);
		void add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2);

		/********************
		Populate lattice
		********************/

		void populate_lattice(std::map<std::string,int> counts);
		void populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, int n_steps);
		void populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, int n_steps);
		void populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>> &q_dict, int n_steps);

		/********************
		Do a uni reaction
		********************/

		void do_uni_rxn(UniReaction *rxn);

		/********************
		Diffuse all the mols and do bimol reactions
		********************/

		void diffuse_mols();

		/********************
		Run simulation
		********************/

		void run(int n_timesteps, bool verbose = true, bool write_counts = false, bool write_nns = false, bool write_latt = false, int write_step = 20, int write_version_no = 0, std::string dir=".");

		/********************
		Write lattice
		********************/

		void write_lattice(int index, int write_version_no, std::string dir);
	};

	/****************************************
	Main simulation IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	// Constructors
	Simulation::Impl::Impl(double dt, int box_length, int dim)
	{
		_dim = dim;
		// Make lattice
		if (_dim == 1) {
			_lattice = new Lattice(box_length);
		} else if (_dim == 2) {
			_lattice = new Lattice(box_length,box_length);
		} else if (_dim == 3) {
			_lattice = new Lattice(box_length,box_length,box_length);
		} else {
			std::cerr << "ERROR: only lattice dimensions 1,2,3 are supported" << std::endl;
			exit(EXIT_FAILURE);
		};

		_uni_next = nullptr;
		_t = 0.0;
		_t_step = 0;
		_dt = dt;
		_box_length = box_length;
	};

	Simulation::Impl::Impl(Impl&& other) {
		_copy(other);
		other._reset();
	};
    Simulation::Impl& Simulation::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
    };

	// Destructor
	Simulation::Impl::~Impl() {
		_clean_up();
	};

	// Helpers
	void Simulation::Impl::_clean_up() {
		if (_lattice != nullptr) {
			delete _lattice;
		};
	};
	void Simulation::Impl::_copy(const Impl& other) {
		_dim = other._dim;
		_box_length = other._box_length;
		_lattice = new Lattice(*(other._lattice));
		_species = other._species;
		_bi_rxns = other._bi_rxns;
		_uni_rxns = other._uni_rxns;
		_t = other._t;
		_t_step = other._t_step;
		_dt = other._dt;
		_t_uni_next = other._t_uni_next;
		if (other._uni_next == nullptr) {
			_uni_next = nullptr;
		} else {
			// Search...
			for (auto it = _uni_rxns.begin(); it != _uni_rxns.end(); it++) {
				if (it->name == other._uni_next->name) {
					_uni_next = &(*it);
					break;
				};
			};
			// Can't find?
			if (_uni_next == nullptr) {
				std::cerr << "Error: can't find uni reaction when copying" << std::endl;
				exit(EXIT_FAILURE);
			};
		};
	};
	void Simulation::Impl::_reset() {
		_dim = 3;
		_lattice = nullptr;
		_species.clear();
		_bi_rxns.clear();
		_uni_rxns.clear();
		_box_length = 0;
		_t = 0.0;
		_t_step = 0;
		_dt = 0.0;
		_t_uni_next = 0.0;
		_uni_next = nullptr;
	};

	/********************
	Add a species
	********************/

	void Simulation::Impl::add_species(std::string name, bool conserved) {
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
	void Simulation::Impl::add_uni_rxn(std::string name, double kr, std::string r) {
		// Find the species
		Species *sr = _find_species(r);
		// Make the reaction
		this->_uni_rxns.push_back(UniReaction(name,kr,sr));
		// Add reaction to all species
		for (auto species = this->_species.begin(); species != this->_species.end(); species++) {
			species->add_rxn(&(this->_uni_rxns.back()));
		};
	};
	void Simulation::Impl::add_uni_rxn(std::string name, double kr, std::string r, std::string p) {
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
	void Simulation::Impl::add_uni_rxn(std::string name, double kr, std::string r, std::string p1, std::string p2) {
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
	void Simulation::Impl::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2) {
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
	void Simulation::Impl::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p) {
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
	void Simulation::Impl::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2) {
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

	void Simulation::Impl::populate_lattice(std::map<std::string,int> counts) {
	    // Go through all species
		Species *s;
		for (auto c: counts) 
		{
			s = _find_species(c.first);
	    	for (auto i=0; i<c.second; i++) 
	    	{
	    		// Make
		    	_lattice->make_mol_random(s);	
	    	};
		};
	};

	void Simulation::Impl::populate_lattice(std::map<std::string,double> &h_dict,std::map<std::string,std::map<std::string,double>> &j_dict, int n_steps) {
		std::map<std::string,std::map<std::string,std::map<std::string,double>>> k_dict;
		std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>> q_dict;		
		populate_lattice(h_dict,j_dict,k_dict,q_dict,n_steps);
	};

	void Simulation::Impl::populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, int n_steps) {
		std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>> q_dict;		
		populate_lattice(h_dict,j_dict,k_dict,q_dict,n_steps);
	};


	void Simulation::Impl::populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>> &q_dict, int n_steps)
	{
		// Start by populating lattice randomly

		// Random number of initial particles (min is 1, max is box vol)
		int n = randI(1, pow(_box_length,_dim));

		// Random initial counts
		int n_possible = pow(_box_length,_dim);
		std::map<std::string,int> counts0;
		for (auto hpr : h_dict) {
			counts0[hpr.first] = randI(0,n_possible);
			n_possible -= counts0[hpr.first];
			if (n_possible < 0) { n_possible = 0; };
		};

		// Random initial lattice
		populate_lattice(counts0);

		// TMP
		// Write the lattice into a temp dir
		// _lattice->write_to_file("lattice_v300/pre.txt");

		// Convert the strings into species
		std::map<Species*, double> h_dict_sp;
		std::map<Species*,std::map<Species*,double>> j_dict_sp;
		std::map<Species*,std::map<Species*,std::map<Species*,double>>> k_dict_sp;
		std::map<Species*,std::map<Species*,std::map<Species*,std::map<Species*,double>>>> q_dict_sp;
		for (auto hpr: h_dict) {
			h_dict_sp[_find_species(hpr.first)] = hpr.second;
		};
		for (auto jpr1: j_dict) {
			for (auto jpr2: jpr1.second) {
				j_dict_sp[_find_species(jpr1.first)][_find_species(jpr2.first)] = jpr2.second;
			};
		};
		for (auto kpr1: k_dict) {
			for (auto kpr2: kpr1.second) {
				for (auto kpr3: kpr2.second) {
					k_dict_sp[_find_species(kpr1.first)][_find_species(kpr2.first)][_find_species(kpr3.first)] = kpr3.second;
				};
			};
		};
		for (auto qpr1: q_dict) {
			for (auto qpr2: qpr1.second) {
				for (auto qpr3: qpr2.second) {
						for (auto qpr4: qpr3.second) {
						q_dict_sp[_find_species(qpr1.first)][_find_species(qpr2.first)][_find_species(qpr3.first)][_find_species(qpr4.first)] = qpr4.second;
					};
				};
			};
		};

		// Now sample
		_lattice->sample(h_dict_sp,j_dict_sp,k_dict_sp,q_dict_sp,n_steps);
	};

	/********************
	Do a uni reaction
	********************/

	void Simulation::Impl::do_uni_rxn(UniReaction *rxn) {

		// Declarations
		std::pair<bool,SiteIt3D> get_it;
		SiteIt3D sit;
		Site3D s,snbr;
		std::pair<bool,Site3D> free_pair;

		// Try to do the reaction at several sites
		// Failure can arise if there is not enough room for the products
		int ctr = 0;
		while (ctr < 20) {

			// Grab a random site
			get_it = _lattice->get_mol_random_it(rxn->r);
			if (!(get_it.first)) {
				// No sites with this species exist; stop
				return;
			};
			sit = get_it.second;
			s = Site3D(sit);

			// Check if there is room for the products
			if (rxn->p.size() == 2) {
				free_pair = _lattice->get_free_neighbor_random(sit);
				if (!(free_pair.first)) {
					// Not enough room for the two products; try again with a different random site
					ctr += 1;
					continue;
				} else {
					snbr = free_pair.second;
				};
			};

			// Remove the reactant
			_lattice->erase_mol_it(sit);

			// Conserve reactants
			if (rxn->r->conserved) {
				_lattice->make_mol_random(rxn->r);
			};

			// Place products, if needed at the neighbor site
			if (rxn->p.size() == 1) {
				_lattice->make_mol(s, rxn->p[0]);
			} else if (rxn->p.size() == 2) {
				_lattice->make_mol(s, rxn->p[0]);
				_lattice->make_mol(snbr, rxn->p[1]);
			};

			// Conserve products
			for (auto p: rxn->p)
			{
				if (p->conserved) {
					_lattice->erase_mol_random(p);
				};
			};

			// Sucess
			return;
		};
	};

	/********************
	Diffuse all the mols and do bimol reactions
	********************/

	void Simulation::Impl::diffuse_mols() 
	{
		// Copy the old map
		Lattice *todo = new Lattice(*_lattice);

		// Clear the current lattice
		_lattice->clear();

		// Declarations needed
		std::pair<bool,SiteIt3D> old_pair;
		SiteIt3D sOldIt,sNewIt;
		Site3D sNew,sOld;
		Mol *mOld,*mColl;
		std::pair<Site3D,std::pair<bool,SiteIt3D>> nbr_pair;
		bool occ_todo, occ_done;
		std::pair<bool,SiteIt3D> coll_done_pair;
		std::pair<bool,BiReaction*> rxn_pair;
		BiReaction *rxn; 

		// Go over all mols to move
		while (todo->size() > 0) {

			if (DIAG_DIFFUSE) { std::cout << std::endl; };

			// Reset
			occ_todo = false;
			occ_done = false;

			// Grab some element
			old_pair = todo->get_mol_random_it();

			sOldIt = old_pair.second;
			mOld = &(sOldIt.it_1D->second);

			if (DIAG_DIFFUSE) { std::cout << "diffuse_mols: got element..." << std::flush; };

			// Move
			nbr_pair = todo->get_neighbor_random(sOldIt);
			sNew = nbr_pair.first;

			if (DIAG_DIFFUSE) { std::cout << "got neighbor..." << std::flush; };

			// Check if occupied in todo lattice
			occ_todo = nbr_pair.second.first;
			if (occ_todo) {
				// Yes; its collided with something in the todo pile
				sNewIt = nbr_pair.second.second;
				mColl = &(sNewIt.it_1D->second);
			} else {
				// No; check if it's collided with something in the done pile
				coll_done_pair = _lattice->get_mol_it(sNew);
				occ_done = coll_done_pair.first;
				if (occ_done) {
					// Yes; its collided with something in the done pile
					sNewIt = coll_done_pair.second;
					mColl = &(sNewIt.it_1D->second);
				};
			};

			if (DIAG_DIFFUSE) { std::cout << "check colls..." << std::flush; };

			// If unoccupied, just commit the move (diffuse)
			if (!occ_todo && !occ_done) {
				// Move
				_lattice->make_mol(sNew, mOld->sp);
				// Continue
				todo->erase_mol_it(sOldIt);
				continue;
			};

			// Collision; check if reaction occurs
			rxn_pair = mOld->check_bi_rxns_mol(mColl);

			if (DIAG_DIFFUSE) { std::cout << "checked rxn..." << std::flush; };

			// No reaction occurred?
			if (!(rxn_pair.first)) {
				// No reaction & don't move
				_lattice->make_mol(Site3D(sOldIt), mOld->sp);
				// Continue
				todo->erase_mol_it(sOldIt);
				continue;
			};

			// Reaction occurred
			rxn = rxn_pair.second;
			rxn->count++;

			if (DIAG_DIFFUSE) { std::cout << "rxn occurred..." << std::flush; };

			// Remove the reactants
			sOld = Site3D(sOldIt); // grab it before it's erased
			todo->erase_mol_it(sOldIt);
			if (occ_todo) { 
				// Note: the old iterator sNewIt has been invalidated :(
				todo->erase_mol(sNew);
			} else if (occ_done) { 
				_lattice->erase_mol(sNewIt); 
			};

			if (DIAG_DIFFUSE) { std::cout << "removed r..." << std::flush; };

			// Conserve reactants
			if (rxn->r1->conserved) {
				_lattice->make_mol_random(rxn->r1);
			};
			if (rxn->r2->conserved) {
				_lattice->make_mol_random(rxn->r2);
			};

			if (DIAG_DIFFUSE) { std::cout << "conserved r..." << std::flush; };

			// Place products, if needed at the old site
			if (rxn->p.size() == 1) {
				_lattice->make_mol(sNew, rxn->p[0]);
			} else if (rxn->p.size() == 2) {
				_lattice->make_mol(sNew, rxn->p[0]);
				_lattice->make_mol(sOld, rxn->p[1]);
			};

			if (DIAG_DIFFUSE) { std::cout << "added p..." << std::flush; };

			// Conserve products
			for (auto p: rxn->p)
			{
				if (p->conserved) {
					_lattice->erase_mol_random(p);
				};
			};

			if (DIAG_DIFFUSE) { std::cout << "conserved p..." << std::flush; };

			// Finished
		};

		// Clear old
		delete todo;
	};


	/********************
	Run simulation
	********************/

	void Simulation::Impl::run(int n_timesteps, bool verbose, bool write_counts, bool write_nns, bool write_latt, int write_step, int write_version_no, std::string dir)
	{
		// Clear data in files if writing
		std::ofstream ofs;
		std::stringstream fname;
		if (write_counts || write_nns || write_latt) {
			// Clear counts
			for (auto s: this->_species) {
				fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/counts/" << s.name << ".txt";
				ofs.open(fname.str(), std::ofstream::out | std::ofstream::trunc);
				ofs.close();
				fname.str("");
			};
			// Clear nns
			auto it1 = this->_species.begin();
			auto it2 = this->_species.begin();
			while (it1 != this->_species.end()) {
				it2 = it1;
				while (it2 != this->_species.end()) {
					fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/nns/" << it1->name << "_" << it2->name << ".txt";
					ofs.open(fname.str(), std::ofstream::out | std::ofstream::trunc);
					ofs.close();
					fname.str("");
					it2++;
				};
				it1++;
			};
			// Clear lattice data
			fname << "exec rm -r ./" << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/lattice/*";
			system(fname.str().c_str());
			fname.str("");
		};

		// Go over all timesteps
		double t_next;
		for (int i_step=0; i_step < n_timesteps; i_step++) 
		{
			// The next time
			t_next = this->_t + this->_dt;

			// Print if needed
			if (this->_t_step % 1 == 0) {
				if (verbose) {
					std::cout << "Time: " << this->_t << " / " << n_timesteps*this->_dt << " " << this->_species << std::endl;
				};
			};

			// Write if needed
			if (this->_t_step % write_step == 0) {
				if (write_counts) {
					// Write counts to file
					for (auto s: this->_species) {
						fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/counts/" << s.name << ".txt";
						ofs.open(fname.str(), std::ofstream::app);
						ofs << this->_t << " " << s.count << "\n";
						ofs.close();
						fname.str("");
					};
				};
				if (write_nns) {
					// Write NNs to file
					auto it1 = this->_species.begin();
					auto it2 = this->_species.begin();
					while (it1 != this->_species.end()) {
						it2 = it1;
						while (it2 != this->_species.end()) {
							fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/nns/" << it1->name << "_" << it2->name << ".txt";
							ofs.open(fname.str(), std::ofstream::app);
							ofs << this->_t << " " << _lattice->get_nn(&(*it1),&(*it2)) << "\n";
							ofs.close();
							fname.str("");
							it2++;
						};
						it1++;
					};
				};
				if (write_latt) {
					// Write the lattice
					write_lattice(int(this->_t_step/write_step), write_version_no, dir);
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

	void Simulation::Impl::write_lattice(int index, int write_version_no, std::string dir)
	{
		std::stringstream fname;
		fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/lattice/" << std::setfill('0') << std::setw(4) << index << ".txt";
		_lattice->write_to_file(fname.str());
		fname.str("");
	};

	/****************************************
	Main simulation - PRIVATE
	****************************************/

	/********************
	Find a species by name
	********************/

	Species* Simulation::Impl::_find_species(std::string name) {
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

	void Simulation::Impl::_schedule_uni() {
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

	/****************************************
	Main simulation IMPL forwards
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	Simulation::Simulation(double dt, int box_length, int dim) : _impl(new Impl(dt,box_length,dim)) {};
	Simulation::Simulation(Simulation&& other) = default; // movable but no copies
    Simulation& Simulation::operator=(Simulation&& other) = default; // movable but no copies
	Simulation::~Simulation() = default;

	/********************
	Add species
	********************/

	void Simulation::add_species(std::string name, bool conserved) {
		_impl->add_species(name,conserved);
	};

	/********************
	 Add a reaction
	********************/

	// Unimolecular rxn
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r) {
		_impl->add_uni_rxn(name,kr,r);
	};
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r, std::string p) {
		_impl->add_uni_rxn(name,kr,r,p);
	};
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r, std::string p1, std::string p2) {
		_impl->add_uni_rxn(name,kr,r,p1,p2);
	};

	// Bimolecular rxn
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2) {
		_impl->add_bi_rxn(name,prob,r1,r2);
	};
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p) {
		_impl->add_bi_rxn(name,prob,r1,r2,p);
	};
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2) {
		_impl->add_bi_rxn(name,prob,r1,r2,p1,p2);
	};

	/********************
	Populate lattice
	********************/

	void Simulation::populate_lattice(std::map<std::string,int> counts) {
		_impl->populate_lattice(counts);
	};
	void Simulation::populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, int n_steps) {
		_impl->populate_lattice(h_dict,j_dict,n_steps);
	};
	void Simulation::populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, int n_steps) {
		_impl->populate_lattice(h_dict,j_dict,k_dict,n_steps);
	};
	void Simulation::populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>> &q_dict, int n_steps) {

		_impl->populate_lattice(h_dict,j_dict,k_dict,q_dict,n_steps);
	};

	/********************
	Run simulation
	********************/

	void Simulation::run(int n_timesteps, bool verbose, bool write_counts, bool write_nns, bool write_latt, int write_step, int write_version_no, std::string dir) {
		_impl->run(n_timesteps,verbose,write_counts,write_nns,write_latt,write_step,write_version_no, dir);
	};

	/********************
	Write lattice
	********************/

	void Simulation::write_lattice(int index, int write_version_no, std::string dir) {
		_impl->write_lattice(index,write_version_no, dir);
	};
};