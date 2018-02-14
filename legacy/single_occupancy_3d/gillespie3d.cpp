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

	double randD(double dMin, double dMax)
	{
	    return dMin + ((double)rand() / RAND_MAX) * (dMax - dMin);
	};
	int randI(int iMin, int iMax)
	{
		return iMin + rand() % (iMax - iMin + 1);
	};

	/********************
	General function to write vector to file
	********************/

	void write_vector_to_file(std::string fname, std::vector<int> v)
	{
		std::ofstream f;
		f.open (fname);
		for (auto i : v) {
			f << i << "\n";
		};
		f.close();
	};

	/****************************************
	Pair of Species
	****************************************/

	SpeciesPair::SpeciesPair(std::string s1In, std::string s2In) {
		if (s1In <= s2In) {
			this->s1 = s1In;
			this->s2 = s2In;
		} else {
			this->s1 = s2In;
			this->s2 = s1In;
		};
	};

	// Comparator
	bool operator <(const SpeciesPair& a, const SpeciesPair& b) {
    	return std::tie(a.s1, a.s2) < std::tie(b.s1, b.s2);
	};

	/****************************************
	Reaction
	****************************************/

	Reaction::Reaction(std::string name) {
		this->_name = name;
	};
	Reaction::Reaction(std::string name, std::string p) {
		this->_name = name;
		this->_p.push_back(p);
	};
	Reaction::Reaction(std::string name, std::string p1, std::string p2) {
		this->_name = name;
		this->_p.push_back(p1);
		this->_p.push_back(p2);
		std::sort(this->_p.begin(),this->_p.end());
	};

	Reaction::~Reaction() {};

	/********************
	Getters/Setters
	********************/

	const std::string Reaction::name() {
		return this->_name;
	};

	const std::vector<std::string> Reaction::get_p() {
		return this->_p;
	};

	/********************
	Get/set/clear conserved
	********************/

	std::vector<std::string> Reaction::get_r_cons() { return this->_r_c; };
	std::vector<std::string> Reaction::get_p_cons() { return this->_p_c; };
	void Reaction::set_r_cons(std::string s) { this->_r_c.push_back(s); };
	void Reaction::set_p_cons(std::string s) { this->_p_c.push_back(s); };
	void Reaction::clear_r_cons() { this->_r_c.clear(); };
	void Reaction::clear_p_cons() { this->_p_c.clear(); };

	/****************************************
	UniReaction
	*****************************************/

	/********************
	Constructor/Destructor
	********************/

	UniReaction::UniReaction(std::string name, double kr, std::string r) : Reaction(name) {
		this->_kr = kr;
		this->_r = r;
	};
	UniReaction::UniReaction(std::string name, double kr, std::string r, std::string p) : Reaction(name, p) {
		this->_kr = kr;
		this->_r = r;
	};
	UniReaction::UniReaction(std::string name, double kr, std::string r, std::string p1, std::string p2) : Reaction(name, p1, p2) {
		this->_kr = kr;
		this->_r = r;
	};
	UniReaction::~UniReaction() {};

	/********************
	Getters
	********************/

	const std::string UniReaction::get_r() {
		return this->_r;
	};

	const double UniReaction::get_rate() {
		return this->_kr;
	}

	/********************
	Parent
	********************/

	// Get the name
	const std::string UniReaction::name() { return Reaction::name(); };

	// Get products
	const std::vector<std::string> UniReaction::get_p() { return Reaction::get_p(); };

	// Get/set/clear conserved
	std::vector<std::string> UniReaction::get_r_cons() { return Reaction::get_r_cons(); };
	std::vector<std::string> UniReaction::get_p_cons() { return Reaction::get_p_cons(); };
	void UniReaction::set_r_cons(std::string s) { Reaction::set_r_cons(s); };
	void UniReaction::set_p_cons(std::string s) { Reaction::set_p_cons(s); };
	void UniReaction::clear_r_cons() { Reaction::clear_r_cons(); };
	void UniReaction::clear_p_cons() { Reaction::clear_p_cons(); };

	/****************************************
	BiReaction
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	BiReaction::BiReaction(std::string name, double prob, std::string r1, std::string r2) : _r(r1,r2), Reaction(name) {
		this->_prob = prob;
	};
	BiReaction::BiReaction(std::string name, double prob, std::string r1, std::string r2, std::string p) : _r(r1,r2), Reaction(name,p) {
		this->_prob = prob;
	};
	BiReaction::BiReaction(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2) : _r(r1,r2), Reaction(name,p1,p2) {
		this->_prob = prob;
	};

	BiReaction::~BiReaction() {};

	/********************
	Getters
	********************/

	const SpeciesPair BiReaction::get_r() {
		return this->_r;
	};

	const double BiReaction::get_prob() {
		return this->_prob;
	}

	/********************
	Parent
	********************/

	// Get the name
	const std::string BiReaction::name() { return Reaction::name(); };

	// Get products
	const std::vector<std::string> BiReaction::get_p() { return Reaction::get_p(); };

	// Get/set/clear conserved
	std::vector<std::string> BiReaction::get_r_cons() { return Reaction::get_r_cons(); };
	std::vector<std::string> BiReaction::get_p_cons() { return Reaction::get_p_cons(); };
	void BiReaction::set_r_cons(std::string s) { Reaction::set_r_cons(s); };
	void BiReaction::set_p_cons(std::string s) { Reaction::set_p_cons(s); };
	void BiReaction::clear_r_cons() { Reaction::clear_r_cons(); };
	void BiReaction::clear_p_cons() { Reaction::clear_p_cons(); };

	/********************
	Struct to hold a lattice Site
	********************/

	// Constructor
	Site::Site(int xIn, int yIn, int zIn) {
		this->x = xIn;
		this->y = yIn;
		this->z = zIn;
	}

	// Comparator
	bool operator <(const Site& a, const Site& b) {
    	return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z);
	};
	bool operator==(const Site& a, const Site& b) {
		return std::tie(a.x, a.y, a.z) == std::tie(b.x, b.y, b.z);
	}; 
	std::ostream& operator<<(std::ostream& os, const Site& s)
	{
	    return os << s.x << " " << s.y << " " << s.z;
	}

	/****************************************
	Lattice
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	// Constructor
	Lattice::Lattice(int box_length)
	{
		this->_box_length = box_length;
	};

	// Destructor
	Lattice::~Lattice() {};

	/********************
	Erase/create sites
	********************/

	void Lattice::erase_site(Site s) 
	{
		lattice_map::iterator it = this->_map.find(s);
		if (it != this->_map.end()) 
		{
			// Update counts
			this->_counts[it->second]--;
			// Erase
			this->_map.erase(it);
		};
		return;
	};

	void Lattice::make_site(Site s, std::string species) 
	{
		// Update counts
		this->_counts[species]++;
		// Insert
		this->_map[s] = species;
		return;
	};

	/********************
	// Get a list of all occupied sites
	********************/

	std::vector<Site> Lattice::get_sites() {
		std::vector<Site> s;
		for (auto it: this->_map) { s.push_back(it.first); };
		return s;
	};

	/********************
	Get a random free site
	********************/

	std::pair<bool,Site> Lattice::get_random_free_site() {
		// Grab any random site
		Site s(randI(1,this->_box_length),randI(1,this->_box_length),randI(1,this->_box_length));

		// Timeout
		int ctr = 0;

		// Occupied?
		bool occ = get_site(s).first;
		std::vector<Site> nbrs;
		while (occ) {
			// Try to get an unoccupied neighbor
			nbrs = get_neighbors_empty(s);
			if (nbrs.size() > 0) {
				// Done
				occ = false;
				s = nbrs[0];
			} else {
				// New site
				s = Site(randI(1,this->_box_length),randI(1,this->_box_length),randI(1,this->_box_length));
				occ = get_site(s).first;
				ctr++;
			};

			if (ctr >= 1000) { 
				std::cout << "ERROR! COULD NOT FIND A FREE SITE\n";
				return std::make_pair(false,s);
			};
		};

		return std::make_pair(true,s);
	};

	/********************
	Get an occupied Site, if it exists
	********************/

	std::pair<bool,std::string> Lattice::get_site(Site s) {
		lattice_map::iterator it;
		it = this->_map.find(s);
		if (it != this->_map.end()) {
			return std::make_pair(true,it->second);
		} else {
			return std::make_pair(false,"");
		};
	};

	/********************
	Find a site by species
	********************/

	std::pair<bool,Site> Lattice::find_site(std::string species) {
		lattice_map::iterator it = this->_map.begin();
		while (it->second != species) {
			it++;
			if (it == this->_map.end()) {
				// Fail
				return std::make_pair(false,Site(0,0,0));
			};
		};
		return std::make_pair(true,it->first);
	};

	/********************
	Get neighboring lattice Sites
	********************/
	std::vector<Site> Lattice::get_neighbors(Site s)
	{
		std::vector<Site> nbrs;

		// 6 are possible
		Site nbr_Site = s;

		// x
		nbr_Site.x += 1;
		if (nbr_Site.x <= this->_box_length) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.x -= 2;
		if (nbr_Site.x >= 1) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.x += 1;
		// y
		nbr_Site.y += 1;
		if (nbr_Site.y <= this->_box_length) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.y -= 2;
		if (nbr_Site.y >= 1) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.y += 1;
		// z
		nbr_Site.z += 1;
		if (nbr_Site.z <= this->_box_length) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.z -= 2;
		if (nbr_Site.z >= 1) {
			nbrs.push_back(nbr_Site);
		};
		nbr_Site.z += 1;

		return nbrs;
	};
	std::vector<Site> Lattice::get_neighbors_empty(Site s)
	{
		// Get all neighbors
		std::vector<Site> nbrs = get_neighbors(s);
		// Discard the occupied ones
		std::vector<Site>::iterator it = nbrs.begin();
		while (it != nbrs.end()) {
			if (get_site(*it).first) {
				it = nbrs.erase(it);
			} else {
				it++;
			};
		};
		return nbrs;
	};
	std::vector<Site> Lattice::get_neighbors_occupied(Site s) {
		// Get all neighbors
		std::vector<Site> nbrs = get_neighbors(s);
		// Discard the occupied ones
		std::vector<Site>::iterator it = nbrs.begin();
		while (it != nbrs.end()) {
			if (!(get_site(*it).first)) {
				it = nbrs.erase(it);
			} else {
				it++;
			};
		};
		return nbrs;
	};

	/********************
	Get number particles, NN
	********************/

	std::map<std::string,int> Lattice::get_n_particles_all() {
		return this->_counts;
	};

	std::pair<bool,int> Lattice::get_n_particles_safe(std::string species) {
		std::map<std::string,int>::iterator it = this->_counts.find(species);
		if (it != this->_counts.end()) {
			return std::make_pair(true,it->second);
		} else {
			return std::make_pair(false,0);
		};
	};

	int Lattice::get_n_particles_total() {
		return this->_map.size();
	};

	int Lattice::get_n_particles(std::string s_name) {
		std::pair<bool,int> ret = get_n_particles_safe(s_name);
		if (ret.first) { 
			return ret.second; 
		} else {
			return 0;
		};
	};

	int Lattice::get_nn() {
		int n=0;
		std::vector<Site> nbrs;
		for (auto it = this->_map.begin(); it != this->_map.end(); it++) {
			// Get neighbors
			nbrs = get_neighbors_occupied(it->first);
			n += nbrs.size();
		};
		return n/2;
	};

	int Lattice::get_nn(SpeciesPair s_pair) {
		int n=0;
		std::vector<Site> nbrs;
		for (auto it = this->_map.begin(); it != this->_map.end(); it++) {
			if (it->second == s_pair.s1 || it->second == s_pair.s2) {
				// Get neighbors
				nbrs = get_neighbors_occupied(it->first);
				for (auto nbr: nbrs) {
					if ((it->second == s_pair.s1 && this->_map[nbr] == s_pair.s2) || (it->second == s_pair.s2 && this->_map[nbr] == s_pair.s1))
					{ 
					 	n += 1; 
					};
				};
			};
		};
		return n/2;
	};

	/********************
	Populate a lattice randomly
	********************/

	void Lattice::populate_lattice(std::map<std::string,int> counts) {
		std::pair<bool,Site> ret(true,Site(0,0,0));

	    // Go through all species
	    for (auto it: counts) {
			// Initialize counts for this species
			this->_counts[it.first] = 0;

	    	// All mols to place of this species
	    	for (auto i=0; i<it.second; i++) {
	    		// Random position
	    		ret = get_random_free_site();
	    		if (!(ret.first)) {
	    			std::cout << "ERROR! COULD NOT RANDOMLY POPULATE LATTICE" << std::endl;
	    			return;
	    		};
				// Place
				this->_map[ret.second] = it.first;
				// Count
				this->_counts[it.first]++;
	    	};
	    };

	    return;
	};

	/********************
	Populate a lattice using annealing
	********************/

	void Lattice::populate_lattice(std::map<std::string,double> &h_dict,std::map<SpeciesPair,double> &j_dict, int n_annealing_steps, bool write_lattice, bool write_statistics)
	{
		// Random number of initial particles (min is 1, max is box vol)
		int n = randI(1, this->_box_length*this->_box_length*this->_box_length);

		// Random initial counts
		int n_possible = pow(this->_box_length,3);
		std::map<std::string,int> counts0;
		for (auto hpr : h_dict) {
			counts0[hpr.first] = randI(0,n_possible);
			n_possible -= counts0[hpr.first];
			if (n_possible < 0) { n_possible = 0; };
		};			

		// Random initial lattice
		populate_lattice(counts0);

		// Write if needed
		if (write_lattice) {
			write_to_file("annealing_0000.txt");
		};

		// Track and later write statistics, if needed
		std::map<std::string,std::vector<int>> stat_n;
		std::map<SpeciesPair,std::vector<int>> stat_nn;
		if (write_statistics) {
			// Go through all species
			for (auto its = h_dict.begin(); its != h_dict.end(); its++) {
				stat_n[its->first].push_back(get_n_particles(its->first));
			};
			for (auto itj = j_dict.begin(); itj != j_dict.end(); itj++) {
				stat_nn[itj->first].push_back(get_nn(itj->first));
			};
		};

		// Do the annealing
		int ix=0,iy=0,iz=0,r_idx=0;
		Site s(ix,iy,iz);
		lattice_map::iterator itl;
		double hOld,hNew,jOld,jNew;
		std::string s_name;
		double energy_diff;
		std::stringstream fname;
		std::map<std::string,double>::iterator ith;
		for (int i=0; i<n_annealing_steps; i++) {

			// Pick a site to flip randomly
			ix = randI(1,this->_box_length);
			iy = randI(1,this->_box_length);
			iz = randI(1,this->_box_length);
			s = Site(ix,iy,iz);

			// Get occupied neighbors 
			std::vector<Site> nbrs = get_neighbors_occupied(s);

			// Check if the site is occupied
			itl = this->_map.find(s);
			if (itl != this->_map.end()) {
				// Occupied, flip down
				hOld = -h_dict[itl->second];
				jOld = 0.0;
				s_name = this->_map[s];
				for (auto nbr: nbrs) {
					jOld -= j_dict[SpeciesPair(s_name,this->_map[nbr])];
				};
				// New couplings
				hNew = 0.0;
				jNew = 0.0;
			} else {
				// Unoccupied, flip up
				hOld = 0.0;
				jOld = 0.0;
				// Random species
				ith = h_dict.begin();
				r_idx = randI(0,h_dict.size()-1); // 0 to size-1
				std::advance(ith, r_idx);
				s_name = ith->first;
				// New couplings
				hNew = -ith->second;
				jNew = 0.0;
				for (auto nbr: nbrs) {
					jNew -= j_dict[SpeciesPair(s_name,this->_map[nbr])];
				};
			};

			// Energy difference
			energy_diff = hNew + jNew - hOld - jOld;
			if (energy_diff < 0.0 || exp(-energy_diff) > randD(0.0,1.0)) {
				// Accept the flip!
				if (itl != this->_map.end()) {
					// Occupied, flip down
					this->_map.erase(itl);
					// Update count
					this->_counts[s_name]--;
				} else {
					// Unoccupied, flip up
					this->_map[s] = s_name;
					// Update count
					this->_counts[s_name]++;
				};
			};

			// Store statistics if needed
			if (write_statistics) {
				// Go through all species
				for (auto its = h_dict.begin(); its != h_dict.end(); its++) {
					stat_n[its->first].push_back(get_n_particles(its->first));
				};
				for (auto itj = j_dict.begin(); itj != j_dict.end(); itj++) {
					stat_nn[itj->first].push_back(get_nn(itj->first));
				};
			};

			// Write the lattice if needed
			if (write_lattice) {
				fname << "annealing_lattice_" << std::setfill('0') << std::setw(4) << i+1 << ".txt";
				write_to_file(fname.str());
				fname.str("");
			};
		};

		// Write statistics if needed
		if (write_statistics) {
			for (auto it: stat_n) {	
				fname << "annealing_stats_" << it.first << ".txt";
				write_vector_to_file(fname.str(),it.second);
				fname.str("");
			};
			for (auto it: stat_nn) {	
				fname << "annealing_stats_" << it.first.s1 << "_" << it.first.s2 << ".txt";
				write_vector_to_file(fname.str(),it.second);
				fname.str("");
			};
		};

		return;
	};

	/********************
	Write lattice to a file
	********************/

	void Lattice::write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		lattice_map::iterator ite = this->_map.end();
		ite--;
		for (auto it = this->_map.begin(); it != this->_map.end(); it++) {
			f << it->first.x << " " << it->first.y << " " << it->first.z << " " << it->second;
			if (it != ite) {
				f << "\n";
			};
		};
		f.close();
	};

	/****************************************
	Main simulation
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	// Constructors
	Simulation::Simulation(Lattice l0, double dt) : _lattice(0)
	{
		this->_lattice = l0;
		this->_uni_next = NULL;
		this->_t = 0.0;
		this->_dt = dt;
	};

	// Destructor
	Simulation::~Simulation() {};

	/********************
	Set conserved quantities
	********************/

	void Simulation::set_conserved(std::string species) {
		this->_conserved_quantities.push_back(species);
	};

	/********************
	 Add a reaction
	********************/

	void Simulation::add_reaction(UniReaction rxn) {		
		this->_unimol_rxn_dict[rxn.get_r()].push_back(rxn);
		// Init rxn,collision count
		this->_rxn_count[rxn.name()].push_back(0);
	};
	void Simulation::add_reaction(BiReaction rxn) {	
		this->_bimol_rxn_dict[rxn.get_r()].push_back(rxn);
		// Init rxn,collision count
		this->_rxn_count[rxn.name()].push_back(0);
		this->_coll_count[rxn.name()].push_back(0);
	};

	/********************
	Diffuse all the mols and do bimol reactions
	********************/
	void Simulation::_diffuse_mols() 
	{
		// List of mols to move
		std::vector<Site> todo = this->_lattice.get_sites();
		// std::cout << "Diffusing " << todo.size() << " mols... ";

		// Shuffle
		std::random_shuffle ( todo.begin(), todo.end() );

		// Iterate over all mols
		std::vector<Site>::iterator it,itn;
		std::vector<Site> nbrs;
		int r_idx,irxn;
		Site s_new(0,0,0), s_old(0,0,0);
		std::pair<bool,std::string> occ_check;
		std::string species_old,species_new;
		std::map<SpeciesPair,std::vector<BiReaction>>::iterator itr;
		double r_prob,cum;
		std::vector<double> probs;
		std::vector<std::string> p;
		int n_moved=0;
		BiReaction rxn("",0.0,"","");
		std::pair<bool,Site> ret(true,Site(0,0,0));

		it = todo.begin();
		while (it != todo.end())
		{
			// Mol
			s_old = *it;
			species_old = this->_lattice.get_site(s_old).second;

			// Move
			nbrs = this->_lattice.get_neighbors(s_old);
			itn = nbrs.begin();
			r_idx = randI(0,nbrs.size()-1); // 0 to size-1
			std::advance(itn, r_idx);
			s_new = *itn;

			// Check if occupied
			occ_check = this->_lattice.get_site(s_new);

			if (! occ_check.first ) {
				// Not occupied; commit the move
				this->_lattice.erase_site(s_old);		
				this->_lattice.make_site(s_new, species_old);
				n_moved++;
				// std::cout << "Moved " << s_old << " to " << s_new << std::endl;
				it = todo.erase(it);				
				continue;
			};			

			// It is occupied - get the species
			species_new = this->_lattice.get_site(s_new).second;

			// std::cout << "It is occupied; checking reactions between: " << species_old << " and " << species_new << "\n";
			// Check for bimol reactions
			itr = this->_bimol_rxn_dict.find(SpeciesPair(species_old,species_new));
			if (itr == this->_bimol_rxn_dict.end())
			{
				// std::cout << "No reaction possible\n";
				// No reaction possible; reject the move
				it = todo.erase(it);
				continue;
			};

			// Count collision
			if (itr->second.size() != 1) { std::cout << "!!!!!!!!!!!!!!!" << std::endl;};
			this->_coll_count[itr->second[0].name()].back()++;

			// Reactions are possible - choose one from probabilities
			probs.clear();
			cum = 0.0;
			for (auto rxn: itr->second) {
				cum += rxn.get_prob();
				probs.push_back(cum);
			};
			probs.push_back(probs.size());
			r_prob = randD(0.0,probs.size()-1.0);
			irxn=0;
			while (probs[irxn] < r_prob) {irxn++;};

			if (irxn == probs.size()-1)
			{
				// All reactions failed; reject the move
				it = todo.erase(it);
				continue;
			};

			// Erase the old sites
			this->_lattice.erase_site(s_old);
			this->_lattice.erase_site(s_new);	

			// We moved
			n_moved++;

			// Place the products
			p = rxn.get_p();
			if (p.size() == 1)
			{
				this->_lattice.make_site(s_new, p[0]);						
			} else if (p.size() == 2)
			{
				this->_lattice.make_site(s_old, p[0]);	
				this->_lattice.make_site(s_new, p[1]);						
			};

			// Enforce conservation laws upon rectants/products
			for (auto r_c: rxn.get_r_cons()) {
				// Add reactants back
				ret = this->_lattice.get_random_free_site();
				if (ret.first) {
					this->_lattice.make_site(ret.second,r_c);
				};
			};
			for (auto p_c: rxn.get_p_cons()) {
				// Remove products
				ret = this->_lattice.find_site(p_c);
				if (ret.first) {
					this->_lattice.erase_site(ret.second);
				};
			};

			// Count reaction
			this->_rxn_count[itr->second[irxn].name()].back()++;

			// Advance iterator
			it = todo.erase(it);
			todo.erase(std::remove(todo.begin(), todo.end(), s_new), todo.end());
			it = todo.begin();	
		};
		// std::cout << "finished - moved " << n_moved << "\n";
	};

	/********************
	Schedule the next uni reaction
	********************/

	void Simulation::_schedule_uni() {
		double props_cum = 0.0;
		std::vector<double> props;
		props.push_back(0.0);
		std::vector<UniReaction*> rxns;

		// Go through all possible reagants, calculate propensities
		int n;
		for (auto u: this->_unimol_rxn_dict) 
		{
			// Count number
			n = this->_lattice.get_n_particles(u.first);

			// Go through all reactions with this reagant
			for (int r = 0; r < u.second.size(); r++)
			{
				props_cum += n * u.second[r].get_rate();
				props.push_back(props_cum);
				rxns.push_back(&this->_unimol_rxn_dict[u.first][r]);
			};
		};

		// Check that at least one reaction is possible
		if (!(props_cum > 0)) {
			this->_t_uni_next = this->_t;
			this->_uni_next = NULL; // Based on NULL, will check again later
			return;
		};

		// Choose a reaction
		double r = randD(0.0,props_cum);
		n = 0;
		while (!(props[n] < r && r < props[n+1])) {
			n++;
		};
		this->_uni_next = rxns[n];

		// Time of next reaction
		this->_t_uni_next = this->_t + log(1.0/randD(0.0,1.0))/props_cum;
	};

	/********************
	Do the scheduled reaction
	********************/

	void Simulation::_do_uni_rxn() {
		std::pair<bool,Site> ret(true,Site(0,0,0));
		std::vector<std::string> p;
		std::vector<Site> nbrs;
		std::string r;
		if (this->_uni_next != NULL) {
			// Grab a random mol
			r = this->_uni_next->get_r();
			ret = this->_lattice.find_site(r);

			if (!(ret.first)) {
				// Could not find a mol...
				this->_uni_next = NULL;
				return;
			};

			// Check if there is room for the products, and place them
			p = this->_uni_next->get_p();
			if (p.size() == 0) {
				this->_lattice.erase_site(ret.second);
			} else if (p.size() == 1) {
				this->_lattice.erase_site(ret.second);
				this->_lattice.make_site(ret.second, p[0]);	
			} else if (p.size() == 2) {
				nbrs = this->_lattice.get_neighbors_empty(ret.second);
				if (nbrs.size() < 1) {
					// Not enough room; fail
					this->_uni_next = NULL;
					return;
				};
				this->_lattice.erase_site(ret.second);
				this->_lattice.make_site(ret.second, p[0]);	
				this->_lattice.make_site(nbrs[0], p[1]);				
			};

			// Enforce conservation laws upon rectants/products
			for (auto r_c: this->_uni_next->get_r_cons()) {
				// Add reactants back
				ret = this->_lattice.get_random_free_site();
				if (ret.first) {
					this->_lattice.make_site(ret.second,r_c);
				};
			};
			for (auto p_c: this->_uni_next->get_p_cons()) {
				// Remove products
				ret = this->_lattice.find_site(p_c);
				if (ret.first) {
					this->_lattice.erase_site(ret.second);
				};
			};

			// Advance time
			this->_t = this->_t_uni_next;

			// Count it
			this->_rxn_count[this->_uni_next->name()].back()++;
		};
	};

	/********************
	Enforce conservation laws upon reactions
	********************/

	void Simulation::_enforce_conservation_in_rxns() {
		std::vector<std::string>::iterator it;

		// Go through all uni reactions
		for (auto rlist = this->_unimol_rxn_dict.begin(); rlist != this->_unimol_rxn_dict.end(); rlist++) {
			for (auto rxn = rlist->second.begin(); rxn != rlist->second.end(); rxn++) {
				// Clear existing
				rxn->clear_p_cons();
				rxn->clear_r_cons();
				// Enforce reactant
				it = std::find(this->_conserved_quantities.begin(),this->_conserved_quantities.end(),rxn->get_r());
				if (it != this->_conserved_quantities.end()) {
					rxn->set_r_cons(rxn->get_r());
				};
				// Enforce products
				for (auto p: rxn->get_p()) {
					it = std::find(this->_conserved_quantities.begin(),this->_conserved_quantities.end(),p);
					if (it != this->_conserved_quantities.end()) {
						rxn->set_p_cons(p);
					};
				};
			};
		};

		// Go through all bi reactions
		SpeciesPair r("","");
		for (auto rlist = this->_bimol_rxn_dict.begin(); rlist != this->_bimol_rxn_dict.end(); rlist++) {
			for (auto rxn = rlist->second.begin(); rxn != rlist->second.end(); rxn++) {
				// Clear existing
				rxn->clear_p_cons();
				rxn->clear_r_cons();
				// Enforce reactants
				r = rxn->get_r();
				it = std::find(this->_conserved_quantities.begin(),this->_conserved_quantities.end(),r.s1);
				if (it != this->_conserved_quantities.end()) {
					rxn->set_r_cons(r.s1);
				};
				it = std::find(this->_conserved_quantities.begin(),this->_conserved_quantities.end(),r.s2);
				if (it != this->_conserved_quantities.end()) {
					rxn->set_r_cons(r.s2);
				};
				// Enforce products
				for (auto p: rxn->get_p()) {
					it = std::find(this->_conserved_quantities.begin(),this->_conserved_quantities.end(),p);
					if (it != this->_conserved_quantities.end()) {
						rxn->set_p_cons(p);
					};
				};
			};
		};
	};

	/********************
	Run simulation
	********************/

	void Simulation::run(int n_timesteps, bool verbose, bool write_statistics)
	{
		// Enforce conservation laws in the reactions first
		_enforce_conservation_in_rxns();

		// Initialize the counts
		for (auto itc : this->_lattice.get_n_particles_all()) {
			this->_n_count[itc.first].push_back(itc.second);
		};

		// Go over all timesteps
		double t_next;
		for (int i_step=0; i_step < n_timesteps; i_step++) 
		{
			// The next time
			t_next = this->_t + this->_dt;

			// Print if needed
			if (verbose) {
				if (int(this->_t) % 20 == 0) {
					std::cout << "Timestep: " << this->_t << " no Z: " << this->_lattice.get_n_particles("Z") << "\n" << std::flush;
				};
			};

			// Advance the reaction,collision counters
			for (auto x: this->_rxn_count) {
				this->_rxn_count[x.first].push_back(x.second.back());
			};
			for (auto x: this->_coll_count) {
				this->_coll_count[x.first].push_back(x.second.back());
			};

			// Do we need to schedule unimolecular reactions?
			// (Either start of sim, or there were none available before)
			if (this->_uni_next == NULL) {
				// Schedule unimolecular reactions
				_schedule_uni();
			};

			// Do uni reactions
			while (this->_uni_next != NULL && this->_t_uni_next < t_next) {
				// Do it
				_do_uni_rxn();
				// Schedule
				_schedule_uni();
			};

			// Diffuse and do bimol reactions
			_diffuse_mols();

			// Store the new counts
			for (auto itc : this->_lattice.get_n_particles_all()) {
				this->_n_count[itc.first].push_back(itc.second);
			};

			// It is now the next time
			this->_t = t_next;

			// Stop and print
			if (this->_t > 0.001) {
				for (auto rc: this->_coll_count) {
					std::cout << "-----" << std::endl;
					std::cout << "Reaction: " << rc.first << std::endl;
					for (auto tmp: this->_bimol_rxn_dict) {
						for (auto tmp2: tmp.second) {
							if (tmp2.name() == rc.first) {
								std::cout << "Expected prob: " << tmp2.get_prob() << std::endl;
							};
						};
					};
					std::cout << "Collisions: " << rc.second.back() << std::endl;
					std::cout << "Reactions: " << this->_rxn_count[rc.first].back() << std::endl;
					std::cout << "Actual prob: " << 1.0 * this->_rxn_count[rc.first].back() / rc.second.back() << std::endl;
				};
				break;
			};
		};

		// Write statistics
		if (write_statistics) {
			write_vector_to_file("X.txt",this->_n_count["X"]);
			write_vector_to_file("Y.txt",this->_n_count["Y"]);
			write_vector_to_file("Z.txt",this->_n_count["Z"]);
			std::stringstream fname;
			for (auto it = this->_rxn_count.begin(); it != this->_rxn_count.end(); it++) {
				fname << it->first << ".txt";
				write_vector_to_file(fname.str(),it->second);
				fname.str("");
			};

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

	/********************
	Get reaction, coll counts
	********************/

	std::map<std::string,std::vector<int>> Simulation::get_rxn_count() {
		return this->_rxn_count;
	};
	std::map<std::string,std::vector<int>> Simulation::get_coll_count() {
		return this->_coll_count;
	};


};