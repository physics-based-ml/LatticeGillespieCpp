#include "lattice1d.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"

/************************************
* Namespace for Gillespie1D
************************************/

namespace Gillespie1D {

	/****************************************
	Structure to hold a lattice site iterator
	****************************************/

	SiteIt::SiteIt() {};
	SiteIt::SiteIt(lattice_map::iterator itIn) 
	{
		this->it = itIn;
	};
	std::ostream& operator<<(std::ostream& os, const SiteIt& sit)
	{
	    return os << sit.it->first;
	};

	/****************************************
	Struct to hold a lattice Site
	****************************************/

	// Constructor
	Site::Site() {};
	Site::Site(int xIn) {
		this->x = xIn;
	};
	Site::Site(SiteIt sit) {
		this->x = sit.it->first;
	};

	// Comparator
	bool operator <(const Site& a, const Site& b) {
    	return a.x < b.x;
	};
	bool operator==(const Site& a, const Site& b) {
		return a.x==b.x;
	}; 
	std::ostream& operator<<(std::ostream& os, const Site& s)
	{
	    return os << s.x;
	};

	/****************************************
	Lattice1D
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	// Constructor
	Lattice1D::Lattice1D(int box_length)
	{
		this->_box_length = box_length;
	};

	// Destructor
	Lattice1D::~Lattice1D() {};

	/********************
	Clear, size
	********************/

	void Lattice1D::clear() { this->_map.clear(); };
	int Lattice1D::size() { return this->_map.size(); };

	/********************
	Make a mol
	********************/

	std::pair<bool,SiteIt> Lattice1D::make_mol(Site s, Species *sp) 
	{
		// Check if the site is empty
		std::pair<bool,SiteIt> spair = get_mol_it(s);
		if (spair.first) {
			// Not empty
			return std::make_pair(false,SiteIt());
		};

		// Update counts
		sp->count++;

		// Make
		lattice_map::iterator it;
		it = this->_map.find(s.x);
		if (it == this->_map.end()) {
			this->_map.insert(std::make_pair(s.x,Mol(sp)));
		} else {
			return std::make_pair(false,SiteIt());
		};

		return std::make_pair(true,SiteIt(it));
	};

	std::pair<bool,SiteIt> Lattice1D::make_mol_random(Species *sp) 
	{
		// Get a random free site
		std::pair<bool,Site> spair = get_free_site();

		if (!(spair.first)) {
			// No free sites at all
			return std::make_pair(false,SiteIt());
		};

		// Make
		std::pair<bool,SiteIt> ret = make_mol(spair.second,sp);
		return ret;
	};

	/********************
	Erase a mol
	********************/

	bool Lattice1D::erase_mol(Site s) 
	{
		// Get an iterator to the site
		std::pair<bool,SiteIt> spair = get_mol_it(s);
		if (spair.first) {
			return erase_mol_it(spair.second);
		} else {
			return false;
		};
	};

	bool Lattice1D::erase_mol_it(SiteIt sit) 
	{
		// Update counts on species
		sit.it->second.sp->count--;

		// Erase
		this->_map.erase(sit.it);
		return true;
	};

	std::pair<bool,Site> Lattice1D::erase_mol_random(Species *sp) 
	{
		// Get an iterator to a random site
		std::pair<bool,SiteIt> spair = get_mol_random_it(sp);
		if (spair.first) {
			Site s = Site(spair.second.it->first);
			bool succ = erase_mol_it(spair.second);
			if (succ) {
				return std::make_pair(true,s);
			};
		};
		return std::make_pair(false,Site());
	};

	/********************
	Get a mol
	********************/

	std::pair<bool,SiteIt> Lattice1D::get_mol_it(Site s) 
	{
		auto it = this->_map.find(s.x);
		if (it != this->_map.end()) {
			return std::make_pair(true,SiteIt(it));
		};
		return std::make_pair(false,SiteIt());	
	};

	std::pair<bool,SiteIt> Lattice1D::get_mol_it(Site s, Species *sp) {
		std::pair<bool,SiteIt> ret = get_mol_it(s);
		if (ret.first) {
			if (ret.second.it->second.sp == sp) {
				return ret;
			};
		};
		return std::make_pair(false,SiteIt());
	};

	std::pair<bool,SiteIt> Lattice1D::get_mol_random_it() {
		// Get random indexes to search
		std::vector<int> idxs = _get_random_idxs();

		// Try all sites
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			spair = get_mol_it(Site(idxs[i1]));
			if (spair.first) 
			{
				return spair;
			};
		};
		return std::make_pair(false,SiteIt());
	};

	std::pair<bool,SiteIt> Lattice1D::get_mol_random_it(Species *sp) 
	{
		// Get random indexes to search
		std::vector<int> idxs = _get_random_idxs();

		// Try all sites
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			spair = get_mol_it(Site(idxs[i1]),sp);
			if (spair.first) 
			{
				return spair;
			};
		};
		return std::make_pair(false,SiteIt());
	};

	/********************
	Get a free site
	********************/

	std::pair<bool,Site> Lattice1D::get_free_site() 
	{
		// Get random indexes to search
		std::vector<int> idxs = _get_random_idxs();

		// Try all sites
		Site s;
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			s = Site(idxs[i1]);
			spair = get_mol_it(s);
			if (!(spair.first)) {
				return std::make_pair(true,s);
			};
		};

		return std::make_pair(false,s);
	};

	/********************
	Get neighbors of a site
	********************/

	std::pair<Site,std::pair<bool,SiteIt>> Lattice1D::get_neighbor_random(Site s)
	{
		// All neighbors
		std::vector<Site> nbrs = _get_all_neighbors(s);

		// Shuffle
		std::random_shuffle(nbrs.begin(),nbrs.end());

		// Random choice
		auto it = nbrs.begin();
		std::advance(it,randI(0,nbrs.size()-1));

		// Check if occupied
		std::pair<bool,SiteIt> occ = get_mol_it(*it);
		if (occ.first) {
			// Yes
			return std::make_pair(Site(occ.second),std::make_pair(true,occ.second));
		} else {
			// No
			return std::make_pair(*it,std::make_pair(false,SiteIt()));
		};
	};

	std::pair<Site,std::pair<bool,SiteIt>> Lattice1D::get_neighbor_random(SiteIt sit)
	{
		return get_neighbor_random(Site(sit));
	};

	std::pair<bool,Site> Lattice1D::get_free_neighbor_random(Site s) 
	{
		// All allowed nbrs
		std::vector<Site> nbrs = _get_all_neighbors(s);

		// Shuffle
		std::random_shuffle(nbrs.begin(),nbrs.end());

		// Check all neighbors
		std::pair<bool,SiteIt> spair;
		for (auto n: nbrs) {
			spair = get_mol_it(n);
			if (!(spair.first)) {
				return std::make_pair(true,n);
			};
		};

		return std::make_pair(false,Site());
	};

	std::pair<bool,Site> Lattice1D::get_free_neighbor_random(SiteIt sit) 
	{
		return get_free_neighbor_random(Site(sit));
	};

	/********************
	Get NN of species
	********************/

	int Lattice1D::get_nn(Species *sa, Species *sb)
	{
		int nn = 0;
		// Go through all sites
		Site s;
		std::pair<bool,SiteIt> spair;
		std::vector<Site> nbrs;
		for (auto i1=1; i1 <= this->_box_length; i1++) {
			s = Site(i1);
			spair = get_mol_it(s,sa);
			if (spair.first) {
				// This one is species A
				// Now search all neigbors
				nbrs = _get_all_neighbors(s);
				// Go through all neighbors
				for (auto nbr: nbrs) {
					spair = get_mol_it(nbr,sb);
					if (spair.first) {
						// Ok!
						nn++;
					};
				};
			};
		};
		// Double counting?
		if (sa == sb) {
			nn /= 2;
		};
		return nn;
	};

	/********************
	Write lattice to a file
	********************/

	void Lattice1D::write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (auto it = this->_map.begin(); it != this->_map.end(); it++) {
			f << it->first << " " << it->second.sp->name << "\n";
		};
		f.close();
	};

	/********************
	Anneal
	********************/

	void Lattice1D::anneal(std::map<Species*,double> &h_dict,std::map<Species*,std::map<Species*,double>> &j_dict, int n_steps) {

		// Declarations

		Site s;

		std::pair<bool,SiteIt> ret_site;
		Species* sp;

		std::vector<Site> nbrs;
		std::vector<Species*> nbrs_occ;
		std::vector<Site>::iterator it_nbr;
		std::pair<bool,SiteIt> ret_nbr;

		std::map<Species*,double>::iterator ith;

		double hOld,jOld,hNew,jNew,energy_diff;

		// Go through the steps
		for (int i=0; i<n_steps; i++) {

			// Pick a site to flip randomly
			s = Site(randI(1,_box_length));

			// Get occupied neighbors 
			nbrs_occ.clear();
			nbrs = _get_all_neighbors(s);
			it_nbr = nbrs.begin();
			while (it_nbr != nbrs.end()) {
				ret_nbr = get_mol_it(*it_nbr);
				if (ret_nbr.first) {
					// Occupied
					nbrs_occ.push_back(ret_nbr.second.it->second.sp);
					it_nbr++;
				} else {
					// Empty
					it_nbr = nbrs.erase(it_nbr);
				};
			};

			// Check if this site is occupied
			ret_site = get_mol_it(s);
			if (ret_site.first) {
				// Occupied, flip down
				sp = ret_site.second.it->second.sp;
				hOld = -h_dict[sp];
				jOld = 0.0;
				for (auto nbr: nbrs_occ) {
					jOld -= j_dict[sp][nbr];
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
				std::advance(ith, randI(0,h_dict.size()-1));
				sp = ith->first;
				// New couplings
				hNew = -ith->second;
				jNew = 0.0;
				for (auto nbr: nbrs_occ) {
					jNew -= j_dict[sp][nbr];
				};
			};

			// Energy difference
			energy_diff = hNew + jNew - hOld - jOld;
			if (energy_diff < 0.0 || exp(-energy_diff) > randD(0.0,1.0)) {
				// Accept the flip!
				if (ret_site.first) {
					// Occupied, flip down
					erase_mol_it(ret_site.second);
				} else {
					// Unoccupied, flip up
					make_mol(s,sp);
				};
			};
		};
	};


	/****************************************
	Lattice1D - PRIVATE
	****************************************/

	/********************
	Get random indexes
	********************/

	std::vector<int> Lattice1D::_get_random_idxs()
	{
		// Vectors of sites
		std::vector<int> x1(this->_box_length);
		std::iota(std::begin(x1), std::end(x1), 1);

		// Shuffle
		std::random_shuffle(x1.begin(),x1.end());

		// Return
		return x1;
	};

	/********************
	Get all neighbors of a site
	********************/

	std::vector<Site> Lattice1D::_get_all_neighbors(Site s)
	{
		std::vector<Site> nbrs;

		// 2 are possible
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

		return nbrs;
	};

};