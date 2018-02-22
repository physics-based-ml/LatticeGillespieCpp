#include "lattice.hpp"
#include <iostream>
#include <fstream>
#include <numeric>

/************************************
* Namespace for Gillespie3D
************************************/

namespace Gillespie3D {

	/****************************************
	Structure to hold a lattice site iterator
	****************************************/

	SiteIt::SiteIt() {};
	SiteIt::SiteIt(lattice_map::iterator itIn, lattice_map_1::iterator it_1In, lattice_map_2::iterator it_2In) 
	{
		this->it = itIn;
		this->it_1 = it_1In;
		this->it_2 = it_2In;
	};
	std::ostream& operator<<(std::ostream& os, const SiteIt& sit)
	{
	    return os << sit.it->first << " " << sit.it_1->first << " " << sit.it_2->first;
	};

	/****************************************
	Struct to hold a lattice Site
	****************************************/

	// Constructor
	Site::Site() {};
	Site::Site(int xIn, int yIn, int zIn) {
		this->x = xIn;
		this->y = yIn;
		this->z = zIn;
	};
	Site::Site(SiteIt sit) {
		this->x = sit.it->first;
		this->y = sit.it_1->first;
		this->z = sit.it_2->first;
	};

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
	};

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
	Clear, size
	********************/

	void Lattice::clear() { this->_map.clear(); };
	int Lattice::size() { return this->_map.size(); };

	/********************
	Make a mol
	********************/

	std::pair<bool,SiteIt> Lattice::make_mol(Site s, Species *sp) 
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
		lattice_map_1::iterator it_1;
		lattice_map_2::iterator it_2;
		it = this->_map.find(s.x);
		if (it == this->_map.end()) {
			auto ret = this->_map.insert(std::make_pair(s.x,lattice_map_1()));
			it = ret.first;
		};
		it_1 = it->second.find(s.y);
		if (it_1 == it->second.end()) {
			auto ret_1 = it->second.insert(std::make_pair(s.y,lattice_map_2()));
			it_1 = ret_1.first;
		};
		auto ret_2 = it_1->second.insert(std::make_pair(s.z,Mol(sp)));
		it_2 = ret_2.first;

		return std::make_pair(true,SiteIt(it,it_1,it_2));
	};

	std::pair<bool,SiteIt> Lattice::make_mol_random(Species *sp) 
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

	bool Lattice::erase_mol(Site s) 
	{
		// Get an iterator to the site
		std::pair<bool,SiteIt> spair = get_mol_it(s);
		if (spair.first) {
			return erase_mol_it(spair.second);
		} else {
			return false;
		};
	};

	bool Lattice::erase_mol_it(SiteIt sit) 
	{
		// Update counts on species
		sit.it_2->second.sp->count--;

		// Erase inner_2 from inner_1
		sit.it_1->second.erase(sit.it_2);
		if (sit.it_1->second.size() == 0)
		{
			// Erase inner_1 from outer
			sit.it->second.erase(sit.it_1);
			if (sit.it->second.size() == 0)
			{
				// Erase outer from map
				this->_map.erase(sit.it);
				return true;
			};
		};
		return true; // Nonsense; it can't fail?!
	};

	std::pair<bool,Site> Lattice::erase_mol_random(Species *sp) 
	{
		// Get an iterator to a random site
		std::pair<bool,SiteIt> spair = get_mol_random_it(sp);
		if (spair.first) {
			Site s = Site(spair.second.it->first,spair.second.it_1->first,spair.second.it_2->first);
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

	std::pair<bool,SiteIt> Lattice::get_mol_it(Site s) 
	{
		auto it = this->_map.find(s.x);
		if (it != this->_map.end()) {
			auto it_1 = it->second.find(s.y);
			if (it_1 != it->second.end()) {
				auto it_2 = it_1->second.find(s.z);
				if (it_2 != it_1->second.end()) {
					return std::make_pair(true,SiteIt(it,it_1,it_2));
				};
			};
		};
		return std::make_pair(false,SiteIt());	
	};

	std::pair<bool,SiteIt> Lattice::get_mol_it(Site s, Species *sp) {
		std::pair<bool,SiteIt> ret = get_mol_it(s);
		if (ret.first) {
			if (ret.second.it_2->second.sp == sp) {
				return ret;
			};
		};
		return std::make_pair(false,SiteIt());
	};

	std::pair<bool,SiteIt> Lattice::get_mol_random_it() {
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = _get_random_idxs();

		// Try all sites
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					spair = get_mol_it(Site(idxs[0][i1],idxs[1][i2],idxs[2][i3]));
					if (spair.first) 
					{
						return spair;
					};
				};
			};
		};
		return std::make_pair(false,SiteIt());
	};

	std::pair<bool,SiteIt> Lattice::get_mol_random_it(Species *sp) 
	{
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = _get_random_idxs();

		// Try all sites
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					spair = get_mol_it(Site(idxs[0][i1],idxs[1][i2],idxs[2][i3]),sp);
					if (spair.first) 
					{
						return spair;
					};
				};
			};
		};
		return std::make_pair(false,SiteIt());
	};

	/********************
	Get a free site
	********************/

	std::pair<bool,Site> Lattice::get_free_site() 
	{
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = _get_random_idxs();

		// Try all sites
		Site s;
		std::pair<bool,SiteIt> spair;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					s = Site(idxs[0][i1],idxs[1][i2],idxs[2][i3]);
					spair = get_mol_it(s);
					if (!(spair.first)) {
						return std::make_pair(true,s);
					};
				};
			};
		};

		return std::make_pair(false,s);
	};

	/********************
	Get neighbors of a site
	********************/

	std::pair<Site,std::pair<bool,SiteIt>> Lattice::get_neighbor_random(Site s)
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

	std::pair<Site,std::pair<bool,SiteIt>> Lattice::get_neighbor_random(SiteIt sit)
	{
		return get_neighbor_random(Site(sit));
	};

	std::pair<bool,Site> Lattice::get_free_neighbor_random(Site s) 
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

	std::pair<bool,Site> Lattice::get_free_neighbor_random(SiteIt sit) 
	{
		return get_free_neighbor_random(Site(sit));
	};

	/********************
	Get NN of species
	********************/

	int Lattice::get_nn(Species *sa, Species *sb)
	{
		int nn = 0;
		// Go through all sites
		Site s;
		std::pair<bool,SiteIt> spair;
		std::vector<Site> nbrs;
		for (auto i1=1; i1 <= this->_box_length; i1++) {
			for (auto i2=1; i2 <= this->_box_length; i2++) {
				for (auto i3=1; i3 <= this->_box_length; i3++) {
					s = Site(i1,i2,i3);
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

	void Lattice::write_to_file(std::string fname) 
	{
		std::ofstream f;
		f.open (fname);
		for (auto it = this->_map.begin(); it != this->_map.end(); it++) {
			for (auto it_1 = it->second.begin(); it_1 != it->second.end(); it_1++) {
				for (auto it_2 = it_1->second.begin(); it_2 != it_1->second.end(); it_2++) {
					f << it->first << " " << it_1->first << " " << it_2->first << " " << it_2->second.sp->name << "\n";
				};
			};
		};
		f.close();
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/

	/********************
	Get random indexes
	********************/

	std::map<int,std::vector<int>> Lattice::_get_random_idxs()
	{
		// Vectors of sites
		std::vector<int> x1(this->_box_length);
		std::iota(std::begin(x1), std::end(x1), 1);
		std::vector<int> x2 = x1, x3 = x1;

		// Shuffle
		std::random_shuffle(x1.begin(),x1.end());
		std::random_shuffle(x2.begin(),x2.end());
		std::random_shuffle(x3.begin(),x3.end());

		// Return
		std::map<int,std::vector<int>> m;
		m[0] = x1;
		m[1] = x2;
		m[2] = x3;
		return m;
	};

	/********************
	Get all neighbors of a site
	********************/

	std::vector<Site> Lattice::_get_all_neighbors(Site s)
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

};