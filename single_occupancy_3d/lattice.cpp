#ifndef LATTICE_h
#define LATTICE_h
#include "lattice.hpp"
#endif

#ifndef IOSTREAM_h
#define IOSTREAM_h
#include <iostream>
#endif

#ifndef FSTREAM_h
#define FSTREAM_h
#include <fstream>
#endif

#ifndef NUMERIC_h
#define NUMERIC_h
#include <numeric>
#endif


/************************************
* Namespace for Gillespie3D
************************************/

namespace Gillespie3D {

	/****************************************
	Struct to hold a lattice Site
	****************************************/

	// Constructor
	Site::Site(int xIn, int yIn, int zIn) {
		this->x = xIn;
		this->y = yIn;
		this->z = zIn;
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
	General functions
	****************************************/

	/********************
	Random numbers
	********************/

	double randD(double dMin, double dMax)
	{
	    return dMin + ((double)rand() / RAND_MAX) * (dMax - dMin);
	};
	int randI(int iMin, int iMax)
	{
		return iMin + rand() % (iMax - iMin + 1);
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
	Erase/create/move mols
	********************/

	void Lattice::erase_mol(Site s) 
	{
		// Get the site
		lattice_map::iterator itmap = this->_map.find(s);
		if (itmap != this->_map.end()) {
			// Update the count
			itmap->second.sp->count--;

			// Erase
			this->_map.erase(itmap);
		};
	};
 
	std::pair<bool,Site> Lattice::erase_mol_random(Species *sp) {
		lattice_map::iterator itmap;

		// Try all sites
		for (auto i1=1; i1 <= this->_box_length; i1++) {
			for (auto i2=1; i2 <= this->_box_length; i2++) {
				for (auto i3=1; i3 <= this->_box_length; i3++) {
					itmap = this->_map.find(Site(i1,i2,i3));
					if (itmap != this->_map.end()) 
					{
						if (itmap->second.sp == sp) {
							// Update counts
							sp->count--;

							// Erase
							this->_map.erase(itmap);
							return std::make_pair(true,Site(i1,i2,i3));
						};
					};
				};
			};
		};

		return std::make_pair(false,Site(0,0,0));
	};

	bool Lattice::make_mol(Site s, Species *sp) 
	{
		// Check if the site is empty
		auto it = this->_map.find(s);
		if (it != this->_map.end()) {
			return false;
		};

		// Update counts
		sp->count++;

		// Insert
		// this->_map[s] = Mol(sp);
		this->_map.insert(std::make_pair(s,Mol(sp)));
		return true;
	};

	std::pair<bool,Site> Lattice::make_mol_random(Species *sp) 
	{
		// Vectors of sites
		std::vector<int> x1(this->_box_length);
		std::iota(std::begin(x1), std::end(x1), 1);
		std::vector<int> x2 = x1, x3 = x1;

		// Shuffle
		std::random_shuffle(x1.begin(),x1.end());
		std::random_shuffle(x2.begin(),x2.end());
		std::random_shuffle(x3.begin(),x3.end());

		// Find a random free site
		lattice_map::iterator itmap;
		int ctr=0;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					ctr+=1;
					itmap = this->_map.find(Site(x1[i1],x2[i2],x3[i3]));
					if (itmap == this->_map.end()) 
					{
						// Empty
						// Update counts
						sp->count++;

						// Insert
						this->_map.insert(std::make_pair(Site(x1[i1],x2[i2],x3[i3]),Mol(sp)));

						return std::make_pair(true,Site(x1[i1],x2[i2],x3[i3]));
					};
				};
			};
		};

		return std::make_pair(false,Site(0,0,0));
	};
	
	/********************
	Populate a lattice
	********************/

	void Lattice::populate_lattice(std::map<Species*,int> counts) {
		std::pair<bool,Site> ret(true,Site(0,0,0));

	    // Go through all species
	    for (auto it: counts) {
	    	// All counts of this species
	    	for (auto i=0; i<it.second; i++) 
	    	{
	    		// Make
		    	make_mol_random(it.first);
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
		for (auto it = this->_map.begin(); it != this->_map.end(); it++) {
			f << it->first.x << " " << it->first.y << " " << it->first.z << " " << it->second.sp->name << "\n";
		};
		f.close();
	};

	/********************
	Do a uni reaction
	********************/

	void Lattice::do_uni_rxn(UniReaction *rxn) {

		// Declarations
		std::pair<bool,lattice_map::iterator> find_site_pair;
		lattice_map::iterator itmap;
		Site s(0,0,0);
		std::pair<bool,Site> check_nbr_pair(true,s);

		// Try to do the reaction at several sites
		// Failure can arise if there is not enough room for the products
		int ctr = 0;
		while (ctr < 20) {

			// Grab a random site
			find_site_pair = _get_random_site_occ(rxn->r);
			if (!(find_site_pair.first)) {
				// No sites with this species exist; stop
				return;
			};
			itmap = find_site_pair.second;
			s = itmap->first;

			// Check if there is room for the products
			if (rxn->p.size() == 2) {
				check_nbr_pair = _get_random_neighbor_free(s);
				if (!(check_nbr_pair.first)) {
					// Not enough room for the two products; try again
					ctr += 1;
					continue;
				};
			};

			// Remove the reactant
			this->_map.erase(itmap);
			rxn->r->count--;

			// Conserve reactants
			if (rxn->r->conserved) {
				make_mol_random(rxn->r);
			};

			// Place products, if needed at the neighbor site
			if (rxn->p.size() == 1) {
				make_mol(s, rxn->p[0]);
			} else if (rxn->p.size() == 2) {
				make_mol(s, rxn->p[0]);
				make_mol(check_nbr_pair.second, rxn->p[1]);
			};

			// Conserve products
			for (auto p: rxn->p)
			{
				if (p->conserved) {
					erase_mol_random(p);
				};
			};

			// Sucess
			return;
		};
	};

	/********************
	Diffuse all the mols and do bimol reactions
	********************/

	void Lattice::diffuse_mols() 
	{
		// Copy the old map
		lattice_map todo = this->_map;

		// Clear the current map
		this->_map.clear();

		// Go over all mols to move
		lattice_map::iterator itmove,itfind;
		Mol *mOld,*mColl;
		Site sOld(0,0,0),sNew(0,0,0);
		bool coll_map,coll_todo;
		std::pair<bool,BiReaction*> check_coll_pair;
		BiReaction *rxn;
		while (todo.size() > 0) {

			// Random element
			itmove = todo.begin();
			std::advance(itmove,randI(0,todo.size()-1));
			mOld = &(itmove->second);
			sOld = itmove->first;

			// Move
			sNew = _get_random_neighbor(sOld);

			// Check collisions
			itfind = todo.find(sNew);
			coll_todo = false;
			coll_map = false;
			if (itfind != todo.end()) {
				coll_todo = true;
				mColl = &(itfind->second);
			} else {
				itfind = this->_map.find(sNew);
				if (itfind != this->_map.end()) {
					coll_map = true;
					mColl = &(itfind->second);
				};
			};

			// No collision
			if (!coll_todo && !coll_map) {
				// Move
				this->_map.insert(std::make_pair(sNew,*mOld));
				// Continue
				todo.erase(itmove);
				continue;
			};

			// Collision; check if reaction occurs
			check_coll_pair = mOld->check_bi_rxns_mol(mColl);

			if (!(check_coll_pair.first)) {
				// No reaction & don't move
				this->_map.insert(std::make_pair(sOld,*mOld));
				// Continue
				todo.erase(itmove);
				continue;
			};

			// Reaction occurred
			rxn = check_coll_pair.second;
			rxn->count++;

			// Remove the reactants
			todo.erase(itmove);
			if (coll_todo) { 
				// The itfind iterator has been invalidated :(
				// Re-do
				itfind = todo.find(sNew);
				todo.erase(itfind);
			};
			if (coll_map) { this->_map.erase(itfind); };

			// Update counts on reactants
			rxn->r1->count--;
			rxn->r2->count--;

			// Conserve reactants
			if (rxn->r1->conserved) {
				make_mol_random(rxn->r1);
			};
			if (rxn->r2->conserved) {
				make_mol_random(rxn->r2);
			};

			// Place products, if needed at the old site
			if (rxn->p.size() == 1) {
				make_mol(sNew, rxn->p[0]);
			} else if (rxn->p.size() == 2) {
				make_mol(sNew, rxn->p[0]);
				make_mol(sOld, rxn->p[1]);
			};

			// Conserve products
			for (auto p: rxn->p)
			{
				if (p->conserved) {
					erase_mol_random(p);
				};
			};

			// Finished
		};
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/

	/********************
	Get an iterator to a random site
	********************/

	std::pair<bool,lattice_map::iterator> Lattice::_get_random_site_occ(Species *sp)
	{
		// Vectors of sites
		std::vector<int> x1(this->_box_length);
		std::iota(std::begin(x1), std::end(x1), 1);
		std::vector<int> x2 = x1, x3 = x1;

		// Shuffle
		std::random_shuffle(x1.begin(),x1.end());
		std::random_shuffle(x2.begin(),x2.end());
		std::random_shuffle(x3.begin(),x3.end());

		// Try all sites
		lattice_map::iterator itmap;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					itmap = this->_map.find(Site(x1[i1],x2[i2],x3[i3]));
					if (itmap != this->_map.end()) 
					{
						// Check species
						if (itmap->second.sp == sp) {
							return std::make_pair(true,itmap);
						};
					};
				};
			};
		};
		return std::make_pair(false,itmap);
	};

	std::pair<bool,Site> Lattice::_get_random_site_free()
	{
		// Vectors of sites
		std::vector<int> x1(this->_box_length);
		std::iota(std::begin(x1), std::end(x1), 1);
		std::vector<int> x2 = x1, x3 = x1;

		// Shuffle
		std::random_shuffle(x1.begin(),x1.end());
		std::random_shuffle(x2.begin(),x2.end());
		std::random_shuffle(x3.begin(),x3.end());

		// Try all sites
		lattice_map::iterator itmap;
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					itmap = this->_map.find(Site(x1[i1],x2[i2],x3[i3]));
					if (itmap == this->_map.end()) 
					{
						return std::make_pair(true,itmap->first);
					};
				};
			};
		};
		return std::make_pair(false,Site(0,0,0));
	};

	/********************
	Get neighbors of a site
	********************/

	Site Lattice::_get_random_neighbor(Site s)
	{
		// All neighbors
		std::vector<Site> nbrs = _get_all_neighbors(s);

		// Shuffle
		std::random_shuffle(nbrs.begin(),nbrs.end());

		// Random choice
		auto it = nbrs.begin();
		std::advance(it,randI(0,nbrs.size()-1));

		return *it;
	};

	std::pair<bool,Site> Lattice::_get_random_neighbor_free(Site s)
	{
		// All allowed nbrs
		std::vector<Site> nbrs = _get_all_neighbors(s);

		// Shuffle
		std::random_shuffle(nbrs.begin(),nbrs.end());

		lattice_map::iterator itmap;
		for (auto n: nbrs) {
			itmap = this->_map.find(n);
			if (itmap == this->_map.end())
			{
				return std::make_pair(true,n);
			};
		};

		return std::make_pair(false,Site(0,0,0));
	};

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