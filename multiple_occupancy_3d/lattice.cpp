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
#include <numeric> // iota
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

	/********************
	Print mvec
	********************/
	
	std::ostream& operator<<(std::ostream& os, const mvec& mv) 
	{
		for (auto m : mv) { os << m.sp->name << " "; };
		return os;
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

	void Lattice::erase_mol(Site s, Species *sp) 
	{
		// Update counts
		sp->count--;

		// Get the site
		lattice_map::iterator ret;
		ret = this->_map.find(s);
		if (ret != this->_map.end()) {
			// Find a mol of the given species
			mvec::iterator it = ret->second.begin();
			while (it != ret->second.end()) {
				if (it->sp == sp) {
					// Erase
					ret->second.erase(it);

					// Possibly erase entry of site
					if (ret->second.size() == 0) {
						this->_map.erase(ret);
					};
				};
				it++;
			};
		};
	};
 
	Site Lattice::erase_mol_random(Species *sp) {
		lattice_map::iterator itmap;

		// Update counts
		sp->count--;

		// Vectors of sites
		std::vector<int> x1(this->_box_length);
		std::iota(std::begin(x1), std::end(x1), 1);
		std::vector<int> x2 = x1, x3 = x1;

		// Shuffle
		std::random_shuffle(x1.begin(),x1.end());
		std::random_shuffle(x2.begin(),x2.end());
		std::random_shuffle(x3.begin(),x3.end());

		// Try all sites
		for (auto i1=0; i1 < this->_box_length; i1++) {
			for (auto i2=0; i2 < this->_box_length; i2++) {
				for (auto i3=0; i3 < this->_box_length; i3++) {
					itmap = this->_map.find(Site(x1[i1],x2[i2],x3[i3]));
					if (itmap != this->_map.end()) 
					{
						// Check if a mol of the desired species exists
						for (auto it = itmap->second.begin(); it != itmap->second.end(); it++) {
							if (it->sp == sp) {
								// Yes; erase
								itmap->second.erase(it);

								// Possibly erase entry of site
								if (itmap->second.size() == 0) {
									this->_map.erase(itmap);
								};

								// Stop
								return Site(x1[i1],x2[i2],x3[i3]);
							};
						};
					};
				};
			};
		};

		return Site(0,0,0);
	};

	void Lattice::make_mol(Site s, Species *sp) 
	{
		// Update counts
		sp->count++;

		// Insert
		this->_map[s].push_back(Mol(sp));
		return;
	};

	Site Lattice::make_mol_random(Species *sp) 
	{
		// Update counts
		sp->count++;

		Site s(randI(1,this->_box_length),randI(1,this->_box_length),randI(1,this->_box_length));

		// Insert
		this->_map[s].push_back(Mol(sp));
		return s;
	};

	/********************
	Get neigbors of a site
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
		lattice_map::iterator ite = this->_map.end();
		ite--;
		for (auto it = this->_map.begin(); it != this->_map.end(); it++) {
			for (auto m : it->second) {
				f << it->first.x << " " << it->first.y << " " << it->first.z << " " << m.sp->name << "\n";
			};
		};
		f.close();
	};

	/********************
	Do a uni reaction
	********************/

	void Lattice::do_uni_rxn(UniReaction *rxn) {

		// Remove reactant
		Site s = erase_mol_random(rxn->r);

		// Conserve reactants
		if (rxn->r->conserved) {
			make_mol_random(rxn->r);
		};

		// Place products (append to the end of vector)
		for (auto p: rxn->p) {
			make_mol(s, p);

			// Conserve products
			if (p->conserved) {
				std::cout << p->name << " " << p->count;
				erase_mol_random(p);
				std::cout << " " << p->count << std::endl;
			};
		};
	};

	/********************
	Check for bimol reactions inside each voxel
	********************/

	void Lattice::check_rxns_in_voxels()
	{
		std::pair<bool,BiReaction*> ret;
		BiReaction *rxn;
		int i1, i2;
		mvec::iterator it;
		bool rxn_occ = false;

		// Fuck it, get a list of all sites
		std::vector<Site> slist;
		for (auto itm : this->_map) { 
			slist.push_back(itm.first); 
		};

		// Go over all lattice sites
		for (auto s: slist) 
		{
			// Iterate over all pairs of mols
			i1=0;
			while (i1 < this->_map[s].size())
			{
				i2=0;
				rxn_occ = false;
				while (i2 < i1 && !rxn_occ)
				{

					// Check for collision
					ret = this->_map[s][i1].check_bi_rxns_mol(this->_map[s][i2]);

					if (ret.first) {

						// The reaction succeeded
						rxn = ret.second;

						// Both i1 > i2 are consumed

						// Remove i1 first
						it = this->_map[s].begin();
						std::advance(it,i1);
						this->_map[s].erase(it);
						// Then i2
						it = this->_map[s].begin();
						std::advance(it,i2);
						this->_map[s].erase(it);

						// Update counts
						rxn->r1->count--;
						rxn->r2->count--;

						// Conserve reactants
						if (rxn->r1->conserved) {
							make_mol_random(rxn->r1);
						};
						if (rxn->r2->conserved) {
							make_mol_random(rxn->r2);
						};
						
						// Place products (append to the end of vector)
						for (auto p: rxn->p) {
							make_mol(s, p);

							// Conserve products
							if (p->conserved) {
								erase_mol_random(p);
							};
						};

						// Count reaction
						rxn->count++;

						// i1 --
						i1--;

						// Restart inner look
						rxn_occ = true;

					} else {
						i2++;
					};
				};
				i1++;
			};
		};

		// Remove empty sites
		auto itm = this->_map.begin();
		while (itm != this->_map.end()) 
		{
			// Check if this lattice site is now empty
			if (itm->second.size() == 0) {
				// Delete it
				itm = this->_map.erase(itm);
			} else {
				itm++;
			};
		};
	};

	/********************
	Diffuse all the mols and do bimol reactions
	********************/

	void Lattice::diffuse_mols() 
	{
		std::vector<Site>::iterator itn;
		std::vector<Site> nbrs;
		mvec::iterator itmol;

		// First check the reactions in the voxels of the current config
		check_rxns_in_voxels();

		// Ensure all the movement flags are off to begin
		lattice_map::iterator itm = this->_map.begin();
		while (itm != this->_map.end())
		{
			itmol = itm->second.begin();
			while (itmol != itm->second.end()) {
				itmol->has_moved = false;
				itmol++;
			};
			itm++;
		};

		// Now everyone moves
		itm = this->_map.begin();
		while (itm != this->_map.end())
		{
			// Get all the neighbors of this site
			nbrs = get_neighbors(itm->first);

			// All mols to move
			itmol = itm->second.begin();
			while (itmol != itm->second.end()) {

				// Already moved?
				if (itmol->has_moved) {
					// Yes
					itmol++;
				} else {
					// No
					
					// Mark that youve moved!
					itmol->has_moved = true;

					// Pick a random site to move to
					itn = nbrs.begin();
					std::advance(itn, randI(0,nbrs.size()-1));

					// std::cout << itm->first << " to " << *itn << std::endl;

					// Insert into new map
					this->_map[*itn].push_back(*itmol);

					// Remove old
					itmol = itm->second.erase(itmol);					
				};
			};

			// Possibly remove the old lattice site completely if empty
			if (itm->second.size() == 0) {
				itm = this->_map.erase(itm);
			} else {
				// Next site
				itm++;
			};
		};

		// Check reactions again!
		check_rxns_in_voxels();
	};

	/****************************************
	Lattice - PRIVATE
	****************************************/


};