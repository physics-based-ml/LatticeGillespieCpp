// string
#ifndef STRING_h
#define STRING_h
#include <string>
#endif

// vector
#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

// map
#ifndef MAP_h
#define MAP_h
#include <map>
#endif

// Other Gillespie3D

#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
#endif

/************************************
* Namespace for Gillespie3D
************************************/

namespace Gillespie3D {

	/****************************************
	Structure to hold a lattice Site
	****************************************/

	struct Site {
		int x;
		int y;
		int z;	

		// Constructor
		Site(int xIn, int yIn, int zIn);	
	};
	// Comparator
	bool operator <(const Site& a, const Site& b);
	bool operator==(const Site& a, const Site& b);
	std::ostream& operator<<(std::ostream& os, const Site& s);

	/****************************************
	General functions
	****************************************/

	typedef std::map<Site, Mol> lattice_map;

	// Random numbers
	double randD(double dMin, double dMax);
	int randI(int iMin, int iMax);

	/****************************************
	Lattice
	****************************************/

	class Lattice
	{
	private:

		// Internal map
		lattice_map _map;

		// Size
		int _box_length;

		/********************
		Get an iterator to a random site
		********************/

		std::pair<bool,lattice_map::iterator> _get_random_site_occ(Species *sp);
		std::pair<bool,Site> _get_random_site_free();

		/********************
		Get neighbors of a site
		********************/

		Site _get_random_neighbor(Site s);
		std::pair<bool,Site> _get_random_neighbor_free(Site s);
		std::vector<Site> _get_all_neighbors(Site s);

	public:

		/********************
		Constructor/Destructor
		********************/

		Lattice(int box_length);
		~Lattice();

		/********************
		Erase/create/move mols
		********************/

		// Erase a mol from a site
		void erase_mol(Site s);
		std::pair<bool,Site> erase_mol_random(Species *sp);

		// Make a new site
		bool make_mol(Site s, Species *sp);
		std::pair<bool,Site> make_mol_random(Species *sp);

		/********************
		Populate a lattice
		********************/

		// Populate lattice randomly
		void populate_lattice(std::map<Species*,int> counts);

		/********************
		Write lattice to a file
		********************/

		// Write lattice to a file
		void write_to_file(std::string fname);

		/********************
		Do a uni reaction
		********************/

		void do_uni_rxn(UniReaction *rxn);

		/********************
		Check for bimol reactions inside each voxel
		********************/

		void check_rxns_in_voxels();

		/********************
		Diffuse all the mols and do bimol reactions
		********************/

		void diffuse_mols();

	};

};