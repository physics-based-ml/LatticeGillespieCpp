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

	typedef std::vector<Mol> mvec;
	typedef std::map<Site, mvec> lattice_map;

	// Random numbers
	double randD(double dMin, double dMax);
	int randI(int iMin, int iMax);

	// Print mvec
	std::ostream& operator<<(std::ostream& os, const mvec& mv);

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

	public:

		/********************
		Constructor/Destructor
		********************/

		Lattice(int box_length);
		~Lattice();

		/********************
		Erase/create/move mols
		********************/

		// Erase a mol of a certain species from a site
		void erase_mol(Site s, Species *sp);
		Site erase_mol_random(Species *sp);

		// Make a new site
		void make_mol(Site s, Species *sp);
		Site make_mol_random(Species *sp);

		/********************
		Get neighbors of a site
		********************/

		// Get neighbors
		std::vector<Site> get_neighbors(Site s);

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