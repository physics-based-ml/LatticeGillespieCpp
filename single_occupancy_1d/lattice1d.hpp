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

// Other Gillespie1D

#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
#endif

/************************************
* Namespace for Gillespie1D
************************************/

namespace Gillespie1D {

	/****************************************
	General functions
	****************************************/

	typedef std::map<int,Mol> lattice_map;

	/****************************************
	Structure to hold a lattice site iterator
	****************************************/

	struct SiteIt {
		lattice_map::iterator it;

		// Constructor
		SiteIt();
		SiteIt(lattice_map::iterator itIn);	
	};
	std::ostream& operator<<(std::ostream& os, const SiteIt& sit);

	/****************************************
	Structure to hold a lattice site
	****************************************/

	struct Site {
		int x;

		// Constructor
		Site();
		Site(int xIn);
		Site(SiteIt sit);
	};
	// Comparator
	bool operator <(const Site& a, const Site& b);
	bool operator==(const Site& a, const Site& b);
	std::ostream& operator<<(std::ostream& os, const Site& s);

	/****************************************
	Lattice1D
	****************************************/

	class Lattice1D
	{
	private:

		// Internal maps
		lattice_map _map;

		// Size
		int _box_length;

		/********************
		Get random indexes
		********************/

		std::vector<int> _get_random_idxs();

		/********************
		Get all neighbors of a site
		********************/

		std::vector<Site> _get_all_neighbors(Site s);

	public:

		/********************
		Constructor/Destructor
		********************/

		Lattice1D(int box_length);
		~Lattice1D();

		/********************
		Clear, size
		********************/

		void clear();
		int size();

		/********************
		Make a mol
		********************/

		std::pair<bool,SiteIt> make_mol(Site s, Species *sp);
		std::pair<bool,SiteIt> make_mol_random(Species *sp);

		/********************
		Erase a mol
		********************/

		bool erase_mol(Site s);
		bool erase_mol_it(SiteIt sit);
		std::pair<bool,Site> erase_mol_random(Species *sp);

		/********************
		Get a mol
		********************/

		std::pair<bool,SiteIt> get_mol_it(Site s);
		std::pair<bool,SiteIt> get_mol_it(Site s, Species *sp);
		std::pair<bool,SiteIt> get_mol_random_it();
		std::pair<bool,SiteIt> get_mol_random_it(Species *sp);

		/********************
		Get a free site
		********************/

		std::pair<bool,Site> get_free_site();

		/********************
		Get neighbors of a site
		********************/

		std::pair<Site,std::pair<bool,SiteIt>> get_neighbor_random(Site s);
		std::pair<Site,std::pair<bool,SiteIt>> get_neighbor_random(SiteIt sit);
		std::pair<bool,Site> get_free_neighbor_random(Site s);
		std::pair<bool,Site> get_free_neighbor_random(SiteIt sit);

		/********************
		Get NN of species
		********************/

		int get_nn(Species *sa, Species *sb);

		/********************
		Write lattice to a file
		********************/

		void write_to_file(std::string fname);

		/********************
		Anneal
		********************/

		void anneal(std::map<Species*,double> &h_dict,std::map<Species*,std::map<Species*,double>> &j_dict, int n_steps);
	};

};