// vector
#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

// list
#ifndef LIST_h
#define LIST_h
#include <list>
#endif

// string
#ifndef STRING_h
#define STRING_h
#include <string>
#endif

// Other Gillespie3D

#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
#endif

#ifndef REACTIONS_h
#define REACTIONS_h
#include "reactions.hpp"
#endif

#ifndef LATTICE_h
#define LATTICE_h
#include "lattice.hpp"
#endif

/************************************
* Namespace for Gillespie3D
************************************/

namespace Gillespie3D {

	/****************************************
	General functions
	****************************************/

	// Write vector to file
	void write_vector_to_file(std::string fname, std::vector<int> v);

	// Print mvec
	std::ostream& operator<<(std::ostream& os, const std::list<Species>& vs);

	/****************************************
	Main simulation class
	****************************************/

	class Simulation
	{
	private:

		// The lattice
		Lattice _lattice;

		// List of species
		std::list<Species> _species;

		// List of bimol rxns
		std::list<BiReaction> _bi_rxns;

		// List of unimol rxns
		std::list<UniReaction> _uni_rxns;

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

	public:

		/********************
		Constructor/Destructor
		********************/

		Simulation(double dt, int box_length);
		~Simulation();

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

		/********************
		Run simulation
		********************/

		void run(int n_timesteps, bool verbose = false, bool write_statistics = false);

		/********************
		Write lattice
		********************/

		void write_lattice(int index);
	};

};