// map
#ifndef MAP_h
#define MAP_h
#include <map>
#endif

// vector
#ifndef VECTOR_h
#define VECTOR_h
#include <vector>
#endif

// string
#ifndef STRING_h
#define STRING_h
#include <string>
#endif

// tuple
#ifndef TUPLE_h
#define TUPLE_h
#include <tuple>
#endif

// utility for pair
#ifndef PAIR_h
#define PAIR_h
#include <utility>
#endif

/************************************
* Namespace for Gillespie2D
************************************/

namespace Gillespie2D {

	/****************************************
	General functions
	****************************************/

	double randD(double dMin, double dMax);
	int randI(int iMin, int iMax);

	// General function to write vector to file
	void write_vector_to_file(std::string fname, std::vector<int> v);

	/****************************************
	Pair of Species
	****************************************/

	struct SpeciesPair {
		std::string s1,s2;
		SpeciesPair(std::string s1In, std::string s2In);
	};
	// Comparator
	bool operator <(const SpeciesPair& a, const SpeciesPair& b);

	/****************************************
	Reaction
	****************************************/

	class Reaction
	{
	private:
		
		// Name
		std::string _name;
		// Products
		std::vector<std::string> _p;
		// Conserved reactants and products
		std::vector<std::string> _r_c;
		std::vector<std::string> _p_c;

	public:

		// Constructor destructor
		Reaction(std::string name);
		Reaction(std::string name, std::string p);
		Reaction(std::string name, std::string p1, std::string p2);
		~Reaction();		

		// Get the name
		const std::string name();

		// Get products
		const std::vector<std::string> get_p();

		// Get/set/clear conserved
		std::vector<std::string> get_r_cons();
		std::vector<std::string> get_p_cons();
		void set_r_cons(std::string s);
		void set_p_cons(std::string s);
		void clear_r_cons();
		void clear_p_cons();
	};

	/****************************************
	Unimolecular reaction
	****************************************/

	class UniReaction : public Reaction
	{
	private:
		
		// Prob
		double _kr;

		// Reactant
		std::string _r;

	public:

		// Constructor destructor
		UniReaction(std::string name, double kr, std::string r);
		UniReaction(std::string name, double kr, std::string r, std::string p);
		UniReaction(std::string name, double kr, std::string r, std::string p1, std::string p2);
		~UniReaction();

		// Get reactants/products
		const std::string get_r();

		// Get the rate
		const double get_rate();

		// Get the name
		const std::string name();

		// Get products
		const std::vector<std::string> get_p();

		// Get/set/clear conserved
		std::vector<std::string> get_r_cons();
		std::vector<std::string> get_p_cons();
		void set_r_cons(std::string s);
		void set_p_cons(std::string s);
		void clear_r_cons();
		void clear_p_cons();
	};

	/****************************************
	Bimolecular reaction
	****************************************/

	class BiReaction : public Reaction
	{
	private:
		
		// Prob
		double _prob;
		// Reactants
		SpeciesPair _r;

	public:

		// Constructor destructor
		BiReaction(std::string name, double prob, std::string r1, std::string r2);
		BiReaction(std::string name, double prob, std::string r1, std::string r2, std::string p);
		BiReaction(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2);
		~BiReaction();

		// Get reactants/products
		const SpeciesPair get_r();

		// Get the rate
		const double get_prob();	

		// Get the name
		const std::string name();

		// Get products
		const std::vector<std::string> get_p();

		// Get/set/clear conserved
		std::vector<std::string> get_r_cons();
		std::vector<std::string> get_p_cons();
		void set_r_cons(std::string s);
		void set_p_cons(std::string s);
		void clear_r_cons();
		void clear_p_cons();
	};

	/****************************************
	Structure to hold a lattice Site
	****************************************/

	struct Site {
		int x;
		int y;

		// Constructor
		Site(int xIn, int yIn);	
	};
	// Comparator
	bool operator <(const Site& a, const Site& b);
	bool operator==(const Site& a, const Site& b);
	std::ostream& operator<<(std::ostream& os, const Site& s);

	/****************************************
	Lattice
	****************************************/

	class Lattice
	{
	private:

		typedef std::map<Site, std::string> lattice_map;

		// Internal unorderd set
		lattice_map _map;

		// Size
		int _box_length;

		// Counts of different species on the lattice
		std::map<std::string,int> _counts;

	public:
		// Constructors
		Lattice(int box_length);

		// Destructor
		~Lattice();

		// Erase a site
		void erase_site(Site s);

		// Make a new site
		void make_site(Site s, std::string species);

		// Get a list of all occupied sites
		std::vector<Site> get_sites();

		// Get a random free site
		std::pair<bool,Site> get_random_free_site();

		// Get an occupied Site, if it exists
		std::pair<bool,std::string> get_site(Site s);

		// Find a site by species
		std::pair<bool,Site> find_site(std::string species);

		// Get neighbors
		std::vector<Site> get_neighbors(Site s);
		std::vector<Site> get_neighbors_empty(Site s);
		std::vector<Site> get_neighbors_occupied(Site s);

		// Get number of particles of a certain species, neighbors
		std::pair<bool,int> get_n_particles_safe(std::string species);
		std::map<std::string,int> get_n_particles_all();
		int get_n_particles(std::string s_name);
		int get_n_particles_total();
		int get_nn();
		int get_nn(SpeciesPair s_pair);

		// Populate lattice randomly
		void populate_lattice(std::map<std::string,int> counts);

		// Function to populate lattice according to coupling strengths
		void populate_lattice(std::map<std::string,double> &h_dict,std::map<SpeciesPair,double> &j_dict, int n_annealing_steps, bool write_lattice = false, bool write_statistics = false);

		// Write lattice to a file
		void write_to_file(std::string fname);
	};

	/****************************************
	Main simulation class
	****************************************/

	class Simulation
	{
	private:

		// The lattice
		Lattice _lattice;

		// Reaction dictionaries
		std::map<SpeciesPair,std::vector<BiReaction>> _bimol_rxn_dict;
		std::map<std::string,std::vector<UniReaction>> _unimol_rxn_dict;

		// Diffuse all the mols and do bimol reactions
		void _diffuse_mols();

		// Conserved species
		std::vector<std::string> _conserved_quantities;

		// Current time
		double _t;

		// Timestep
		double _dt;

		// Time and value of the next unimolecular reaction
		double _t_uni_next;
		UniReaction *_uni_next; // Null indicates there is none scheduled

		// Schedule the next uni reaction
		void _schedule_uni();

		// Do the scheduled reaction
		void _do_uni_rxn();

		// Counter reactions
		std::map<std::string,std::vector<int>> _rxn_count;
		// Counter collisions
		std::map<std::string,std::vector<int>> _coll_count;
		// Counter counts
		std::map<std::string,std::vector<int>> _n_count;

		// Enforce conservation laws upon reactions
		void _enforce_conservation_in_rxns();

	public:

		// Constructor
		Simulation(Lattice l0, double dt);
		// Destructor
		~Simulation();

		// Add a reaction
		void add_reaction(UniReaction rxn);
		void add_reaction(BiReaction rxn);

		// Set conserved quantities
		void set_conserved(std::string species);

		// Run simulation
		void run(int n_timesteps, bool verbose = false, bool write_statistics = false);

		// Write lattice
		void write_lattice(int index);

		// Get reaction, coll counts
		std::map<std::string,std::vector<int>> get_rxn_count();
		std::map<std::string,std::vector<int>> get_coll_count();
	};

};