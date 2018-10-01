#include <string>
#include <vector>
#include <map>

/************************************
* Namespace for lattg
************************************/

namespace lattg {

	// Forwards
	struct Species;
	struct Mol;

	/****************************************
	General functions
	****************************************/

	typedef std::map<int,std::map<int,std::map<int,Mol>>> lattice_map_3D;
	typedef std::map<int,std::map<int,Mol>> lattice_map_2D;
	typedef std::map<int,Mol> lattice_map_1D;

	/****************************************
	Structure to hold a lattice site iterator
	****************************************/

	struct SiteIt3D {
		lattice_map_3D::iterator it_3D;
		lattice_map_2D::iterator it_2D;
		lattice_map_1D::iterator it_1D;	

		// Constructor
		SiteIt3D();
		SiteIt3D(lattice_map_3D::iterator it_3DIn, lattice_map_2D::iterator it_2DIn, lattice_map_1D::iterator it_1DIn);	
	};
	std::ostream& operator<<(std::ostream& os, const SiteIt3D& sit);

	/****************************************
	Structure to hold a lattice site
	****************************************/

	struct Site3D {
		int x;
		int y;
		int z;	

		// Constructor
		Site3D();
		Site3D(int xIn, int yIn, int zIn);
		Site3D(SiteIt3D sit);
	};
	// Comparator
	bool operator <(const Site3D& a, const Site3D& b);
	bool operator==(const Site3D& a, const Site3D& b);
	std::ostream& operator<<(std::ostream& os, const Site3D& s);

	/****************************************
	Lattice
	****************************************/

	class Lattice
	{
	protected:

		// Internal maps
		lattice_map_3D _map;

		// Size in each dim
		int _box_length_x;
		int _box_length_y;
		int _box_length_z;

		// No dims
		int _dim;

		// Steps to search for neighbors
		std::vector<Site3D> _steps_nbrs;
		std::vector<std::pair<Site3D,Site3D>> _steps_triplet_paths;
		std::vector<std::vector<Site3D>> _steps_quartic_paths;

	private:

		// Check if a site is in the latt
		bool _check_if_in_latt(Site3D s);

		// Helpers
		void _clean_up();
		void _copy(const Lattice& other);
		void _reset();

	public:

		/********************
		Constructor/Destructor
		********************/

		Lattice(int box_length);
		Lattice(int box_length_1, int box_length_2);
		Lattice(int box_length_1, int box_length_2, int box_length_3);
		Lattice(const Lattice& other);
		Lattice(Lattice&& other);
		Lattice& operator=(const Lattice& other);
	    Lattice& operator=(Lattice&& other);
		~Lattice();

		/********************
		Clear, size
		********************/

		void clear();
		int size();

		/********************
		Make a mol
		********************/

		std::pair<bool,SiteIt3D> make_mol(Site3D s, Species *sp);
		std::pair<bool,SiteIt3D> replace_mol(Site3D s, Species *sp);
		std::pair<bool,SiteIt3D> make_mol_random(Species *sp);
		std::pair<bool,SiteIt3D> make_mol_at_empty(Site3D s, Species *sp);

		/********************
		Erase a mol
		********************/

		bool erase_mol(Site3D s);
		bool erase_mol_it(SiteIt3D sit);
		std::pair<bool,Site3D> erase_mol_random(Species *sp);

		/********************
		Get a mol
		********************/

		std::pair<bool,SiteIt3D> get_mol_it(Site3D s);
		std::pair<bool,SiteIt3D> get_mol_it(Site3D s, Species *sp);
		std::pair<bool,SiteIt3D> get_mol_random_it();
		std::pair<bool,SiteIt3D> get_mol_random_it(Species *sp);

		/********************
		Get a free site
		********************/

		std::pair<bool,Site3D> get_free_site();

		/********************
		Get neighbors of a site
		********************/

		std::pair<Site3D,std::pair<bool,SiteIt3D>> get_neighbor_random(Site3D s);
		std::pair<Site3D,std::pair<bool,SiteIt3D>> get_neighbor_random(SiteIt3D sit);
		std::pair<bool,Site3D> get_free_neighbor_random(Site3D s);
		std::pair<bool,Site3D> get_free_neighbor_random(SiteIt3D sit);

		/********************
		Get NN of species
		********************/

		int get_nn(Species *sa, Species *sb);

		/********************
		Write/Read lattice to a file
		********************/

		void write_to_file(std::string fname);
		void read_from_file(std::string fname, std::map<std::string,Species*> sp_map);

		/********************
		Sample
		********************/

		void sample(std::map<Species*,double> &h_dict,std::map<Species*, std::map<Species*,double>> &j_dict, std::map<Species*, std::map<Species*, std::map<Species*,double>>> &k_dict, std::map<Species*, std::map<Species*, std::map<Species*,std::map<Species*, double>>>> &q_dict, int n_steps);

		/********************
		Sample probabilities/propensities
		********************/

		// Sample an unnormalized probability vector
		int sample_prop_vec(std::vector<double> &props);
		
		// Sample a vector of propensities (cumulative probabilities)
		int sample_prob_vec(std::vector<double> &probs);

		/********************
		Get random indexes
		********************/

		std::map<int,std::vector<int>> get_random_idxs();

		/********************
		Get all neighbors of a site
		********************/

		std::vector<Site3D> get_all_neighbors(Site3D s);
		std::vector<std::pair<Site3D,Site3D>> get_all_triplet_considerations(Site3D s);
		std::vector<std::vector<Site3D>> get_all_quartic_considerations(Site3D s);
	};

};