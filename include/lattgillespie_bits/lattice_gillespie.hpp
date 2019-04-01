#include <string>
#include <map>

/************************************
* Namespace for lattg
************************************/

namespace lattg {

	/****************************************
	Main simulation class
	****************************************/

	class Simulation
	{
	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor/Destructor
		********************/

		Simulation(double dt, int box_length, int dim=3);
		Simulation(Simulation&& other); // movable but no copies
	    Simulation& operator=(Simulation&& other); // movable but no copies
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
		void populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, int n_steps);
		void populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, int n_steps);
		void populate_lattice(std::map<std::string,double> &h_dict, std::map<std::string,std::map<std::string,double>> &j_dict, std::map<std::string, std::map<std::string,std::map<std::string,double>>> &k_dict, std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double>>>> &q_dict, int n_steps);

        /********************
        Set directory
         ********************/

        void set_write_dir(std::string dir);
        
		/********************
		Run simulation
		********************/

		void run(int n_timesteps, bool verbose = true, bool write_counts = false, bool write_nns = false, bool write_latt = false, int write_step = 20, int write_version_no = 0, std::string dir=".");

		/********************
		Write/Read lattice
		********************/

		void write_lattice(int index, int write_version_no, std::string dir=".");
		void read_lattice(std::string fname);
	};
};
