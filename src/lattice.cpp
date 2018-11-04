#include "lattice.hpp"

// Other headers
#include "species.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include "math.h"
#include <random>
#include <sstream>
#include <algorithm>
#include <chrono>

/************************************
* Namespace for lattg
************************************/

namespace lattg {

	/****************************************
	Function to find all possible triplet path steps for a given dimension
	****************************************/

	std::vector<std::pair<Site3D,Site3D>> find_triplet_path_steps(int dim) {
		// All sites of nbrhood of size 2
		/*
		std::vector<Site3D> all_steps;
		Site3D s(0,0,0);
		if (dim == 1) {
			for (int iz=-2; iz<=2; iz++) {
				s.z = iz;
				all_sites.push_back(s);
			};
		} else if (dim == 2) {
			for (int iy=-2; iy<=2; iy++) {
				for (int iz=-2; iz<=2; iz++) {
					s.y = iy;
					s.z = iz;
					all_sites.push_back(s);
				};
			};
		} else if (dim == 3) {
			for (int ix=-2; ix<=2; ix++) {
				for (int iy=-2; iy<=2; iy++) {
					for (int iz=-2; iz<=2; iz++) {
						s.x = ix;
						s.y = iy;
						s.z = iz;
						all_sites.push_back(s);
					};
				};
			};
		};
		*/

		// Pick all possible triplets
		// ...

		std::vector<std::pair<Site3D,Site3D>> ret;
		ret.push_back(std::make_pair(Site3D(0,0,-2),Site3D(0,0,-1)));
		ret.push_back(std::make_pair(Site3D(0,0,-1),Site3D(0,0,1)));
		ret.push_back(std::make_pair(Site3D(0,0,1),Site3D(0,0,2)));

		return ret;
	};

	std::vector<std::vector<Site3D>> find_quartic_path_steps(int dim) {
		std::vector<std::vector<Site3D>> ret;
		ret.push_back(std::vector<Site3D>({Site3D(0,0,-3),Site3D(0,0,-2),Site3D(0,0,-1)}));
		ret.push_back(std::vector<Site3D>({Site3D(0,0,-2),Site3D(0,0,-1),Site3D(0,0,1)}));
		ret.push_back(std::vector<Site3D>({Site3D(0,0,-1),Site3D(0,0,1),Site3D(0,0,2)}));
		ret.push_back(std::vector<Site3D>({Site3D(0,0,1),Site3D(0,0,2),Site3D(0,0,3)}));
		return ret;
	};

	/****************************************
	Structure to hold a lattice site iterator
	****************************************/

	SiteIt3D::SiteIt3D() {};
	SiteIt3D::SiteIt3D(lattice_map_3D::iterator it_3DIn, lattice_map_2D::iterator it_2DIn, lattice_map_1D::iterator it_1DIn) 
	{
		this->it_3D = it_3DIn;
		this->it_2D = it_2DIn;
		this->it_1D = it_1DIn;
	};
	std::ostream& operator<<(std::ostream& os, const SiteIt3D& sit)
	{
	    return os << sit.it_3D->first << " " << sit.it_2D->first << " " << sit.it_1D->first;
	};

	/****************************************
	Struct to hold a lattice Site3D
	****************************************/

	// Constructor
	Site3D::Site3D() {};
	Site3D::Site3D(int xIn, int yIn, int zIn) {
		this->x = xIn;
		this->y = yIn;
		this->z = zIn;
	};
	Site3D::Site3D(SiteIt3D sit) {
		this->x = sit.it_3D->first;
		this->y = sit.it_2D->first;
		this->z = sit.it_1D->first;
	};

	// Comparator
	bool operator <(const Site3D& a, const Site3D& b) {
    	return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z);
	};
	bool operator==(const Site3D& a, const Site3D& b) {
		return std::tie(a.x, a.y, a.z) == std::tie(b.x, b.y, b.z);
	}; 
	std::ostream& operator<<(std::ostream& os, const Site3D& s)
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
	Lattice::Lattice(int box_length) {
		_dim = 1;
		_box_length_x = 1;
		_box_length_y = 1;
		_box_length_z = box_length;
		// Nbrs
		_steps_nbrs.push_back(Site3D(0,0,-1));
		_steps_nbrs.push_back(Site3D(0,0,1));
		// Next nbr paths
		_steps_triplet_paths = find_triplet_path_steps(_dim);
		_steps_quartic_paths = find_quartic_path_steps(_dim);
	};
	Lattice::Lattice(int box_length_1, int box_length_2) {
		_dim = 2;
		_box_length_x = 1;
		_box_length_y = box_length_1;
		_box_length_z = box_length_2;
		// Nbrs
		_steps_nbrs.push_back(Site3D(0,0,-1));
		_steps_nbrs.push_back(Site3D(0,0,1));
		_steps_nbrs.push_back(Site3D(0,-1,0));
		_steps_nbrs.push_back(Site3D(0,1,0));
	};
	Lattice::Lattice(int box_length_1, int box_length_2, int box_length_3) {
		_dim = 3;
		_box_length_x = box_length_1;
		_box_length_y = box_length_2;
		_box_length_z = box_length_3;
		// Nbrs
		_steps_nbrs.push_back(Site3D(0,0,-1));
		_steps_nbrs.push_back(Site3D(0,0,1));
		_steps_nbrs.push_back(Site3D(0,-1,0));
		_steps_nbrs.push_back(Site3D(0,1,0));
		_steps_nbrs.push_back(Site3D(-1,0,0));
		_steps_nbrs.push_back(Site3D(1,0,0));
	};
	Lattice::Lattice(const Lattice& other) {
		_copy(other);
	};
	Lattice::Lattice(Lattice&& other) {
		_copy(other);
		other._reset();
	};
	Lattice& Lattice::operator=(const Lattice& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
	};
    Lattice& Lattice::operator=(Lattice&& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
			other._reset();
		};
		return *this;
    };

	// Destructor
	Lattice::~Lattice() {
		_clean_up();
	};

	// Helpers
	void Lattice::_clean_up() {
		// Nothing...
	};
	void Lattice::_copy(const Lattice& other) {
		_box_length_x = other._box_length_x;
		_box_length_y = other._box_length_y;
		_box_length_z = other._box_length_z;
		_dim = other._dim;
		_map = other._map;
		_steps_nbrs = other._steps_nbrs;
		_steps_triplet_paths = other._steps_triplet_paths;
		_steps_quartic_paths = other._steps_quartic_paths;
	};
	void Lattice::_reset() {
		_box_length_x = 0;
		_box_length_y = 0; 
		_box_length_z = 0;
		_dim = 0;
		_map.clear();
		_steps_nbrs.clear();
		_steps_triplet_paths.clear();
		_steps_quartic_paths.clear();
	};

	/********************
	Clear, size
	********************/

	void Lattice::clear() { this->_map.clear(); };
	int Lattice::size() { 
		int ctr=0;
		for (auto it1 = this->_map.begin(); it1 != this->_map.end(); it1++) {
			for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
				for (auto it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
					ctr++;
				};
			};
		};
		return ctr; 
	};

	/********************
	Make a mol
	********************/

	std::pair<bool,SiteIt3D> Lattice::make_mol(Site3D s, Species *sp) 
	{
		// Check if the site is empty
		std::pair<bool,SiteIt3D> spair = get_mol_it(s);
		if (spair.first) {
			// Not empty
			return std::make_pair(false,SiteIt3D());
		};

		return make_mol_at_empty(s,sp);
	};

	std::pair<bool,SiteIt3D> Lattice::replace_mol(Site3D s, Species *sp) 
	{
		// Check if the site is empty
		std::pair<bool,SiteIt3D> spair = get_mol_it(s);
		if (spair.first) {
			// Not empty - erase
			erase_mol_it(spair.second);
		};

		// Now make it
		return make_mol_at_empty(s,sp);
	};

	std::pair<bool,SiteIt3D> Lattice::make_mol_random(Species *sp) 
	{
		// Get a random free site
		std::pair<bool,Site3D> spair = get_free_site();

		if (!(spair.first)) {
			// No free sites at all
			return std::make_pair(false,SiteIt3D());
		};

		// Make
		std::pair<bool,SiteIt3D> ret = make_mol_at_empty(spair.second,sp);
		return ret;
	};

	std::pair<bool,SiteIt3D> Lattice::make_mol_at_empty(Site3D s, Species *sp) {

		// Update counts
		sp->count++;

		// Make
		lattice_map_3D::iterator it_3D;
		lattice_map_2D::iterator it_2D;
		lattice_map_1D::iterator it_1D;
		it_3D = this->_map.find(s.x);
		if (it_3D == this->_map.end()) {
			auto ret = this->_map.insert(std::make_pair(s.x,lattice_map_2D()));
			it_3D = ret.first;
		};
		it_2D = it_3D->second.find(s.y);
		if (it_2D == it_3D->second.end()) {
			auto ret_1 = it_3D->second.insert(std::make_pair(s.y,lattice_map_1D()));
			it_2D = ret_1.first;
		};
		auto ret_2 = it_2D->second.insert(std::make_pair(s.z,Mol(sp)));
		it_1D = ret_2.first;

		return std::make_pair(true,SiteIt3D(it_3D,it_2D,it_1D));
	};

	/********************
	Erase a mol
	********************/

	bool Lattice::erase_mol(Site3D s) 
	{
		// Get an iterator to the site
		std::pair<bool,SiteIt3D> spair = get_mol_it(s);
		if (spair.first) {
			return erase_mol_it(spair.second);
		} else {
			return false;
		};
	};

	bool Lattice::erase_mol_it(SiteIt3D sit) 
	{
		// Update counts on species
		sit.it_1D->second.sp->count--;

		// Erase inner_2 from inner_1
		sit.it_2D->second.erase(sit.it_1D);
		if (sit.it_2D->second.size() == 0)
		{
			// Erase inner_1 from outer
			sit.it_3D->second.erase(sit.it_2D);
			if (sit.it_3D->second.size() == 0)
			{
				// Erase outer from map
				this->_map.erase(sit.it_3D);
				return true;
			};
		};
		return true; // Nonsense; it can't fail?!
	};

	std::pair<bool,Site3D> Lattice::erase_mol_random(Species *sp) 
	{
		// Get an iterator to a random site
		std::pair<bool,SiteIt3D> spair = get_mol_random_it(sp);
		if (spair.first) {
			Site3D s = Site3D(spair.second.it_3D->first,spair.second.it_2D->first,spair.second.it_1D->first);
			bool succ = erase_mol_it(spair.second);
			if (succ) {
				return std::make_pair(true,s);
			};
		};
		return std::make_pair(false,Site3D());
	};

	/********************
	Get a mol
	********************/

	std::pair<bool,SiteIt3D> Lattice::get_mol_it(Site3D s) 
	{
		auto it_3D = this->_map.find(s.x);
		if (it_3D != this->_map.end()) {
			auto it_2D = it_3D->second.find(s.y);
			if (it_2D != it_3D->second.end()) {
				auto it_1D = it_2D->second.find(s.z);
				if (it_1D != it_2D->second.end()) {
					return std::make_pair(true,SiteIt3D(it_3D,it_2D,it_1D));
				};
			};
		};
		return std::make_pair(false,SiteIt3D());	
	};

	std::pair<bool,SiteIt3D> Lattice::get_mol_it(Site3D s, Species *sp) {
		std::pair<bool,SiteIt3D> ret = get_mol_it(s);
		if (ret.first) {
			if (ret.second.it_1D->second.sp == sp) {
				return ret;
			};
		};
		return std::make_pair(false,SiteIt3D());
	};

	std::pair<bool,SiteIt3D> Lattice::get_mol_random_it() {
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = get_random_idxs();

		// Try all sites
		std::pair<bool,SiteIt3D> spair;
		for (auto i1=0; i1 < this->_box_length_x; i1++) {
			for (auto i2=0; i2 < this->_box_length_y; i2++) {
				for (auto i3=0; i3 < this->_box_length_z; i3++) {
					spair = get_mol_it(Site3D(idxs[0][i1],idxs[1][i2],idxs[2][i3]));
					if (spair.first) 
					{
						return spair;
					};
				};
			};
		};
		return std::make_pair(false,SiteIt3D());
	};

	std::pair<bool,SiteIt3D> Lattice::get_mol_random_it(Species *sp) 
	{
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = get_random_idxs();

		// Try all sites
		std::pair<bool,SiteIt3D> spair;
		for (auto i1=0; i1 < this->_box_length_x; i1++) {
			for (auto i2=0; i2 < this->_box_length_y; i2++) {
				for (auto i3=0; i3 < this->_box_length_z; i3++) {
					spair = get_mol_it(Site3D(idxs[0][i1],idxs[1][i2],idxs[2][i3]),sp);
					if (spair.first) 
					{
						return spair;
					};
				};
			};
		};
		return std::make_pair(false,SiteIt3D());
	};
	/********************
	Get a free site
	********************/

	std::pair<bool,Site3D> Lattice::get_free_site() 
	{
		// Get random indexes to search
		std::map<int,std::vector<int>> idxs = get_random_idxs();

		// Try all sites
		Site3D s;
		std::pair<bool,SiteIt3D> spair;
		for (auto i1=0; i1 < this->_box_length_x; i1++) {
			for (auto i2=0; i2 < this->_box_length_y; i2++) {
				for (auto i3=0; i3 < this->_box_length_z; i3++) {
					s = Site3D(idxs[0][i1],idxs[1][i2],idxs[2][i3]);
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

	std::pair<Site3D,std::pair<bool,SiteIt3D>> Lattice::get_neighbor_random(Site3D s)
	{
		// All neighbors
		std::vector<Site3D> nbrs = get_all_neighbors(s);

		// Shuffle
		std::random_shuffle(nbrs.begin(),nbrs.end());

		// Random choice
		auto it = nbrs.begin();
		std::advance(it,randI(0,nbrs.size()-1));

		// Check if occupied
		std::pair<bool,SiteIt3D> occ = get_mol_it(*it);
		if (occ.first) {
			// Yes
			return std::make_pair(Site3D(occ.second),std::make_pair(true,occ.second));
		} else {
			// No
			return std::make_pair(*it,std::make_pair(false,SiteIt3D()));
		};
	};

	std::pair<Site3D,std::pair<bool,SiteIt3D>> Lattice::get_neighbor_random(SiteIt3D sit)
	{
		return get_neighbor_random(Site3D(sit));
	};

	std::pair<bool,Site3D> Lattice::get_free_neighbor_random(Site3D s) 
	{
		// All allowed nbrs
		std::vector<Site3D> nbrs = get_all_neighbors(s);

		// Shuffle
		std::random_shuffle(nbrs.begin(),nbrs.end());

		// Check all neighbors
		std::pair<bool,SiteIt3D> spair;
		for (auto n: nbrs) {
			spair = get_mol_it(n);
			if (!(spair.first)) {
				return std::make_pair(true,n);
			};
		};

		return std::make_pair(false,Site3D());
	};

	std::pair<bool,Site3D> Lattice::get_free_neighbor_random(SiteIt3D sit) 
	{
		return get_free_neighbor_random(Site3D(sit));
	};

	/********************
	Get NN of species
	********************/

	int Lattice::get_nn(Species *sa, Species *sb)
	{
		int nn = 0;
		// Go through all sites
		Site3D s;
		std::pair<bool,SiteIt3D> spair;
		std::vector<Site3D> nbrs;
		for (auto i1=1; i1 <= this->_box_length_x; i1++) {
			for (auto i2=1; i2 <= this->_box_length_y; i2++) {
				for (auto i3=1; i3 <= this->_box_length_z; i3++) {
					s = Site3D(i1,i2,i3);
					spair = get_mol_it(s,sa);
					if (spair.first) {
						// This one is species A
						// Now search all neigbors
						nbrs = get_all_neighbors(s);
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
		// Check that we can open it
		if (!f.is_open()) {
			std::cerr << ">>> Error: Lattice::write_to_file <<< could not open: " << fname << " for writing" << std::endl;
			exit(EXIT_FAILURE);
		};

		for (auto it_3D = this->_map.begin(); it_3D != this->_map.end(); it_3D++) {
			for (auto it_2D = it_3D->second.begin(); it_2D != it_3D->second.end(); it_2D++) {
				for (auto it_1D = it_2D->second.begin(); it_1D != it_2D->second.end(); it_1D++) {
					if (_dim == 1) {
						f << it_1D->first << " " << it_1D->second.sp->name << "\n";
					} else if (_dim == 2) {
						f << it_2D->first << " " << it_1D->first << " " << it_1D->second.sp->name << "\n";
					} else if (_dim == 3) {
						f << it_3D->first << " " << it_2D->first << " " << it_1D->first << " " << it_1D->second.sp->name << "\n";
					};
				};
			};
		};
		f.close();
	};

	void Lattice::read_from_file(std::string fname, std::map<std::string,Species*> sp_map) 
	{
		// Clear current
		clear();

		std::ifstream f;
		f.open(fname);
		std::string site_x="",site_y="",site_z="";
		std::string species="";
		std::string line;
		std::istringstream iss;
		if (f.is_open()) { // make sure we found it
			while (getline(f,line)) {
				if (line == "") { continue; };
				iss = std::istringstream(line);
				if (_dim == 1) {
					iss >> site_z;
				} else if (_dim == 2) {
					iss >> site_y;
					iss >> site_z;
				} else if (_dim == 3) {
					iss >> site_x;
					iss >> site_y;
					iss >> site_z;
				};
			    iss >> species;
			    // Find species
			    auto sp = sp_map.find(species);
			    if (sp == sp_map.end()) {
			    	std::cerr << "Error could not find species while reading lattice" << std::endl;
			    	exit(EXIT_FAILURE);
			    };
		    	// Add
		    	std::pair<bool,SiteIt3D> ret;
		    	if (_dim == 1) {
			    	ret = make_mol(Site3D(1,1,atoi(site_z.c_str())),sp->second);
			    } else if (_dim == 2) {
			    	ret = make_mol(Site3D(1,atoi(site_y.c_str()),atoi(site_z.c_str())),sp->second);
			    } else if (_dim == 3) {
			    	ret = make_mol(Site3D(atoi(site_x.c_str()),atoi(site_y.c_str()),atoi(site_z.c_str())),sp->second);
			    };
		    	if (!ret.first) {
		    		std::cerr << "Error making mol?" << std::endl;
		    		exit(EXIT_FAILURE);
		    	};
		    	// Reset
		    	site_x=""; site_y=""; site_z="";
		    	species="";
			};
		};
		f.close();
	};

	/********************
	Sample
	********************/

	void Lattice::sample(std::map<Species*,double> &h_dict,std::map<Species*, std::map<Species*,double>> &j_dict, std::map<Species*, std::map<Species*, std::map<Species*,double>>> &k_dict, std::map<Species*, std::map<Species*, std::map<Species*,std::map<Species*,double>>>> &q_dict, int n_steps) {

		if (_dim != 1 && ( k_dict.size() > 0 || q_dict.size() > 0) ) {
			std::cerr << "ERROR! Triplet/quartic sampling for dims > 1 not yet supported" << std::endl;
			exit(EXIT_FAILURE);
		};

		// Construct a vec of all possible species
		std::vector<Species*> sp_vec;
		for (auto sp_pair: h_dict) {
			sp_vec.push_back(sp_pair.first);
		};

		// Declarations
		Site3D s;
		double energy;
		std::vector<double> props; // propensities
		std::vector<double> probs; // probs
		int i_chosen;
		std::vector<Site3D> nbrs;
		std::vector<std::pair<Site3D,Site3D>> nnbrs;
		std::vector<std::vector<Site3D>> nnnbrs;
		std::pair<bool,SiteIt3D> ret_nbr1,ret_nbr2,ret_nbr3;

		// Go over all steps
		for (int i_step=0; i_step<n_steps; i_step++) {

			// Go through all lattice sites
			for (int i=1; i<=_box_length_x; i++) {
				for (int j=1; j<=_box_length_y; j++) {
					for (int k=1; k<=_box_length_z; k++) {

						// The site
						s = Site3D(i,j,k);

						// Clear propensities
						probs.clear();
						props.clear();
						props.push_back(0.0);

						// Propensity for no spin is exp(0) = 1
						probs.push_back(1.0);
						props.push_back(1.0);
						
						// Go through all possible species this could be, calculate propensities
						for (auto sp_new: sp_vec) {

							// Bias
							energy = h_dict[sp_new];

							// NNs for J 
							nbrs = get_all_neighbors(s);
							for (auto nbr: nbrs) {
								ret_nbr1 = get_mol_it(nbr);
								if (ret_nbr1.first) {
									// Occupied
									energy += j_dict[sp_new][ret_nbr1.second.it_1D->second.sp];
								};
							};

							// Triplets for K
							if (k_dict.size() != 0) {
								nnbrs = get_all_triplet_considerations(s);
								for (auto nbr_pair: nnbrs) {
									ret_nbr1 = get_mol_it(nbr_pair.first);
									if (ret_nbr1.first) {
										ret_nbr2 = get_mol_it(nbr_pair.second);
										if (ret_nbr2.first) {
											// Both occupied
											energy += k_dict[sp_new][ret_nbr1.second.it_1D->second.sp][ret_nbr2.second.it_1D->second.sp];
										};
									};
								};
							};

							// Quartics for Q
							if (q_dict.size() != 0) {
								nnnbrs = get_all_quartic_considerations(s);
								for (auto nbrs: nnnbrs) {
									ret_nbr1 = get_mol_it(nbrs[0]);
									if (ret_nbr1.first) {
										ret_nbr2 = get_mol_it(nbrs[1]);
										if (ret_nbr2.first) {
											ret_nbr3 = get_mol_it(nbrs[2]);
											if (ret_nbr3.first) {
												// Both occupied
												energy += q_dict[sp_new][ret_nbr1.second.it_1D->second.sp][ret_nbr2.second.it_1D->second.sp][ret_nbr3.second.it_1D->second.sp];
											};
										};
									};
								};
							};

							// Append prob
							probs.push_back(exp(energy));
							props.push_back(props.back()+exp(energy));
						};

						// Normalize probs
						double tot=0.0;
						for (auto pr: probs) {
							tot += pr;
						};
						for (int i=0; i<probs.size(); i++) {
							probs[i] /= tot;
							// std::cout << probs[i] << " ";
						};
						// std::cout << std::endl;

						// Sample RV
						i_chosen = sample_prop_vec(props);

						if (i_chosen==0) {
							// Flip down (new spin = 0)
							erase_mol(s);
							// std::cout << " flipped down" << std::endl;
						} else {
							// Make the appropriate species at this site
							replace_mol(s,sp_vec[i_chosen-1]);
							// std::cout << " flipped up" << std::endl;
						};

					};
				};
			};
		};
	};

	/********************
	Sample probabilities/propensities
	********************/

	// Sample an unnormalized probability vector
	int Lattice::sample_prob_vec(std::vector<double> &probs) {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::discrete_distribution<int> dd(probs.begin(),probs.end());
		return dd(generator);
	};

	// Sample a vector of propensities (cumulative probabilities)
	int Lattice::sample_prop_vec(std::vector<double> &props) {
		// Sample RV
		double r = randD(0.0,props.back());

		// Find interval
		for (int i=0; i<props.size()-1; i++) {
			if (props[i] <= r && r <= props[i+1]) {
				return i;
			};
		};
		return 0; // never get here
	};

	/********************
	Get random indexes
	********************/

	std::map<int,std::vector<int>> Lattice::get_random_idxs()
	{
		// Vectors of sites
		std::vector<int> x1(this->_box_length_x);
		std::iota(std::begin(x1), std::end(x1), 1);
		std::vector<int> x2(this->_box_length_y);
		std::iota(std::begin(x2), std::end(x2), 1);
		std::vector<int> x3(this->_box_length_z);
		std::iota(std::begin(x3), std::end(x3), 1);

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
	Check if a site is in the lattice
	********************/

	// Check if a site is in the latt
	bool Lattice::_check_if_in_latt(Site3D s) {
		if (s.x <= this->_box_length_x && s.x >= 1 && s.y <= this->_box_length_y && s.y >= 1 && s.z <= this->_box_length_z && s.z >= 1) {
			return true;
		} else {
			return false;
		};
	};

	/********************
	Get all neighbors of a site
	********************/

	std::vector<Site3D> Lattice::get_all_neighbors(Site3D s)
	{
		std::vector<Site3D> nbrs;
		Site3D nbr = s;

		for (auto step: _steps_nbrs) {
			nbr = s;
			nbr.x += step.x;
			nbr.y += step.y;
			nbr.z += step.z;
			// Check
			if (_check_if_in_latt(nbr)) {
				nbrs.push_back(nbr);
			};
		};
		return nbrs;
	};

	std::vector<std::pair<Site3D,Site3D>> Lattice::get_all_triplet_considerations(Site3D s)
	{
		std::vector<std::pair<Site3D,Site3D>> nnbrs;
		Site3D nbr1 = s;
		Site3D nbr2 = s;

		for (auto step: _steps_triplet_paths) {
			nbr1 = s;
			nbr1.x += step.first.x;
			nbr1.y += step.first.y;
			nbr1.z += step.first.z;
			// Check first
			if (_check_if_in_latt(nbr1)) {
				nbr2 = s;
				nbr2.x += step.second.x;
				nbr2.y += step.second.y;
				nbr2.z += step.second.z;
				// Check second
				if (_check_if_in_latt(nbr2)) {
					nnbrs.push_back(std::make_pair(nbr1,nbr2));
				};
			};
		};

		return nnbrs;
	};

	std::vector<std::vector<Site3D>> Lattice::get_all_quartic_considerations(Site3D s)
	{
		std::vector<std::vector<Site3D>> nnnbrs;
		Site3D nbr1 = s;
		Site3D nbr2 = s;
		Site3D nbr3 = s;

		for (auto step: _steps_quartic_paths) {
			nbr1 = s;
			nbr1.x += step[0].x;
			nbr1.y += step[0].y;
			nbr1.z += step[0].z;
			// Check first
			if (_check_if_in_latt(nbr1)) {
				nbr2 = s;
				nbr2.x += step[1].x;
				nbr2.y += step[1].y;
				nbr2.z += step[1].z;
				// Check second
				if (_check_if_in_latt(nbr2)) {
					nbr3 = s;
					nbr3.x += step[2].x;
					nbr3.y += step[2].y;
					nbr3.z += step[2].z;
					// Check third
					if (_check_if_in_latt(nbr3)) {
						nnnbrs.push_back(std::vector<Site3D>({nbr1,nbr2,nbr3}));
					};
				};
			};
		};

		return nnnbrs;
	};



};
