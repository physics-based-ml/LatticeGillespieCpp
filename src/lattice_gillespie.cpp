#include "../include/lattgillespie_bits/lattice_gillespie.hpp"

// Other headers
#include "lattice.hpp"
#include "species.hpp"
#include "reactions.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <math.h>
#include <fstream>
#include <vector>

// Diagnostic flags
#define DIAG_DIFFUSE 0

/************************************
* Namespace for lattg
************************************/

namespace lattg {

	/****************************************
	General functions
	****************************************/

	// Write vector to file
	void write_vector_to_file(std::string fname, std::vector<int> v)
	{
		std::ofstream f;
		f.open (fname);
		for (auto i : v) {
			f << i << "\n";
		};
		f.close();
	};

	// Print mvec
	std::ostream& operator<<(std::ostream& os, const std::vector<Species*>& vs)
	{
		for (auto v : vs) { os << v->name << ": " << v->count << " "; };
		return os;
	};

	/****************************************
	Main simulation IMPLEMENTATION HEADER
	****************************************/
    
	class Simulation::Impl
	{
	private:

		// The dimensionality of the lattice
		int _dim;

		// The lattice
		Lattice *_lattice;

		// Vector of species
		std::vector<Species*> _species;

		// Vector of bimol rxns
		std::vector<BiReaction*> _bi_rxns;

		// Vector of unimol rxns
		std::vector<UniReaction*> _uni_rxns;

		// Box length
		int _box_length;

		// Current time
		double _t;
		int _t_step;

		// Timestep
		double _dt;

        // Writing dir
        std::string _dir_write;
        
		// Time and value of the next unimolecular reaction
		double _t_uni_next;
		UniReaction *_uni_next; // Null indicates there is none scheduled

		/********************
		Find a species by name
		********************/

		Species* _find_species(std::string name) const;

		/********************
		Schedule the next uni reaction
		********************/

		void _schedule_uni();

		// Constructor helpers
		void _clean_up();
        void _move(Impl &other);

	public:

		/********************
		Constructor/Destructor
		********************/

		Impl(double dt, int box_length, int dim);
		Impl(Impl&& other);
	    Impl& operator=(Impl&& other);
		~Impl();

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

        std::vector<std::string> _vec_from_dict_keys(const dict1 &dict) const;
        sdict1* _convert_dict_to_sdict(const dict1 &h_dict) const;
        sdict2* _convert_dict_to_sdict(const dict2 &j_dict) const;
        sdict3* _convert_dict_to_sdict(const dict3 &k_dict) const;
        sdict4* _convert_dict_to_sdict(const dict4 &q_dict) const;

        void populate_random(std::vector<std::string> sps);
		void populate_lattice(std::map<std::string,int> counts);
        void populate_lattice_1d(const dict1 &h_dict, int n_steps);
		void populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, int n_steps);
		void populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, const dict3 &k_dict, int n_steps);
		void populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, const dict3 &k_dict, const dict4 &q_dict, int n_steps);

		/********************
		Do a uni reaction
		********************/

		void do_uni_rxn(UniReaction *rxn);

		/********************
		Diffuse all the mols and do bimol reactions
		********************/

		void diffuse_mols(bool periodic);

        /********************
         Set directory
         ********************/
        
        void set_write_dir(std::string dir);
        
		/********************
		Run simulation
		********************/

		void run(int n_timesteps, bool verbose = true, bool write_counts = false, bool write_nns = false, bool write_latt = false, int write_step = 20, int write_version_no = 0, std::string dir=".", bool periodic_bc=false);

		/********************
		Write/Read lattice
		********************/

		void write_lattice(int index, int write_version_no, std::string dir);
		void read_lattice(std::string fname);
	};

	/****************************************
	Main simulation IMPLEMENTATION DEFINITIONS
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	// Constructors
	Simulation::Impl::Impl(double dt, int box_length, int dim)
	{
		_dim = dim;
		// Make lattice
		if (_dim == 1) {
			_lattice = new Lattice(box_length);
		} else if (_dim == 2) {
			_lattice = new Lattice(box_length,box_length);
		} else if (_dim == 3) {
			_lattice = new Lattice(box_length,box_length,box_length);
		} else {
			std::cerr << "ERROR: only lattice dimensions 1,2,3 are supported" << std::endl;
			exit(EXIT_FAILURE);
		};

		_uni_next = nullptr;
		_t = 0.0;
		_t_step = 0;
		_dt = dt;
		_box_length = box_length;
	};

	Simulation::Impl::Impl(Impl&& other) {
		_move(other);
	};
    Simulation::Impl& Simulation::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };

	// Destructor
	Simulation::Impl::~Impl() {
		_clean_up();
	};

	// Helpers
	void Simulation::Impl::_clean_up() {
		if (_lattice != nullptr) {
			delete _lattice;
		};
        
        for (auto sp: _species) {
            if (sp) {
                delete sp;
                sp = nullptr;
            };
        };

        for (auto b: _bi_rxns) {
            if (b) {
                delete b;
                b = nullptr;
            };
        };

        for (auto u: _uni_rxns) {
            if (u) {
                delete u;
                u = nullptr;
            };
        };
	};
	void Simulation::Impl::_move(Impl& other) {
		_dim = other._dim;
		_box_length = other._box_length;
		_lattice = new Lattice(*(other._lattice));
		_species = other._species;
		_bi_rxns = other._bi_rxns;
		_uni_rxns = other._uni_rxns;
		_t = other._t;
		_t_step = other._t_step;
		_dt = other._dt;
		_t_uni_next = other._t_uni_next;
        _dir_write = other._dir_write;
        _uni_next = other._uni_next;
        
        // Clear the other
        other._dim = 3;
        other._lattice = nullptr;
        other._species.clear();
        other._bi_rxns.clear();
        other._uni_rxns.clear();
        other._box_length = 0;
        other._t = 0.0;
        other._t_step = 0;
        other._dt = 0.0;
        other._t_uni_next = 0.0;
        other._dir_write = "";
        other._uni_next = nullptr;
	};

	/********************
	Add a species
	********************/

	void Simulation::Impl::add_species(std::string name, bool conserved) {
		// Add
		_species.push_back(new Species(name,conserved));
        for (auto u: _uni_rxns) {
			_species.back()->add_rxn(u);
		};
        for (auto b: _bi_rxns) {
			_species.back()->add_rxn(b);
		};
	};

	/********************
	 Add a reaction
	********************/

	// Unimolecular rxn
	void Simulation::Impl::add_uni_rxn(std::string name, double kr, std::string r) {
		// Find the species
		Species *sr = _find_species(r);
		// Make the reaction
		_uni_rxns.push_back(new UniReaction(name,kr,sr));
		// Add reaction to all species
        for (auto species: _species) {
			species->add_rxn(_uni_rxns.back());
		};
	};
	void Simulation::Impl::add_uni_rxn(std::string name, double kr, std::string r, std::string p) {
        // Check against empty
        if (r=="" && p=="") {
            std::cerr << "Error add_uni_rxn: reactant and product cannot be nullptr" << std::endl;
            exit(EXIT_FAILURE);
        };
        
		// Find the species
		Species *sr = _find_species(r);
		Species *sp = _find_species(p);
		// Make the reaction
		_uni_rxns.push_back(new UniReaction(name,kr,sr,sp));
		// Add reaction to all species
        for (auto species: _species) {
			species->add_rxn(_uni_rxns.back());
		};
	};
	void Simulation::Impl::add_uni_rxn(std::string name, double kr, std::string r, std::string p1, std::string p2) {
        if (r=="") {
            std::cerr << "Error add_uni_rxn: Only 0->A is supported, not 0->A+B" << std::endl;
            exit(EXIT_FAILURE);
        };

		// Find the species
		Species *sr = _find_species(r);
		Species *sp1 = _find_species(p1);
		Species *sp2 = _find_species(p2);
		// Make the reaction
		_uni_rxns.push_back(new UniReaction(name,kr,sr,sp1,sp2));
		// Add reaction to all species
        for (auto species: _species) {
			species->add_rxn(_uni_rxns.back());
		};
	};

	// Bimolecular rxn
	void Simulation::Impl::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2) {
		// Find the species
		Species *sr1 = _find_species(r1);
		Species *sr2 = _find_species(r2);
		// Make the reaction
		_bi_rxns.push_back(new BiReaction(name,prob,sr1,sr2));
		// Add reaction to all species
        for (auto species: _species) {
			species->add_rxn(_bi_rxns.back());
		};
	};
	void Simulation::Impl::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p) {
		// Find the species
		Species *sr1 = _find_species(r1);
		Species *sr2 = _find_species(r2);
		Species *sp = _find_species(p);
		// Make the reaction
		_bi_rxns.push_back(new BiReaction(name,prob,sr1,sr2,sp));
		// Add reaction to all species
        for (auto species: _species) {
			species->add_rxn(_bi_rxns.back());
		};
	};
	void Simulation::Impl::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2) {
		// Find the species
		Species *sr1 = _find_species(r1);
		Species *sr2 = _find_species(r2);
		Species *sp1 = _find_species(p1);
		Species *sp2 = _find_species(p2);
		// Make the reaction
		_bi_rxns.push_back(new BiReaction(name,prob,sr1,sr2,sp1,sp2));
		// Add reaction to all species
        for (auto species: _species) {
			species->add_rxn(_bi_rxns.back());
		};
	};

	/********************
	Populate lattice
	********************/

    void Simulation::Impl::populate_random(std::vector<std::string> sps) {
        // Random initial counts
        int n_possible = pow(_box_length,_dim);
        std::map<std::string,int> counts0;
        for (auto sp : sps) {
            counts0[sp] = randI(0,n_possible);
            n_possible -= counts0[sp];
            if (n_possible < 0) { n_possible = 0; };
        };
        
        // Random initial lattice
        populate_lattice(counts0);
    };
    
	void Simulation::Impl::populate_lattice(std::map<std::string,int> counts) {
	    // Go through all species
		Species *s;
		for (auto c: counts) 
		{
			s = _find_species(c.first);
	    	for (auto i=0; i<c.second; i++) 
	    	{
	    		// Make
		    	_lattice->make_mol_random(s);	
	    	};
		};
	};

    std::vector<std::string> Simulation::Impl::_vec_from_dict_keys(const dict1 &dict) const {
        std::vector<std::string> sps;
        for (auto pr: dict) {
            sps.push_back(pr.first);
        };
        return sps;
    };
    sdict1* Simulation::Impl::_convert_dict_to_sdict(const dict1 &h_dict) const {
        sdict1* h_dict_sp = new sdict1();
        for (auto const &hpr: h_dict) {
            (*h_dict_sp)[_find_species(hpr.first)] = hpr.second;
        };
        return h_dict_sp;
    };
    sdict2* Simulation::Impl::_convert_dict_to_sdict(const dict2 &j_dict) const {
        sdict2* j_dict_sp = new sdict2();
        for (auto jpr1: j_dict) {
            for (auto jpr2: jpr1.second) {
                (*j_dict_sp)[_find_species(jpr1.first)][_find_species(jpr2.first)] = jpr2.second;
            };
        };
        return j_dict_sp;
    };
    sdict3* Simulation::Impl::_convert_dict_to_sdict(const dict3 &k_dict) const {
        sdict3* k_dict_sp = new sdict3();
        for (auto kpr1: k_dict) {
            for (auto kpr2: kpr1.second) {
                for (auto kpr3: kpr2.second) {
                    (*k_dict_sp)[_find_species(kpr1.first)][_find_species(kpr2.first)][_find_species(kpr3.first)] = kpr3.second;
                };
            };
        };
        return k_dict_sp;
    };
    sdict4* Simulation::Impl::_convert_dict_to_sdict(const dict4 &q_dict) const {
        sdict4* q_dict_sp = new sdict4();
        for (auto qpr1: q_dict) {
            for (auto qpr2: qpr1.second) {
                for (auto qpr3: qpr2.second) {
                    for (auto qpr4: qpr3.second) {
                        (*q_dict_sp)[_find_species(qpr1.first)][_find_species(qpr2.first)][_find_species(qpr3.first)][_find_species(qpr4.first)] = qpr4.second;
                    };
                };
            };
        };
        return q_dict_sp;
    };

    void Simulation::Impl::populate_lattice_1d(const dict1 &h_dict, int n_steps) {
        // Start by populating lattice randomly
        populate_random(_vec_from_dict_keys(h_dict));
        
        // Run
        auto h_dict_sp = _convert_dict_to_sdict(h_dict);
        _lattice->sample_1d(h_dict_sp,n_steps);
        delete h_dict_sp;
    };
    
	void Simulation::Impl::populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, int n_steps) {
        // Start by populating lattice randomly
        populate_random(_vec_from_dict_keys(h_dict));
        
        // Run
        auto h_dict_sp = _convert_dict_to_sdict(h_dict);
        auto j_dict_sp = _convert_dict_to_sdict(j_dict);
        _lattice->sample_1d(h_dict_sp,j_dict_sp,n_steps);
        delete h_dict_sp;
        delete j_dict_sp;
    };

	void Simulation::Impl::populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, const dict3 &k_dict, int n_steps) {
        // Start by populating lattice randomly
        populate_random(_vec_from_dict_keys(h_dict));
        
        // Run
        auto h_dict_sp = _convert_dict_to_sdict(h_dict);
        auto j_dict_sp = _convert_dict_to_sdict(j_dict);
        auto k_dict_sp = _convert_dict_to_sdict(k_dict);
        _lattice->sample_1d(h_dict_sp,j_dict_sp,k_dict_sp,n_steps);
        delete h_dict_sp;
        delete j_dict_sp;
        delete k_dict_sp;
	};


	void Simulation::Impl::populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, const dict3 &k_dict, const dict4 &q_dict, int n_steps)
	{
        // Start by populating lattice randomly
        populate_random(_vec_from_dict_keys(h_dict));
        
        // Run
        auto h_dict_sp = _convert_dict_to_sdict(h_dict);
        auto j_dict_sp = _convert_dict_to_sdict(j_dict);
        auto k_dict_sp = _convert_dict_to_sdict(k_dict);
        auto q_dict_sp = _convert_dict_to_sdict(q_dict);
        _lattice->sample_1d(h_dict_sp,j_dict_sp,k_dict_sp,q_dict_sp,n_steps);
        delete h_dict_sp;
        delete j_dict_sp;
        delete k_dict_sp;
        delete q_dict_sp;
	};

	/********************
	Do a uni reaction
	********************/

	void Simulation::Impl::do_uni_rxn(UniReaction *rxn) {

		// Declarations
		std::pair<bool,SiteIt3D> get_it;
        std::pair<bool,Site3D> get_empty;
		SiteIt3D sit;
		Site3D s,snbr;
		std::pair<bool,Site3D> free_pair;

		// Try to do the reaction at several sites
		// Failure can arise if there is not enough room for the products
		int ctr = 0;
		while (ctr < 20) {

            if (rxn->r) {
                
                // Grab a random site with this species
                get_it = _lattice->get_mol_random_it(rxn->r);
                if (!(get_it.first)) {
                    return; // Impossible; stop
                };
                sit = get_it.second;
                s = Site3D(sit);
                
                // Check if there is room for the products
                if (rxn->p.size() == 2) {
                    free_pair = _lattice->get_free_neighbor_random(sit);
                    if (!(free_pair.first)) {
                        // Not enough room for the two products; try again with a different random site
                        ctr += 1;
                        continue;
                    } else {
                        snbr = free_pair.second;
                    };
                };
                
                // Remove the reactant
                _lattice->erase_mol_it(sit);

                // Conserve reactants
                if (rxn->r->conserved) {
                    _lattice->make_mol_random(rxn->r);
                };
                
                // Place products, if needed at the neighbor site
                if (rxn->p.size() == 1) {
                    _lattice->make_mol(s, rxn->p[0]);
                } else if (rxn->p.size() == 2) {
                    _lattice->make_mol(s, rxn->p[0]);
                    _lattice->make_mol(snbr, rxn->p[1]);
                };

            } else {
                
                // Reaction is type 0->A
                // Get an empty site!
                get_empty = _lattice->get_free_site();
                if (!(get_empty.first)) {
                    return; // Impossible; stop
                };
                s = get_empty.second;
                
                // Place products; only 1 is allowed anyways
                _lattice->make_mol(s, rxn->p[0]);
            };

			// Conserve products
			for (auto p: rxn->p)
			{
				if (p->conserved) {
					_lattice->erase_mol_random(p);
				};
			};

			// Sucess
			return;
		};
	};

	/********************
	Diffuse all the mols and do bimol reactions
	********************/

	void Simulation::Impl::diffuse_mols(bool periodic)
	{
		// Copy the old map
		Lattice *todo = new Lattice(*_lattice);

		// Clear the current lattice
		_lattice->clear();

		// Declarations needed
		std::pair<bool,SiteIt3D> old_pair;
		SiteIt3D sOldIt,sNewIt;
		Site3D sNew,sOld;
		Mol *mOld,*mColl;
		std::pair<Site3D,std::pair<bool,SiteIt3D>> nbr_pair;
		bool occ_todo, occ_done;
		std::pair<bool,SiteIt3D> coll_done_pair;
		std::pair<bool,BiReaction*> rxn_pair;
		BiReaction *rxn; 

		// Go over all mols to move
		while (todo->size() > 0) {

			if (DIAG_DIFFUSE) { std::cout << std::endl; };

			// Reset
			occ_todo = false;
			occ_done = false;

			// Grab some element
			old_pair = todo->get_mol_random_it();

			sOldIt = old_pair.second;
			mOld = &(sOldIt.it_1D->second);

			if (DIAG_DIFFUSE) { std::cout << "diffuse_mols: got element..." << std::flush; };

			// Move
			nbr_pair = todo->get_neighbor_random(sOldIt, periodic);
			sNew = nbr_pair.first;

			if (DIAG_DIFFUSE) { std::cout << "got neighbor..." << std::flush; };

			// Check if occupied in todo lattice
			occ_todo = nbr_pair.second.first;
			if (occ_todo) {
				// Yes; its collided with something in the todo pile
				sNewIt = nbr_pair.second.second;
				mColl = &(sNewIt.it_1D->second);
			} else {
				// No; check if it's collided with something in the done pile
				coll_done_pair = _lattice->get_mol_it(sNew);
				occ_done = coll_done_pair.first;
				if (occ_done) {
					// Yes; its collided with something in the done pile
					sNewIt = coll_done_pair.second;
					mColl = &(sNewIt.it_1D->second);
				};
			};

			if (DIAG_DIFFUSE) { std::cout << "check colls..." << std::flush; };

			// If unoccupied, just commit the move (diffuse)
			if (!occ_todo && !occ_done) {
				// Move
				_lattice->make_mol(sNew, mOld->sp);
				// Continue
				todo->erase_mol_it(sOldIt);
				continue;
			};

			// Collision; check if reaction occurs
			rxn_pair = mOld->check_bi_rxns_mol(mColl);

			if (DIAG_DIFFUSE) { std::cout << "checked rxn..." << std::flush; };

			// No reaction occurred?
			if (!(rxn_pair.first)) {
				// No reaction & don't move
				_lattice->make_mol(Site3D(sOldIt), mOld->sp);
				// Continue
				todo->erase_mol_it(sOldIt);
				continue;
			};

			// Reaction occurred
			rxn = rxn_pair.second;
			rxn->count++;

			if (DIAG_DIFFUSE) { std::cout << "rxn occurred..." << std::flush; };

			// Remove the reactants
			sOld = Site3D(sOldIt); // grab it before it's erased
			todo->erase_mol_it(sOldIt);
			if (occ_todo) { 
				// Note: the old iterator sNewIt has been invalidated :(
				todo->erase_mol(sNew);
			} else if (occ_done) { 
				_lattice->erase_mol(sNewIt); 
			};

			if (DIAG_DIFFUSE) { std::cout << "removed r..." << std::flush; };

			// Conserve reactants
			if (rxn->r1->conserved) {
				_lattice->make_mol_random(rxn->r1);
			};
			if (rxn->r2->conserved) {
				_lattice->make_mol_random(rxn->r2);
			};

			if (DIAG_DIFFUSE) { std::cout << "conserved r..." << std::flush; };

			// Place products, if needed at the old site
			if (rxn->p.size() == 1) {
				_lattice->make_mol(sNew, rxn->p[0]);
			} else if (rxn->p.size() == 2) {
				_lattice->make_mol(sNew, rxn->p[0]);
				_lattice->make_mol(sOld, rxn->p[1]);
			};

			if (DIAG_DIFFUSE) { std::cout << "added p..." << std::flush; };

			// Conserve products
			for (auto p: rxn->p)
			{
				if (p->conserved) {
					_lattice->erase_mol_random(p);
				};
			};

			if (DIAG_DIFFUSE) { std::cout << "conserved p..." << std::flush; };

			// Finished
		};

		// Clear old
		delete todo;
	};

    /********************
     Set directory
     ********************/
    
    void Simulation::Impl::set_write_dir(std::string dir) {
        _dir_write = dir;
    };
    
	/********************
	Run simulation
	********************/

	void Simulation::Impl::run(int n_timesteps, bool verbose, bool write_counts, bool write_nns, bool write_latt, int write_step, int write_version_no, std::string dir, bool periodic_bc)
	{
		// Clear data in files if writing
		std::ofstream ofs;
		std::stringstream fname;
		if (write_counts || write_nns || write_latt) {
			// Clear counts
			for (auto s: _species) {
				fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/counts/" << s->name << ".txt";
				ofs.open(fname.str(), std::ofstream::out | std::ofstream::trunc);
				ofs.close();
				fname.str("");
			};
			// Clear nns
            for (auto i1=0; i1<_species.size(); i1++) {
                for (auto i2=i1; i2<_species.size(); i2++) {
                    fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/nns/" << _species[i1]->name << "_" << _species[i2]->name << ".txt";
                    ofs.open(fname.str(), std::ofstream::out | std::ofstream::trunc);
                    ofs.close();
                    fname.str("");
                };
            };
            
			// Clear lattice data
			fname << "exec rm -r ./" << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/lattice/*";
			system(fname.str().c_str());
			fname.str("");
		};

		// Go over all timesteps
		double t_next;
		for (int i_step=0; i_step < n_timesteps; i_step++) 
		{
			// The next time
			t_next = _t + _dt;

			// Print if needed
			if (_t_step % 1 == 0) {
				if (verbose) {
					std::cout << "Time: " << _t << " / " << n_timesteps*_dt << " " << _species << std::endl;
				};
			};

			// Write if needed
			if (_t_step % write_step == 0) {
				if (write_counts) {
					// Write counts to file
					for (auto s: _species) {
						fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/counts/" << s->name << ".txt";
						ofs.open(fname.str(), std::ofstream::app);
						ofs << _t << " " << s->count << "\n";
						ofs.close();
						fname.str("");
					};
				};
				if (write_nns) {
					// Write NNs to file
                    for (auto i1=0; i1<_species.size(); i1++) {
                        for (auto i2=i1; i2<_species.size(); i2++) {
							fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/nns/" << _species[i1]->name << "_" << _species[i2]->name << ".txt";
							ofs.open(fname.str(), std::ofstream::app);
							ofs << _t << " " << _lattice->get_nn(_species[i1],_species[i2]) << "\n";
							ofs.close();
							fname.str("");
						};
					};
				};
				if (write_latt) {
					// Write the lattice
					write_lattice(int(_t_step/write_step), write_version_no, dir);
				};
			};

			// Do we need to schedule unimolecular reactions?
			// (Either start of sim, or there were none available before)
			if (_uni_next == nullptr) {
				// Schedule unimolecular reactions
				_schedule_uni();
			};

			// Do uni reactions
			while (_uni_next != nullptr && _t_uni_next < t_next) {
				// Do it
				do_uni_rxn(_uni_next);
				// Advance time
				_t = _t_uni_next;
				// Schedule
				_schedule_uni();
			};

			// Diffuse and do bimol reactions
			diffuse_mols(periodic_bc);

			// It is now the next time
			_t = t_next;
			_t_step++;

		};
	};

	/********************
	Write/Read lattice
	********************/

	void Simulation::Impl::write_lattice(int index, int write_version_no, std::string dir)
	{
		std::stringstream fname;
		fname << dir << "/lattice_v" << std::setfill('0') << std::setw(3) << write_version_no << "/lattice/" << std::setfill('0') << std::setw(4) << index << ".txt";
		_lattice->write_to_file(fname.str());
		fname.str("");
	};

	void Simulation::Impl::read_lattice(std::string fname) {
		std::map<std::string,Species*> sp_map;
        for (auto s: _species) {
			sp_map[s->name] = s;
		};
		_lattice->read_from_file(fname, sp_map);
	};

	/****************************************
	Main simulation - PRIVATE
	****************************************/

	/********************
	Find a species by name
	********************/

	Species* Simulation::Impl::_find_species(std::string name) const {
		// Go through the species
        for (auto s: _species) {
			if (s->name == name) {
				return s;
			};
		};
		return nullptr;
	};

	/********************
	Schedule the next uni reaction
	********************/

	void Simulation::Impl::_schedule_uni() {
		double props_cum = 0.0;
		std::vector<double> props;
		props.push_back(0.0);

		// Go through all possible reagants, calculate propensities
		for (auto u: _uni_rxns)
		{
            if (u->r) {
                props_cum += u->r->count * u->kr;
            } else {
                // 0 -> A ; use A
                props_cum += u->p[0]->count * u->kr;
            };
			props.push_back(props_cum);
		};

		// Check that at least one reaction is possible
		if (!(props_cum > 0)) {
			_t_uni_next = _t;
			_uni_next = nullptr; // Based on nullptr, will check again later
			return;
		};

		// Choose a reaction
		double r = randD(0.0,props_cum);
        int i = 0;
		int n = 0;
		while (!(props[n] < r && r < props[n+1])) {
			n++;
			i += 1;
		};
		_uni_next = _uni_rxns[i];

		// Time of next reaction
		_t_uni_next = _t + log(1.0/randD(0.0,1.0))/props_cum;
	};

	/****************************************
	Main simulation IMPL forwards
	****************************************/

	/********************
	Constructor/Destructor
	********************/

	Simulation::Simulation(double dt, int box_length, int dim) : _impl(new Impl(dt,box_length,dim)) {};
	Simulation::Simulation(Simulation&& other) = default; // movable but no copies
    Simulation& Simulation::operator=(Simulation&& other) = default; // movable but no copies
	Simulation::~Simulation() = default;

	/********************
	Add species
	********************/

	void Simulation::add_species(std::string name, bool conserved) {
		_impl->add_species(name,conserved);
	};

	/********************
	 Add a reaction
	********************/

	// Unimolecular rxn
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r) {
		_impl->add_uni_rxn(name,kr,r);
	};
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r, std::string p) {
		_impl->add_uni_rxn(name,kr,r,p);
	};
	void Simulation::add_uni_rxn(std::string name, double kr, std::string r, std::string p1, std::string p2) {
		_impl->add_uni_rxn(name,kr,r,p1,p2);
	};

	// Bimolecular rxn
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2) {
		_impl->add_bi_rxn(name,prob,r1,r2);
	};
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p) {
		_impl->add_bi_rxn(name,prob,r1,r2,p);
	};
	void Simulation::add_bi_rxn(std::string name, double prob, std::string r1, std::string r2, std::string p1, std::string p2) {
		_impl->add_bi_rxn(name,prob,r1,r2,p1,p2);
	};

	/********************
	Populate lattice
	********************/

    void Simulation::populate_random(std::vector<std::string> sps) {
        _impl->populate_random(sps);
    };
	void Simulation::populate_lattice(std::map<std::string,int> counts) {
		_impl->populate_lattice(counts);
	};
	void Simulation::populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, int n_steps) {
		_impl->populate_lattice_1d(h_dict,j_dict,n_steps);
	};
	void Simulation::populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, const dict3 &k_dict, int n_steps) {
		_impl->populate_lattice_1d(h_dict,j_dict,k_dict,n_steps);
	};
	void Simulation::populate_lattice_1d(const dict1 &h_dict, const dict2 &j_dict, const dict3 &k_dict, const dict4 &q_dict, int n_steps) {

		_impl->populate_lattice_1d(h_dict,j_dict,k_dict,q_dict,n_steps);
	};
    
    /********************
     Set directory
     ********************/
    
    void Simulation::set_write_dir(std::string dir) {
        _impl->set_write_dir(dir);
    };

	/********************
	Run simulation
	********************/

	void Simulation::run(int n_timesteps, bool verbose, bool write_counts, bool write_nns, bool write_latt, int write_step, int write_version_no, std::string dir, bool periodic_bc) {
		_impl->run(n_timesteps,verbose,write_counts,write_nns,write_latt,write_step,write_version_no, dir, periodic_bc);
	};

	/********************
	Write/Read lattice
	********************/

	void Simulation::write_lattice(int index, int write_version_no, std::string dir) {
		_impl->write_lattice(index,write_version_no, dir);
	};
	void Simulation::read_lattice(std::string fname) {
		_impl->read_lattice(fname);
	};


};
