#ifndef SPECIES_h
#define SPECIES_h
#include "species.hpp"
#endif

// Other Gillespie3D

#ifndef LATTICE_h
#define LATTICE_h
#include "lattice.hpp"
#endif

/************************************
* Namespace for Gillespie3D
************************************/

namespace Gillespie3D {

	/****************************************
	Species
	****************************************/
	
	Species::Species(std::string nameIn, bool conservedIn) {
		this->name = nameIn;
		this->conserved = conservedIn;
		this->count = 0;
	};

	// Comparator
	bool operator <(const Species& a, const Species& b) {
		return a.name < b.name;
	};

	/********************
	Add a reaction, if appropriate
	********************/

	void Species::add_rxn(BiReaction* rxn) {
		if (rxn->r1 == this) 
		{ 
			this->bi_rxns[rxn->r2].push_back(rxn); 
		} else if (rxn->r2 == this) 
		{ 
			this->bi_rxns[rxn->r1].push_back(rxn); 
		};
	};
	void Species::add_rxn(UniReaction* rxn) {
		if (rxn->r == this) { this->uni_rxns.push_back(rxn); };
	};

	/********************
	Check if any reactions are possible; if so, return a random one
	********************/

	std::pair<bool,BiReaction*> Species::check_bi_rxns_mol(Mol &other) {
		std::map<Species*,std::vector<BiReaction*>>::iterator itm;
		std::vector<BiReaction*>::iterator itb;

		// Check all reactions on this species
		itm = this->bi_rxns.find(other.sp);
		if (itm != this->bi_rxns.end()) {
			// There is a reaction - check if any succeed
			for (auto ib=0; ib < itm->second.size(); ib++)
			{
				// Check probability
				if (randD(0.0,1.0) < itm->second[ib]->prob)
				{
					return std::make_pair(true,itm->second[ib]);
				};
			};
		};
		return std::make_pair(false,nullptr);
	};

	std::pair<bool,std::pair<Mol*,BiReaction*>> Species::check_bi_rxns_mols(mvec &other) {
		std::pair<bool,BiReaction*> ret;

		// Go through all other mols
		for (auto im = 0; im < other.size(); im++) {
			ret = check_bi_rxns_mol(other[im]);
			if (ret.first) {
				return std::make_pair(true,std::make_pair(&(other[im]),ret.second));
			};
		};
		return std::make_pair(false,std::make_pair(nullptr,nullptr));
	};

	/****************************************
	Mol
	****************************************/

	Mol::Mol(Species *spIn) {
		this->sp = spIn;
		this->has_moved = false;
	};

	// Comparator
	bool operator <(const Mol& a, const Mol& b) {
    	return a.sp->name < b.sp->name;
	};

	/********************
	Check if any reactions are possible; if so, return a random one
	********************/

	std::pair<bool,BiReaction*> Mol::check_bi_rxns_mol(Mol &other) {
		return this->sp->check_bi_rxns_mol(other);
	};
	std::pair<bool,std::pair<Mol*,BiReaction*>> Mol::check_bi_rxns_mols(mvec &other) {
		return this->sp->check_bi_rxns_mols(other);
	};
};