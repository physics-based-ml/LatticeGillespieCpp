#include "species.hpp"

/************************************
* Namespace for Gillespie1D
************************************/

namespace Gillespie1D {

	/****************************************
	General functions
	****************************************/

	/********************
	Random numbers
	********************/

	double randD(double dMin, double dMax)
	{
	    return dMin + ((double)rand() / RAND_MAX) * (dMax - dMin);
	};
	int randI(int iMin, int iMax)
	{
		return iMin + rand() % (iMax - iMin + 1);
	};

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

	std::pair<bool,BiReaction*> Species::check_bi_rxns_mol(Mol *other) {
		std::map<Species*,std::vector<BiReaction*>>::iterator itm;
		std::vector<BiReaction*>::iterator itb;

		// Check all reactions on this species
		itm = this->bi_rxns.find(other->sp);
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

	/****************************************
	Mol
	****************************************/

	Mol::Mol(Species *spIn) {
		this->sp = spIn;
	};

	// Comparator
	bool operator <(const Mol& a, const Mol& b) {
    	return a.sp->name < b.sp->name;
	};

	/********************
	Check if any reactions are possible; if so, return a random one
	********************/

	std::pair<bool,BiReaction*> Mol::check_bi_rxns_mol(Mol *other) {
		return this->sp->check_bi_rxns_mol(other);
	};

};