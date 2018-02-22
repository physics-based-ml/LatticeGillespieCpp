#include "reactions.hpp"

/************************************
* Namespace for Gillespie3D
************************************/

namespace Gillespie3D {

	/****************************************
	Unireaction
	****************************************/

	UniReaction::UniReaction(std::string nameIn, double krIn, Species* rIn)
	{
		this->name = nameIn;
		this->kr = krIn;
		this->r = rIn;
		this->count = 0;
	};
	UniReaction::UniReaction(std::string nameIn, double krIn, Species* rIn, Species *pIn)
	{
		this->name = nameIn;
		this->kr = krIn;
		this->r = rIn;
		this->p.push_back(pIn);
		this->count = 0;
	};
	UniReaction::UniReaction(std::string nameIn, double krIn, Species* rIn, Species *p1In, Species *p2In)
	{
		this->name = nameIn;
		this->kr = krIn;
		this->r = rIn;
		this->p.push_back(p1In);
		this->p.push_back(p2In);
		this->count = 0;
	};

	/****************************************
	BiReaction
	****************************************/

	BiReaction::BiReaction(std::string nameIn, double probIn, Species* r1In, Species* r2In)
	{
		this->name = nameIn;
		this->prob = probIn;
		this->r1 = r1In;
		this->r2 = r2In;
		this->count = 0;
	};
	BiReaction::BiReaction(std::string nameIn, double probIn, Species* r1In, Species* r2In, Species* pIn)
	{
		this->name = nameIn;
		this->prob = probIn;
		this->r1 = r1In;
		this->r2 = r2In;
		this->p.push_back(pIn);
		this->count = 0;
	};
	BiReaction::BiReaction(std::string nameIn, double probIn, Species* r1In, Species* r2In, Species* p1In, Species* p2In)
	{
		this->name = nameIn;
		this->prob = probIn;
		this->r1 = r1In;
		this->r2 = r2In;
		this->p.push_back(p1In);
		this->p.push_back(p2In);
		this->count = 0;
	};

};