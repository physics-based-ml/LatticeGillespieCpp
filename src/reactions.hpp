#include <string>
#include <vector>

/************************************
* Namespace for lattg
************************************/

namespace lattg {

	/****************************************
	Necessary declarations
	****************************************/

	struct Species;

	/****************************************
	Unireaction
	****************************************/

	struct UniReaction
	{
		std::string name;
		std::vector<Species*> p;
		Species* r;
		double kr;
		int count;

		UniReaction(std::string nameIn, double krIn, Species* rIn);
		UniReaction(std::string nameIn, double krIn, Species* rIn, Species* pIn);
		UniReaction(std::string nameIn, double krIn, Species* rIn, Species* p1In, Species* p2In);
	};

	/****************************************
	Bireaction
	****************************************/

	struct BiReaction
	{
		std::string name;
		std::vector<Species*> p;
		Species *r1;
		Species *r2;
		double prob;
		int count;

		BiReaction(std::string nameIn, double probIn, Species* r1In, Species* r2In);
		BiReaction(std::string nameIn, double probIn, Species* r1In, Species* r2In, Species* pIn);
		BiReaction(std::string nameIn, double probIn, Species* r1In, Species* r2In, Species* p1In, Species* p2In);
	};

};