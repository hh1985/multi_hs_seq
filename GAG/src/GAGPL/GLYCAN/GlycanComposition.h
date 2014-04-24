/********************************************************************
	created:	2013/04/04
	created:	4:4:2013   9:18
	filename: 	GlycanComposition.h
	file path:	GAGPL/GLYCAN
	file base:	GlycanComposition
	file ext:	h
	author:		Han Hu (HH), hh.earlydays@gmail.com
	
	purpose:	general class for storing glycan composition.
*********************************************************************/

#ifndef GAG_GLYCANCOMPOSITION_INC
#define GAG_GLYCANCOMPOSITION_INC

//#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <string>
#include <map>

namespace gag
{
	class GlycanComposition
	{
	public:
		GlycanComposition() 
			: fgt(FunctionalGroupTable::Instance()){}
			
	public:
		// This function allows user to specify the number of composition for given monosaccharide residue. If the residue specified has been given, just update the number.
		void addGlycanComposition(const std::string& mono, int number);
		void addModificationComposition(const std::string& mod_symbol, int number);
		inline std::map<std::string, int>& getGlycanComposition()
		{
			return mono_compo;
		}
		//GlycanSequence getHSBackboneSequence(const std::string init_mono = "GlcA");
		// Get the modification map of the composition. Notice that in the future the functional group needs to extract the information of 
		inline std::map<std::string, int>& getModifications()
		{
			return mod_compo;
		}

	public:
		FunctionalGroupTable& fgt;
		std::map<std::string, int> mono_compo;
		std::map<std::string, int> mod_compo;

	};
}

#endif /* GAG_GLYCANCOMPOSITION_INC */