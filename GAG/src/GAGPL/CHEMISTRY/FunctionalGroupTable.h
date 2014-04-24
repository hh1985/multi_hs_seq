/*
 * =====================================================================================
 *
 *       Filename:  FunctionalGroupTable.h
 *
 *    Description:  Functional group table loaded from xml drive.
 *
 *        Version:  1.0
 *        Created:  04/29/2012  9:51:00 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */


#ifndef	GAG_FUNCTIONALGROUPTABLE_H
#define GAG_FUNCTIONALGROUPTABLE_H

#include <GAGPL/MISC/ConfigLoader.h>
#include <GAGPL/CHEMISTRY/FunctionalGroup.h>
//#include <GAGPL/CHEMISTRY/Modification.h>
#include <boost/noncopyable.hpp>

namespace gag
{
	struct FunctionalGroupProtoType
	{
		std::string _symbol;
		std::string _name;
		std::string _composition_string;
		size_t _start;
		size_t _end;
		// Simply store the name of the functional groups.
		// The pair represents core and sub-functional group.
		std::vector<std::pair<std::string, std::vector<std::string> > > _sites;
	};

	typedef std::map<std::string, FunctionalGroupProtoType> FunctionalGroupReference;

	class FunctionalGroupTable: public ConfigLoader, private boost::noncopyable
	{
		private:
			//boost::property_tree::ptree params;
			std::map<std::string, FunctionalGroup> functionalgroups;

			// Map FunctionalGroupReference object to functionalgroups.
			void buildtree(const FunctionalGroupReference& fg_ref);
			// If this functional group has been defined in functionalgroups, return it directly. Otherwise, create it from the scratch.
			FunctionalGroup createFunctionalGroup(const FunctionalGroupReference& fg_ref, const std::string& symbol); 

		protected:
			FunctionalGroupTable() {
				load();
			}

		public:
			
			static FunctionalGroupTable& Instance();
			bool containFunctionalGroup(const std::string& symbol);

			void load(const std::string& filename = "./config/commongroup.xml");

			//// The modification will locate to the right position based on its own definition.
			//void appendModificationToFunctionalGroup(FunctionalGroup& fg, Modification& mod);

			//// Specify the position of the modification. This function might not apply to some of the modifications (e.g. 2-AB labeling). The copy of mod will be modified by idx.
			//void appendModificationToFunctionalGroup(FunctionalGroup& fg, Modification mod, const size_t idx);

			// Conversion from string to FunctionalGroupProtoType, and to FunctionalGroup.
			FunctionalGroup getFunctionalGroupBySymbol(const std::string&) const;
			FunctionalGroup getFunctionalGroupByName(const std::string&);		
			std::vector<std::string> getAllFunctionalGroupSymbols() const;
	};
}



#endif   /* ----- #ifndef GAG_FUNCTIONALGROUPTABLE_H ----- */
