/*
 * =====================================================================================
 *
 *       Filename:  InternalSite.h
 *
 *    Description:  The internal site inside a monosaccharide unit.
 *
 *        Version:  1.0
 *        Created:  04/30/2012  3:50:55 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef  GAG_INTERNALSITE_H
#define  GAG_INTERNALSITE_H

#include <GAGPL/CHEMISTRY/FunctionalGroup.h>
#include <set>

namespace gag
{

	class InternalSite: public Unit 
	{
		private:
			// This information should come from the xml file.
			// But the CompositionGroup should be a copy instead of singleton.
			// The only reason using FunctionalGroup instead of Composition is that
			// The name can be used for verification purposes.
			std::multiset<FunctionalGroup> _groups;
		
		public:
	
		InternalSite() {}
		InternalSite(Composition& compo)
			: Unit(compo)
		{}
		InternalSite(std::string str)
			: Unit(str)
		{}
		// No modification of the composition.
		inline void addFunctionalGroup(const FunctionalGroup& fg)
		{
			// TBD: Exception of repeating keys.
			_groups.insert(fg);
		}
		inline void removeFunctionalGroup(const FunctionalGroup& fg)
		{
			_groups.erase(fg);
		}

		void addFunctionalGroup(const std::string& str);

		// Specify the original fg, the fg that will add to it, and the lost composition of all.
		// High level.
		void modify(FunctionalGroup& ori, FunctionalGroup& sub, Composition& loss);
		// Specify the original fg, and the overall composition that will add to it.
		// Low level.
		//void add(FunctionalGroup& ori, const Composition& compo);

		// The problem is: how to specify which functionalgroup it is going to remove.
		// If there is repeating fgs.
		//void remove(FunctionalGroup& ori, const Composition& compo);
		//void remove(FunctionalGroup& ori);
		void replace(FunctionalGroup& ori, FunctionalGroup& replace);

		// Get the tracked functional group and take it as the core for addition.
		void addFunctionalGroup(std::vector<FunctionalGroup>& fg_chain, FunctionalGroup& fg);
		void removeFunctionalGroup(std::vector<FunctionalGroup>& fg_chain);

		// This function is used to specify the atom shift for B/Y cleavage.
		Composition getCleavageShift(FunctionalGroup& fg, const size_t idx);
		// Simply get the first non-H functional group.
		Composition getCleavageShift(const size_t idx);

		inline std::multiset<FunctionalGroup>& getFunctionalGroups()
		{
			return _groups;
		}

		bool hasFunctionalGroup(FunctionalGroup& fg);
		bool hasFunctionalGroup(const std::string& str);
		
	};
}

#endif   /* ----- #ifndef GAG_INTERNALSITE_H ---- */
