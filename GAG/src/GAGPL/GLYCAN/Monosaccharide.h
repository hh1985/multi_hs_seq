/*
 * =====================================================================================
 *
 *       Filename:  Monosaccharide.h
 *
 *    Description:  The definition of monosaccaride unit.
 *
 *        Version:  1.0
 *        Created:  04/30/2012  9:50:40 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */


#ifndef  GAG_MONOSACCHARIDE_H
#define  GAG_MONOSACCHARIDE_H

//#include <GAGPL/CHEMISTRY/Unit.h>
//#include <GAGPL/GLYCAN/InternalSite.h>
//#include <set>
#include <GAGPL/CHEMISTRY/FunctionalGroup.h>

namespace gag
{
	typedef FunctionalGroup Monosaccharide;
	//class Monosaccharide: public Unit
	//{
	//	private:
	//		std::string symbol;
	//		std::string name;

	//		size_t ring_start;
	//		size_t ring_end;
	//		// I double the use of InternalSite here.
	//		std::vector<InternalSite> internalsites;

	//	public:

	//		Monosaccharide() {}

	//		// str -- composition.
	//		Monosaccharide(const std::string& str, const std::string sb, const std::string nm, const size_t start, const size_t end)
	//			: Unit(str), symbol(sb), name(nm), ring_start(start), ring_end(end),internalsites()
	//		{}

	//		inline const size_t& getRingStart() const
	//		{
	//			return ring_start;
	//		}
	//		inline const size_t& getRingEnd() const
	//		{
	//			return ring_end;
	//		}
	//		inline const std::string& getSymbol() const 
	//		{
	//			return symbol;
	//		}

	//		inline const std::string& getName() const
	//		{
	//			return name;
	//		}

	//		// CarbonID
	//		inline std::pair<size_t, size_t> getRingPosByCarbonID()
	//		{
	//			return std::pair<size_t, size_t>(ring_start, ring_end);
	//		}
	//		
	//		inline std::pair<size_t, size_t> getRingPosByRingID()
	//		{
	//			return std::pair<size_t, size_t>(1, ring_end-ring_start+1);
	//		}
	//		
	//		inline std::vector<InternalSite>& getInternalSites()
	//		{
	//			return internalsites;
	//		}

	//		// The Carbon ID starts from 1.
	//		inline InternalSite& getInternalSiteByCarbonID(const size_t c_id)
	//		{
	//			return internalsites.at(c_id);
	//		}
	//		// The Ring ID starts from 0.
	//		InternalSite& getInternalSiteByRingID(const size_t r_id);
	//		
	//		Composition getCleavageShift();

	//		size_t getRingID(const size_t c_id);
	//		size_t getCarbonID(const size_t r_id);

	//		inline void addInternalSite(const InternalSite& is)
	//		{
	//			internalsites.push_back(is);
	//		}

	//		bool hasFunctionalGroup(const size_t& pos, FunctionalGroup& fg);
	//		// No reference is allowed for temporary variable.
	//		Composition getSubCompositionByCarbonID(const size_t& s1, const size_t& s2);
	//		Composition getSubCompositionByRingID(const size_t& s1, const size_t& s2);

	//		/* It might be better to also add the functions into glycan sequence object.*/
	//		
	//		/* Just remember that FunctionalGroup in essence is a named composition with one or more functional groups */

	//		// pos is Carbon ID, and ori is the original functional group on the ring while 
	//		// plus is the added functional group and minus is the reduced functional group.
	//		void add(size_t pos, const FunctionalGroup& ori, const Composition& plus, const Composition& minus);
	//		void add(size_t pos, const std::string& ori, const std::string& plus, const std::string& minus);
	//		// No ori information needed.
	//		void addComposition(size_t pos, const std::string& plus);

	//		//void addFunctionalGroupByChain(FunctionalGroupChain& chain);

	//		// Specify the position of the site and the type of lost functionalgroup.
	//		// pos is Carbon ID.
	//		void remove(size_t pos, FunctionalGroup& fg, const Composition& lost);
	//		void remove(size_t pos, const std::string& ori, const std::string& lost);
	//		void remove(std::map<size_t, FunctionalGroup>& substract);
	//		// The functional group will be removed.
	//		// pos is Carbon ID.
	//		void remove(size_t pos, FunctionalGroup& fg);
	//		void removeFunctionalGroupByChain(FunctionalGroupChain& chain);
	//		
	//		void replace(size_t pos, FunctionalGroup& ori, FunctionalGroup& re);


	//};
}


#endif   /* ----- #ifndef GAG_MONOSACCHARIDE_H----- */
