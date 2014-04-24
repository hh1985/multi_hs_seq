/*
 * =====================================================================================
 *
 *       Filename:  FunctionGroup.h
 *
 *    Description:  FunctionGroup
 *
 *        Version:  1.0
 *        Created:  4/23/2012 3:11:45 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef  GAG_FUNCTIONALGROUP_H
#define  GAG_FUNCTIONALGROUP_H

#include <GAGPL/CHEMISTRY/Composition.h>
#include <GAGPL/CHEMISTRY/Unit.h>

namespace gag
{
	// Forward declaration.
	class FunctionalGroup;

	// No need to maintain an overall copy of composition on the site structure.
	struct Site
	{
		Element core;
		std::multimap<std::string, FunctionalGroup> fg_map;

		Composition getSiteComposition();
		std::string getSiteCompositionString();

		inline bool attachFunctionalGroup(const std::string& fg_symbol) const
		{
			return fg_map.find(fg_symbol) != fg_map.end() ? true : false;
		}

		void printStructure();

	};

	// Recursively recording the position and corresponding functional group.
	typedef std::vector<std::pair<size_t, FunctionalGroup> > FunctionalGroupChain;
	// The definition of functional group should be something similar to InternalSites.
	class FunctionalGroup: public Unit
	{
		// This should be ordered from NRE to RE.
		//typedef std::vector<Composition> Sites;

		// The structure of FunctionalGroup is simple. All the reference and lookup function is moved to FunctionalGroupTable. It is reasonable, since in most cases, the only thing we need is the composition of the functional group.
	public:

		std::string _symbol;
		std::string _name;
		// If the functional group has underwent any modification. The basic idea is that the original functional group should be kept for the sake of reference.
		std::string _status;
		std::vector<Site> _sites;

		// By default, the ring_start and ring_end are not important. If there are any cleavages related to cross-ring cleavage, this information will be very helpful in order to decide the cleavage position.
		size_t ring_start;
		size_t ring_end;

		// When start and end are both 0, it just means one site available.
		FunctionalGroup(const size_t start = 0, const size_t end = 0) {}
		FunctionalGroup(const std::string& symbol, const size_t start =0, const size_t end=0)
			: _symbol(symbol), ring_start(start), ring_end(end)
		{}

		FunctionalGroup(const std::string& symbol, const std::string& name, const size_t start =0, const size_t end=0)
			: _symbol(symbol), _name(name), ring_start(start), ring_end(end)
		{}
		FunctionalGroup(const std::string& compo, const std::string& symbol, const std::string& name, const size_t start = 0, const size_t end=0)
			: Unit(compo), _symbol(symbol), _name(name), ring_start(start), ring_end(end)
		{}
		
		inline const std::string& getFunctionalGroupSymbol() const
		{
			return _symbol;
		}
		inline const std::string& getFunctionalGroupName() const
		{
			return _name;
		}
		inline void setFunctionalGroupSymbol(const std::string& symbol)
		{
			_symbol = symbol;
		}
		inline void setFunctionalGroupName(const std::string& name)
		{
			_name = name;
		}
		inline void setStatus(const std::string& status)
		{
			_status = status;
		}
		inline const std::string& getStatus() const
		{
			return _status;
		}	
		inline bool hasDependency()
		{
			return !_sites.empty();
		}
		inline std::vector<Site>& getInternalSites()
		{
			return _sites;
		}
		inline const std::vector<Site>& getInternalSites() const
		{
			return _sites;
		}
		inline Site& getInternalSiteByCarbonID(const size_t c_id)
		{
			return _sites.at(c_id);
		}
		Site& getInternalSiteByRingID(const size_t r_id)
		{
			if(r_id == 0)
				return _sites.at(0);
			else if(r_id > 0 && (ring_start + r_id - 1 <= ring_end))
				return _sites.at(ring_start + r_id - 1);
			else
				throw std::runtime_error("Unqualified ring id.");
		}

		// Convert a functional group to site. Cool function!!!
		Site getConvertedSite(const size_t idx = 0);

		// Present pointed sub-functional group divide the sub functional group from the original one. If core is true, take the functional group as the last one, otherwise, put it in the first..
		void reorganizeSites(const std::string& fg_str, const size_t idx, bool core = false);

		// Recursively explore the functional group on the sub-tree.
		bool containFunctionalGroup(const std::string& fg, const size_t idx = 0) const;
		bool containFunctionalGroup(const FunctionalGroup& fg, const size_t idx = 0) const;
		// Non-recursive way of exploring the functional group.
		bool attachFunctionalGroup(const std::string& fg_symbol, const size_t idx = 0) const;
		bool attachFunctionalGroup(const FunctionalGroup& fg, const size_t idx = 0) const;

		// Add sub-functional group into the core. 
		void addFunctionalGroup(FunctionalGroup& sub_fg, const size_t idx = 0);
		//void addFunctionalGroup(const std::string& fg_str, const size_t idx = 0);
		void addFunctionalGroupByChain(FunctionalGroupChain& chain, size_t idx = 0);

		// Remove sub-functional group from ori. The removed one should be on the terminal
		void removeFunctionalGroup(const std::string& fg_str, const size_t idx=0);
		void removeFunctionalGroup(FunctionalGroup& fg, const size_t idx = 0);
		void removeFunctionalGroupByChain(FunctionalGroupChain& chain, size_t idx = 0);
		bool containFunctionalGroupByChain(FunctionalGroupChain& chain, size_t idx = 0) const;

		void replaceFunctionalGroupByChain(FunctionalGroupChain& chain, FunctionalGroup& fg_new);

		// The pointed sub-functionalgroup will be returned for further modification.
		FunctionalGroup& getSubFunctionalGroup(const std::string& fg_str, const size_t idx = 0);
		const FunctionalGroup& getSubFunctionalGroup(const std::string& fg_str, const size_t idx = 0) const;

		Composition getCompositionByID(const size_t idx);
		
		FunctionalGroup& operator=(const FunctionalGroup& fg);

		void printFunctionalGroupInformation();

		friend bool operator<(const FunctionalGroup& left, const FunctionalGroup& right)
		{
			return left._symbol < right._symbol;
			//return &left < &right; // Pointer comparison.
		}

		inline const size_t& getRingStart() const
		{
			return ring_start;
		}
		inline const size_t& getRingEnd() const
		{
			return ring_end;
		}
		inline const std::string& getSymbol() const 
		{
			return _symbol;
		}

		inline const std::string& getName() const
		{
			return _name;
		}

		// CarbonID
		inline std::pair<size_t, size_t> getRingPosByCarbonID()
		{
			return std::pair<size_t, size_t>(ring_start, ring_end);
		}

		inline std::pair<size_t, size_t> getRingPosByRingID()
		{
			return std::pair<size_t, size_t>(1, ring_end-ring_start+1);
		}

		size_t getRingID(const size_t c_id);
		size_t getCarbonID(const size_t r_id);

		Composition getSubCompositionByCarbonID(const size_t& s1, const size_t& s2);
		Composition getSubCompositionByRingID(const size_t& s1, const size_t& s2);
		// To be compatible with addFunctionalGroup and removeFunctionalGroup functions.
		void addFunctionalGroup(size_t pos, FunctionalGroup& fg, const Composition& plus, const Composition& minus);
		
		void removeFunctionalGroup(size_t pos, FunctionalGroup& fg, const Composition& lost);

		void printStructure();
		
	};
}

#endif   /* ----- #ifndef GAG_FUNCTIONALGROUP_H_INC  ----- */
