/********************************************************************
	created:	2012/06/26
	created:	26:6:2012   15:12
	filename: 	Modification.h
	file path:	GAG\src\GAGPL\CHEMISTRY
	file base:	Modification
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef GAG_MODIFICATION_H
#define GAG_MODIFICATION_H

#include <set>
#include <map>
#include <vector>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

namespace gag
{
	// Behavior of operation on the atom.
	enum AtomOperation {Addition, Removement, N};

	// The core index, operation and corresponding composition.
	typedef std::vector<boost::tuple<size_t, AtomOperation, std::string> > OperationSet;

	// The operation and corresponding composition. This will work on the given core.
	typedef std::vector<std::pair<AtomOperation, std::string> > PairSet;
	
	// The class contains the modification location and the exact atom on the location. It is useful when dealing with N-sulfation/O-sulfation.
	struct MacroPosition
	{
		size_t branch_id;
		size_t mono_id;

    //MacroPosition() {}

		MacroPosition(const size_t& b_id, const size_t& m_id)
			: branch_id(b_id), mono_id(m_id)
		{}

		inline bool operator<(const MacroPosition& mac_pos) const
		{
			return branch_id < mac_pos.branch_id || (branch_id == mac_pos.branch_id && mono_id < mac_pos.mono_id);
		}

		inline bool operator==(const MacroPosition& mac_pos) const
		{
			return (branch_id == mac_pos.branch_id && mono_id == mac_pos.mono_id);
		}
		inline bool operator!=(const MacroPosition& mac_pos) const
		{
			return (branch_id != mac_pos.branch_id || mono_id != mac_pos.mono_id);
		}
	};

	struct ModificationPosition
	{	
		MacroPosition macro_pos;

		size_t site_id;

		std::string atom;

		inline bool operator==(const ModificationPosition& mod) const
		{
			return (this->macro_pos.branch_id == mod.macro_pos.branch_id && this->macro_pos.mono_id == mod.macro_pos.mono_id && this->site_id == mod.site_id && this->atom == mod.atom);
		}

		// The sorting order is branch id, mono id, site id, and atom in the last.
		inline bool operator<(const ModificationPosition& mod) const
		{
			return (macro_pos < mod.macro_pos || (macro_pos == mod.macro_pos && site_id < mod.site_id) || (macro_pos == mod.macro_pos && site_id == mod.site_id && atom < mod.atom));
		}

    //ModificationPosition& operator=(const ModificationPosition& mod_pos)
    //{
    //  macro_pos = mod_pos.macro_pos;
    //  site_id = mod_pos.site_id;
    //  atom = mod_pos.atom;
    //  return *this;
    //}

		inline size_t getBranchID() const
		{
			return macro_pos.branch_id;
		}
		inline size_t getMonosaccharideID() const
		{
			return macro_pos.mono_id;
		}

		// Only use the first and second id information.
		ModificationPosition(size_t b_id, size_t m_id)
			: macro_pos(b_id, m_id)
		{}

		ModificationPosition(size_t b_id, size_t m_id, size_t s_id, std::string a)
			: macro_pos(b_id, m_id), site_id(s_id), atom(a)
		{}
		
		ModificationPosition(size_t b_id, size_t m_id, size_t s_id)
			: macro_pos(b_id, m_id), site_id(s_id)
		{}

    /*const ModificationPosition& operator=(const ModificationPosition& mod_pos)
    {
      macro_pos = mod_pos.macro_pos;
      atom = mod_pos.atom;
      site_id = mod_pos.site_id;
      return *this;
    }*/

		std::string printString() const;

	};

	typedef std::set<ModificationPosition> ModificationSites;
  typedef std::map<std::string, ModificationSites> ModificationSequence;

  // A - B.
	ModificationSites getSiteDifference(const ModificationSites& m1, const ModificationSites& m2);
	ModificationSites getSiteIntersection(const ModificationSites& m1, const ModificationSites& m2);
  // If ms1 includes ms2.
    bool containSubset(const ModificationSites& ms1, const ModificationSites& ms2);

	void printModificationSites(const ModificationSites& mod_sites);
    std::ostream& operator<<(std::ostream& os, const ModificationSites& mod_sites);

  std::string modificationString(const ModificationSites& mod_sites);

	struct Reaction
	{
		// Multiple core.
		size_t position;

		// Operation on the core. There is no position information.
		PairSet core_operation;

		// Operation on the CORE of the SUB functional group. It will ONLY apply if there is such kind of functional group.
		std::multimap<std::string, OperationSet> sub_fg_operation;

		// Operation on the reactant.
		std::pair<std::string, OperationSet> reactant_operation;
	};

	// Functional group name (usually monosaccharide unit) and modification sites.
	typedef std::map<std::string, std::vector<size_t> > ModificationRule;

	class Modification
	{
	private:
		std::string mod_name;
		
		// The reason that I use multiple shortcuts (symbols) is that I need to parse different types of user input.
		std::set<std::string> symbols;
		
		std::vector<Reaction> mod_steps;

		// Functional group name and modification sites.
		ModificationRule mod_rule; 

	public:
		//void addOperation(GlycanSequence& gs, Operation& op);

		inline const std::string& getSymbol()
		{
			return *(symbols.begin());
		}
		
		inline const std::string& getModificationName() const
		{
			return mod_name;
		}

		inline const std::vector<Reaction>& getModificationReactions() const
		{
			return mod_steps;
		}
		inline std::vector<Reaction>& getModificationReactions()
		{
			return mod_steps;
		}

		inline std::set<std::string>& getSymbols()
		{
			return symbols;
		}

		inline void addReaction(const Reaction& reaction)
		{
			mod_steps.push_back(reaction);
		}

		inline void addModificationRule(const std::string& fg_symbol, const std::vector<size_t>& sites)
		{
			mod_rule.insert(std::make_pair(fg_symbol, sites));
		}

		std::vector<size_t> getModificationSites(const std::string& fg_symbol);
		// Apply the reactions storing in mod_steps to the input functional group.
		// void applyReactions(FunctionalGroup& fg);

		// Constructor.
		Modification() {}
		Modification(const std::string& name, const std::set<std::string>& sc_set, const std::vector<Reaction>& reactions, const ModificationRule& rules)
			: mod_name(name), symbols(sc_set), mod_steps(reactions), mod_rule(rules)
		{}
		Modification(const std::string& name)
			: mod_name(name)
		{}

	};

	
	// Recording the status of modification sites on the sequence.
	
	using namespace ::boost;
	using namespace ::boost::multi_index;

	struct ModificationInfo
	{
		// The id of modification position should be controlled by createModificationPosition.
		ModificationPosition mod_pos;
		std::string mod_symbol;
		int mod_status;

		inline bool operator<(const ModificationInfo& mod_info) const
		{
			return (mod_pos < mod_info.mod_pos || (mod_pos == mod_info.mod_pos && mod_symbol < mod_info.mod_symbol));
		}

		inline size_t getBranchID() const
		{
			return mod_pos.macro_pos.branch_id;
		}

		//ModificationInfo() {}
		ModificationInfo(const ModificationPosition& pos, const std::string& mod_symbol, int status)
			: mod_pos(pos), mod_symbol(mod_symbol), mod_status(status) 
		{}

	};

	// It is better to use pointer, since the information between the sequence and fragments need to be synchronized.
	//typedef boost::make_shared<ModificationInfo> ModificationInfoPtr;

	struct mod_symbol{};
	struct mod_status{};
	struct mod_position{};
	struct mod_bid{};
	struct mod_compo_key{};

	typedef multi_index_container<
		ModificationInfo,
		indexed_by<
			ordered_non_unique< // There might be multiple modifications mapped to 
			// one site.
			tag<mod_position>, member<ModificationInfo, ModificationPosition, &ModificationInfo::mod_pos>
			>,
			ordered_non_unique<
			tag<mod_bid>, const_mem_fun<ModificationInfo, size_t, &ModificationInfo::getBranchID>
			>,
			ordered_non_unique<
				tag<mod_compo_key>,
				composite_key<
					ModificationInfo,
					member<ModificationInfo, std::string, &ModificationInfo::mod_symbol>,
					member<ModificationInfo, int, &ModificationInfo::mod_status>
				> 
			>
		>
	> ModificationContainer;

	//typedef ModificationContainer::index<mod_symbol>::type mod_by_symbol;
	typedef ModificationContainer::index<mod_position>::type ModificationByPosition;
	typedef ModificationContainer::index<mod_bid>::type ModificationByBranchID;
	typedef ModificationContainer::index<mod_compo_key>::type ModificationByCompositeKey;

}


#endif /* GAG_MODIFICATION_H */