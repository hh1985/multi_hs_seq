/********************************************************************
	created:	2013/01/23
	created:	23:1:2013   13:32
	filename: 	Modifier.h
	file path:	GAG\src\GAGPL\CHEMISTRY
	file base:	Modifier
	file ext:	h
	author:		Han Hu
	
	purpose:	Higher level of managing modification stuff. Especially the
	          communication between mod table and fg table. 
*********************************************************************/

#ifndef GAG_MODIFIER_H
#define GAG_MODIFIER_H

#include <GAGPL/CHEMISTRY/ModificationTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/GLYCAN/Monosaccharide.h>
//#include <GAGPL/GLYCAN/GlycanSequence.h>

namespace gag
{
	class Modifier
	{
	public:
		

		//ModificationMap mod_map;

		Modifier()
			: mod_table(ModificationTable::Instance()), fun_table(FunctionalGroupTable::Instance())
		{}

		Modifier(Modifier& modifier)
			: mod_table(ModificationTable::Instance()), fun_table(FunctionalGroupTable::Instance()), mod_pool(modifier.getModificationContainer())
		{}

		// The modification with self-owned site position.
		void modifyFunctionalGroup(FunctionalGroup& fg, const std::string& mod_symbol);
		void modifyFunctionalGroup(FunctionalGroup& fg, const Modification& mod);

		// The modification with user-specified site position.
		void modifyFunctionalGroup(FunctionalGroup& fg, const std::string& mod_symbol, const size_t index);

		// The function is used to estimate the potential of modification, but do not operation on the functional group, if the site has been occupied, the modification is not appliable.
		bool isFunctionalGroupModifiable(const FunctionalGroup& fg, Modification& mod);
		bool isFunctionalGroupModifiable(const FunctionalGroup& fg, const std::string& mod_symbol, const size_t index);

		// The function doesn't consider about the occupation of modification site.
		std::vector<size_t> getModificationSetByPosition(const std::string& fg_symbol, const std::string& mod_symbol);
		// The function is able to deal with the modified functional group.
		std::vector<size_t> getModificationSetByPosition(const FunctionalGroup& fg, const std::string& mod_symbol);
		size_t getModificationSiteNum(const std::string& mod_symbol, int flag) const;
		size_t getModificationSiteNum(const std::string& mod_symbol) const;

		inline size_t getSize() const
		{
			return mod_pool.size();
		}

		void modifyModificationStatus( const std::string& mod_symbol, const ModificationPosition& pos, int status = 0);
		void modifyModificationStatus(ModificationByPosition::iterator pos_iter, int new_status);
		// Search the modification information on specified ring.
		std::pair<ModificationByPosition::iterator, ModificationByPosition::iterator> getModificationSiteIDRange(const MacroPosition& mac_pos);

		inline void addModification(const ModificationPosition& pos, const std::string& mod_symbol, int status)
		{
			mod_pool.insert(ModificationInfo(pos, mod_symbol, status));
		}

		inline ModificationContainer& getModificationContainer()
		{
			return mod_pool;
		}

		ModificationSites getModificationSitesBySymbol(const std::string& mod_symbol, int flag);
		ModificationSites getModificationSitesBySymbol(const std::string& mod_symbol);
	
		int getModificationStatusByPosition(const ModificationPosition& mod_pos, const std::string& mod_symbol);

		ModificationContainer mod_pool;
	private:
		ModificationTable& mod_table;
		FunctionalGroupTable& fun_table;
	};
}



#endif /* GAG_MODIFIER_H */