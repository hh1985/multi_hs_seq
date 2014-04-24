#include <GAGPL/GLYCAN/InternalSite.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>

namespace gag
{
	void InternalSite::replace(FunctionalGroup& ori, FunctionalGroup& re)
	{
		std::multiset<FunctionalGroup>::iterator iter = _groups.find(ori);
		if(iter != _groups.end())
			*iter = re;
		// update InternalSite composition.
		compo.deduct(ori.getComposition());
		compo.add(re.getComposition());
	}

	Composition InternalSite::getCleavageShift(FunctionalGroup& fg, const size_t idx)
	{
		std::multiset<FunctionalGroup>::iterator iter = _groups.find(fg);
		return (iter != _groups.end()) ? iter->getCompositionByID(idx) : Composition();
	}

	Composition InternalSite::getCleavageShift(const size_t idx)
	{
		std::multiset<FunctionalGroup>::iterator iter = _groups.begin();
		Composition temp_compo;
		for(; iter != _groups.end(); iter++)
		{
			if(iter->getCompositionString() != "H") {
				temp_compo = iter->getCompositionByID(idx);
				break;
			}

		}
		return temp_compo;
	}

	bool InternalSite::hasFunctionalGroup( FunctionalGroup& fg )
	{
		std::multiset<FunctionalGroup>::iterator iter = _groups.find(fg);
		if(iter != _groups.end())
			return true;
		else 
			return false;
	}

	bool InternalSite::hasFunctionalGroup( const std::string& str )
	{
		FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
		FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol(str);
		return this->hasFunctionalGroup(fg);
	}

	void InternalSite::addFunctionalGroup( const std::string& str )
	{
		FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
		FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol(str);
		_groups.insert(fg);
	}

	void InternalSite::addFunctionalGroup(std::vector<FunctionalGroup>& fg_chain, FunctionalGroup& fg)
	{
		if(fg_chain.empty()) {
			this->addFunctionalGroup(fg);
		} else {
			std::multiset<FunctionalGroup>::iterator iter = _groups.find(fg_chain.front());
			if(iter != _groups.end()) {
				FunctionalGroup& next_fg = *iter;
				for(size_t i=1; i<fg_chain.size(); i++) {
					next_fg = next_fg.getSubFunctionalGroup(fg_chain.at(i)._symbol);
				}
				next_fg.addFunctionalGroup(fg);
			} else {
				// Throw an exception.
			}
		}
	}

	void InternalSite::removeFunctionalGroup(std::vector<FunctionalGroup>& fg_chain)
	{
		std::multiset<FunctionalGroup>::iterator iter = _groups.find(fg_chain.front());
		if(iter != _groups.end()) {
			FunctionalGroup& next_fg = *iter;
			FunctionalGroup& last_fg = *iter;
			for(size_t i=1; i<fg_chain.size()-1; i++) {
				next_fg = last_fg.getSubFunctionalGroup(fg_chain.at(i)._symbol);
			}
			last_fg.removeFunctionalGroup(next_fg);
		} else {
			// Throw an exception.
		}
	}

}