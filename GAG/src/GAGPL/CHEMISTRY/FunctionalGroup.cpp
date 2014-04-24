#include <GAGPL/CHEMISTRY/FunctionalGroup.h>
#include <boost/type_traits.hpp>
#include <iostream>

namespace gag
{
	
	Site FunctionalGroup::getConvertedSite(const size_t idx)
	{
		// Get the first site, retrieve the core and satellite functional groups.
		// If the satellite functional group has no sites, return Site();
		return _sites.empty()? Site(): _sites.at(idx);
	}

	void FunctionalGroup::addFunctionalGroup(FunctionalGroup& sub_fg, const size_t idx /* = 0 */)
	{
		Site& temp_st = _sites.at(idx);
		temp_st.fg_map.insert(std::make_pair(sub_fg._symbol, sub_fg));
		// Update composition.
		this->add(sub_fg.getComposition());
	}

	void FunctionalGroup::removeFunctionalGroup( const std::string& fg_str, const size_t idx/*=0*/ )
	{
		Site& temp_st = _sites.at(idx);
		std::multimap<std::string, FunctionalGroup>::iterator iter = temp_st.fg_map.find(fg_str);
		if(iter != temp_st.fg_map.end()) {
			this->deduct(iter->second.getComposition());
			temp_st.fg_map.erase(iter);	
			//return true;
		} else {
			throw std::runtime_error("The to-be-removed functional group doesn't exist!");
		}
	}

	void FunctionalGroup::removeFunctionalGroup(FunctionalGroup& sub_fg, const size_t idx /* = 0 */)
	{
		this->removeFunctionalGroup(sub_fg._symbol, idx);
	}

	FunctionalGroup& FunctionalGroup::operator =(const FunctionalGroup& rhs)
	{
		if(this != &rhs){
			_symbol = rhs._symbol;
			_name = rhs._name;
			_sites = rhs._sites;
		}
		return *this;
	}

	void FunctionalGroup::printFunctionalGroupInformation()
	{
		std::cout << "Symbol: " << _symbol << std::endl;
		std::cout << "Name: " << _name << std::endl;
		std::cout << "Composition: " << this->getComposition().getCompositionString() << std::endl;

		std::cout << "Ring: Start -- " << ring_start << " End -- " << ring_end << std::endl;

		for(std::vector<Site>::iterator iter = _sites.begin(); iter != _sites.end(); iter++)
		{
			std::cout << "Site: " << std::endl;
			std::cout << "--Core: " << (*iter).core.symbol << std::endl;
			std::cout << "--FunctionalGroup:" << std::endl;;
			for(std::multimap<std::string, FunctionalGroup>::iterator ptr_iter = (*iter).fg_map.begin(); ptr_iter != (*iter).fg_map.end(); ptr_iter++)
			{
				std::cout << "  Symbol: " << ptr_iter->second.getFunctionalGroupSymbol();
				std::cout << "  Composition: " << ptr_iter->second.getCompositionString() << std::endl;
				std::cout << " Ring: Start -- " << ptr_iter->second.ring_start << " End -- " << ptr_iter->second.ring_end << std::endl;
			}
			std::cout << std::endl;
			
		}		
		std::cout << std::endl;

	}

	
	FunctionalGroup& FunctionalGroup::getSubFunctionalGroup( const std::string& fg_str, const size_t idx /*= 0*/ )
	{
		Site& temp_st = _sites.at(idx);
		std::multimap<std::string, FunctionalGroup>::iterator iter = temp_st.fg_map.find(fg_str);
		if(iter != temp_st.fg_map.end())
			return iter->second;
		else 
			throw std::runtime_error("Cannot get requested sub-functionalgroup!");
	}

	const FunctionalGroup& FunctionalGroup::getSubFunctionalGroup( const std::string& fg_str, const size_t idx /*= 0*/ ) const
	{
		const Site& temp_st = _sites.at(idx);
		std::multimap<std::string, FunctionalGroup>::const_iterator iter = temp_st.fg_map.find(fg_str);
		if(iter != temp_st.fg_map.end())
			return iter->second;
		else 
			throw std::runtime_error("Cannot get requested sub-functionalgroup!");
	}

	void FunctionalGroup::reorganizeSites(const std::string& fg_str, const size_t idx, bool core /*=false */ )
	{
		if(!this->containFunctionalGroup(fg_str, idx))
			return;
		FunctionalGroup fg(this->getSubFunctionalGroup(fg_str, idx));
		this->removeFunctionalGroup(fg_str, idx);
		Site st = fg.getConvertedSite();
		if(core) {
			fg._sites.push_back(st);			
		} else {
			std::vector<Site>::iterator iter = _sites.begin();
			_sites.insert(iter, st);
		}
	}

	bool FunctionalGroup::containFunctionalGroup(const std::string& fg_str, const size_t idx /* = 0 */) const
	{
		const Site& temp_st = _sites.at(idx);
		std::multimap<std::string, FunctionalGroup>::const_iterator iter = temp_st.fg_map.find(fg_str);
		bool status = false;
		if(iter == temp_st.fg_map.end()) { // Not found. 
			// For each site in the subsequent functional group, find check if there is any one.
			for(std::multimap<std::string, FunctionalGroup>::const_iterator iter = temp_st.fg_map.begin(); iter != temp_st.fg_map.end(); iter++){
				size_t site_idx = iter->second.getInternalSites().size();
				for(size_t i = 0; i < site_idx; i++)
				{
					if(iter->second.containFunctionalGroup(fg_str, i)) {
						status = true;
						break;
					}
				}
				if(status)
					break;
			}

		} else {
			return true;
		}

		return status;
	}

	bool FunctionalGroup::containFunctionalGroup(const FunctionalGroup& fg, const size_t idx) const
	{
		return this->containFunctionalGroup(fg._symbol, idx);
	}

	bool FunctionalGroup::attachFunctionalGroup( const std::string& fg_symbol, const size_t idx) const
	{
		const Site& temp_st = _sites.at(idx);

		return temp_st.attachFunctionalGroup(fg_symbol);
	}
	bool FunctionalGroup::attachFunctionalGroup( const FunctionalGroup& fg, const size_t idx) const
	{
		const Site& temp_st = _sites.at(idx);

		return temp_st.attachFunctionalGroup(fg._symbol);
	}


	void FunctionalGroup::addFunctionalGroupByChain(FunctionalGroupChain& chain, size_t idx /* = 0 */ )
	{
		if(idx == chain.size()-1) { // Do the modification.
			this->addFunctionalGroup(chain.at(idx).second, chain.at(idx).first);
		} else if(idx < chain.size()-1) {
			FunctionalGroup& sub_fg = this->getSubFunctionalGroup(chain.at(idx).second.getFunctionalGroupSymbol(), chain.at(idx).first);
			sub_fg.addFunctionalGroupByChain(chain, idx+1);
			// *** Recursively update the composition. ***
			compo.add(chain.back().second.getComposition());
		} else if(idx > chain.size()-1) {
			throw std::runtime_error("Index of the chain is illegal!");
		}

	}

	void FunctionalGroup::removeFunctionalGroupByChain(FunctionalGroupChain& chain, size_t idx /*= 0*/ )
	{
		if(idx == chain.size()-1) { // Do the modification.
			this->removeFunctionalGroup(chain.at(idx).second, chain.at(idx).first);
		} else if(idx < chain.size()-1) {
			FunctionalGroup& sub_fg = this->getSubFunctionalGroup(chain.at(idx).second.getFunctionalGroupSymbol(), chain.at(idx).first);
			sub_fg.removeFunctionalGroupByChain(chain, idx+1);
			// *** Recursively update the composition. ***
			compo.deduct(chain.back().second.getComposition());
		} else if(idx > chain.size()-1) {
			throw std::runtime_error("Index of the chain is illegal!");
		}
	}


	bool FunctionalGroup::containFunctionalGroupByChain( FunctionalGroupChain& chain, size_t idx /*= 0*/ ) const
	{
		if(idx == chain.size()-1) { // Do the modification.
			return this->attachFunctionalGroup(chain.at(idx).second, chain.at(idx).first);
		} else if(idx < chain.size()-1) {
			const FunctionalGroup& sub_fg = this->getSubFunctionalGroup(chain.at(idx).second.getFunctionalGroupSymbol(), chain.at(idx).first);
			return sub_fg.containFunctionalGroupByChain(chain, idx+1);
		} else {
			throw std::runtime_error("Index of the chain is illegal!");
			return false;
		}
	}

	void FunctionalGroup::replaceFunctionalGroupByChain( FunctionalGroupChain& chain, FunctionalGroup& fg_new )
	{
		this->removeFunctionalGroupByChain(chain);
		// Update chain information.
		size_t position = chain.back().first;
		chain.pop_back();
		chain.push_back(std::make_pair(position, fg_new));
		
		this->addFunctionalGroupByChain(chain);
	}

	Composition FunctionalGroup::getCompositionByID( const size_t idx )
	{
		if(_sites.size() != 1) {
			Site& temp_st = _sites.at(idx);
			return temp_st.getSiteComposition();
		} else {
			this->reorganizeSites("H",0,true);
			return _sites.at(0).getSiteComposition();
		}
	}

	size_t FunctionalGroup::getRingID( const size_t c_id )
	{
		if(c_id == 0)
			return 0;
		else if(c_id <= ring_start)
			return 1;
		else if(c_id <= ring_end)
			return (c_id - ring_start + 1);
		else if(c_id > ring_end)
			return (ring_end - ring_start + 1);
		else 
			throw std::runtime_error("Unqualified carbon ID");
	}

	size_t FunctionalGroup::getCarbonID( const size_t r_id )
	{
		if(r_id == 0)
			return 0;
		else if(r_id >= 1 && (ring_start + r_id -1 <= ring_end))
			return (ring_start + r_id -1);
		else if(ring_start + r_id -1 > ring_end)
			throw std::runtime_error("Unqualified ring ID");

		return 9999;
	}

	Composition FunctionalGroup::getSubCompositionByCarbonID( const size_t& s1, const size_t& s2 )
	{
		Composition chain_compo;

		if(s1 > s2)
			throw std::runtime_error("Unqualified carbon ID pair!");

		for(std::vector<Site>::iterator iter = _sites.begin() + s1; iter!= _sites.begin() + s2 + 1; iter++)
		{
			chain_compo.add(iter->getSiteComposition());
		}

		return chain_compo;
	}

	Composition FunctionalGroup::getSubCompositionByRingID( const size_t& s1, const size_t& s2 )
	{
		if(s1 > s2 || ring_start + s2 -1 > ring_end)
			throw std::runtime_error("Unqualified ring ID");

		Composition ring_compo;
		if(s1 == 0 && s2 == 0) // The whole ring. 
			ring_compo = (*this).getSubCompositionByCarbonID(0, _sites.size()-1);
		else if(s1 == 0 && ring_start+s2-1 == ring_end) 
			ring_compo = (*this).getSubCompositionByCarbonID(1, _sites.size()-1);
		else if(s1 == 0 && ring_start+s2-1 < ring_end) // Start from Carbon ID 1.
			ring_compo = (*this).getSubCompositionByCarbonID(1, ring_start + s2 -1);
		else if(ring_start+s2-1 < ring_end)
			ring_compo = (*this).getSubCompositionByCarbonID(ring_start+s1, ring_start+s2-1);
		else if(ring_start+s2-1 == ring_end)
			ring_compo = (*this).getSubCompositionByCarbonID(ring_start+s1, _sites.size()-1);
		else 
			throw std::runtime_error("Strange things happened!");
		//ring_compo.updateString();
		return ring_compo;
	}

	
	void FunctionalGroup::addFunctionalGroup( size_t pos, FunctionalGroup& fg, const Composition& plus, const Composition& minus )
	{
		std::multimap<std::string, FunctionalGroup>::iterator iter = _sites.at(pos).fg_map.find(fg.getSymbol());
		if(iter != _sites.at(pos).fg_map.end())
		{
			(*iter).second.add(plus);
			(*iter).second.deduct(minus);

			compo.add(plus);
			compo.deduct(minus);
		}
	}

	void FunctionalGroup::removeFunctionalGroup( size_t pos, FunctionalGroup& fg, const Composition& lost )
	{
		Site& st = _sites.at(pos);
		std::multimap<std::string, FunctionalGroup>::iterator iter = st.fg_map.find(fg.getSymbol());
		
		if(iter != st.fg_map.end())
		{
			if(!(lost == (*iter).second.getComposition())){
				(*iter).second.deduct(lost);
				//_sites.at(pos).deduct(lost);
				compo.deduct(lost);
			}
			else {
				st.fg_map.erase(iter);
				//_sites.at(pos).deduct(fg.getComposition());
				compo.deduct(fg.getComposition());
			}
		}
	}

	void FunctionalGroup::printStructure()
	{
		std::cout << "Name: " << this->getSymbol() << std::endl;
		std::cout << "Monosaccharide Composition: " << this->getCompositionString() << std::endl;
		std::vector<Site>& sites = this->getInternalSites();

		for(std::vector<Site>::iterator iter = sites.begin(); iter != sites.end(); iter++)
		{
			iter->printStructure();
		}
	}

	Composition Site::getSiteComposition()
	{
		Composition temp_compo(core.symbol);
		std::multimap<std::string, FunctionalGroup>::iterator iter = fg_map.begin();
		for(; iter != fg_map.end(); iter++)
		{
			temp_compo.add(iter->second.getComposition());
		}
		return temp_compo;
	}
	std::string Site::getSiteCompositionString()
	{
		return this->getSiteComposition().getCompositionString();
	}

	void Site::printStructure()
	{
		std::cout << "Site Composition: " << this->getSiteCompositionString() << std::endl;
		std::multimap<std::string, FunctionalGroup>& fg = this->fg_map;
		for(std::multimap<std::string, FunctionalGroup>::iterator iter = fg.begin(); iter != fg.end(); iter++)
		{
			std::cout << (*iter).second.getCompositionString() << std::endl;
		}
	}

}