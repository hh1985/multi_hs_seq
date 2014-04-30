/*
 * =====================================================================================
 *
 *       Filename:  Fragmentation.cpp
 *
 *    Description:  The implement file for class Fragment.
 *
 *        Version:  1.0
 *        Created:  05/ 8/2012  2:24:18 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <algorithm>
#include <set>
#include <string>
#include <sstream>
#include <iostream>

namespace gag
{
  // if true, simply remove the corresponding part instead of recalculation.
	// Only works for A type.
	bool Fragment::isInternalCleavage(const FragmentPosition& fp)
	{
		Branch& bc = glyco_seq->getBranchByID(fp.getBranchID());
		if(fp.getMonosaccharideID() != 0) { // Within branch.
			if(bc.getLinkages().size() < fp.getMonosaccharideID()+1) // Reducing end;
				return true;
			if((bc.getLinkages().at(fp.getMonosaccharideID()-1).end > bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_second)&& bc.getLinkages().at(fp.getMonosaccharideID()).start <= bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_first)) ||
				(bc.getLinkages().at(fp.getMonosaccharideID()-1).end <= bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_second)&& bc.getLinkages().at(fp.getMonosaccharideID()).start > bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_first)))
				return true;
			else 
				return false;
		} else {
			BranchMap::right_const_iterator right_iter = glyco_seq->getBranchLinks().right.find(fp.getBranchID());
			if(right_iter != glyco_seq->getBranchLinks().right.end()){ // No leaf
				if((glyco_seq->getBranchByID(fp.getBranchID()-1).getLinkages().back().end > bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_second) && bc.getLinkages().at(fp.getMonosaccharideID()).start <= bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_first)) ||
					(glyco_seq->getBranchByID(fp.getBranchID()-1).getLinkages().back().end <= bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_second) && bc.getLinkages().at(fp.getMonosaccharideID()).start > bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_first)))
					return true;
				else 
					return false;
			} else { // Leaf.
				if(bc.getLinkages().at(fp.getMonosaccharideID()).start > bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_first) && 
					bc.getLinkages().at(fp.getMonosaccharideID()).start < bc.getUnitByID(fp.getMonosaccharideID()).getCarbonID(fp.xring_second))
					return true;
				else 
					return false;
			}
		}

	}
	// Record cleavage type and fragmentation position.
	// 1. Each time when the function is called, the validity of the cleavage should be checked.
	// 2. The composition should be updated immediately.
	void Fragment::setFragmentation(const std::string& type, const FragmentPosition& site)
	{
		//CleavageCollection::iterator iter1 = cleavage_sites.find(type);
		//if(iter1 != cleavage_sites.end()) {// The type of cleavage has been created.	
		//	(*iter1).second.push_back(site);
		//} else {
		//	std::vector<FragmentPosition> site_vec;
		//	site_vec.push_back(site);
		//	cleavage_sites.insert(std::make_pair(type, site_vec));
		//}
		cleavage_sites.insert(std::make_pair(type, site));
		
		/* Update the composition. */
		Branch& bc = glyco_seq->getBranchByID(site.getBranchID());
		Composition ori_compo = glyco_seq->getComposition();	
		Composition tmp_compo2 = this->getFragmentComposition(site, false);
		Composition tmp_compo1 = this->getFragmentComposition(site);

		if(type == "Y" || type == "Z") { // The RE part will be kept.
			compo.deduct(tmp_compo1); 
			this->cleanModificationStatus(site);
			if(type == "Y") compo.add("O");
		} else if(type == "B" || type == "C") {				// The NRE part will be kept.
			ori_compo.deduct(tmp_compo1);
			compo.deduct(ori_compo);
			this->cleanModificationStatus(site, false);
			if(type == "B") compo.deduct("O");
		} else if(type == "X") {
			this->cleanModificationStatus(site);
			if(bc.getUnitByID(site.getMonosaccharideID()).getCarbonID(site.xring_first) < bc.getUnitByID(site.getMonosaccharideID()).getRingStart()){  // 0-2
				compo.deduct(tmp_compo2);
				//this->cleanModificationStatus(site, false);
			} else { // 1-5.
				compo.deduct(tmp_compo1);
				//this->cleanModificationStatus(site);
			}
		} else if(type == "A") {
			this->cleanModificationStatus(site, false);
			if(bc.getUnitByID(site.getMonosaccharideID()).getCarbonID(site.xring_first) < bc.getUnitByID(site.getMonosaccharideID()).getRingStart()){
				ori_compo.deduct(tmp_compo2);
				compo.deduct(ori_compo);
				
			} else {
				ori_compo.deduct(tmp_compo1);
				compo.deduct(ori_compo);
				//this->cleanModificationStatus(site);
			}
		} 
		
		////FragmentationParams& fp = ft.getFragmentationParams(type);
		FragmentationTable& ft = FragmentationTable::Instance();
		CompositionShift& cs = ft.getCleavageShift(type); 

		// Just mechanically add or remove the composition shift at this moment.
		for(CompositionShift::iterator iter = cs.begin(); iter != cs.end(); iter++) {
			if(!iter->second.empty()) {
				if(iter->first == 1)
					compo.add(iter->second);
				else if(iter->first == -1)
					if(!(compo < iter->second))
						compo.deduct(iter->second);
			}
		}
	}

	void Fragment::setFragmentation( const std::string& type, const size_t& index, const std::string leaf, const size_t& xring_first, const size_t& xring_second )
	{
		if(type == "A" || type == "B" || type == "C")
			this->setFragmentation(type, FragmentPosition(0, index-1, xring_first, xring_second));
		else if(type == "X" || type == "Y" || type == "Z")
			this->setFragmentation(type, FragmentPosition(0, glyco_seq->getBranchByID(0).getUnitNum()-index-1, xring_first, xring_second));
		else {
			throw std::runtime_error("The cleavage is not covered yet!");
		}
	}

  void Fragment::setFragmentation( const std::string& type, const size_t& index, const size_t& xring_first, const size_t& xring_second )
  {
    if(type == "A" || type == "B" || type == "C")
      this->setFragmentation(type, FragmentPosition(0, index-1, xring_first, xring_second));
    else if(type == "X" || type == "Y" || type == "Z")
      this->setFragmentation(type, FragmentPosition(0, glyco_seq->getBranchByID(0).getUnitNum()-index-1, xring_first, xring_second));
    else {
      throw std::runtime_error("The cleavage is not covered yet!");
    }
  }

	// Generate generic segment from NRE for each fragmentation site.
	// The start id can be and end id can be the 
	Composition Fragment::getFragmentComposition(const FragmentPosition& fp, bool cw)
	{
		Branch& bc = glyco_seq->getBranchByID(fp.getBranchID());
		Monosaccharide& ms = bc.getUnitByID(fp.getMonosaccharideID());

		Composition tmp_compo;
		Composition& cp1 = ms.getSubCompositionByRingID(fp.xring_first, fp.xring_second);
		if(cw == true){
			tmp_compo.add(cp1);
			if(fp.getMonosaccharideID() != 0) {
				
				if(bc.getLinkages().at(fp.getMonosaccharideID()-1).end <= ms.getCarbonID(fp.xring_second) || fp.xring_second == 0) {
					tmp_compo.add(glyco_seq->getSubTreeComposition(fp.getBranchID()));					tmp_compo.add(bc.getSubComposition(0,fp.getMonosaccharideID()-1));
				}
				
			} else {
				BranchMap::right_const_iterator right_iter = glyco_seq->getBranchLinks().right.find(fp.getBranchID());
				BOOST_FOREACH(BranchMap::right_reference right_ref, glyco_seq->getBranchLinks().right.equal_range(fp.getBranchID()))
				{
					size_t hinge_id = glyco_seq->getBranchByID(right_ref.second).getLinkages().back().end;

					if(hinge_id > ms.getCarbonID(fp.xring_first) && hinge_id <= ms.getCarbonID(fp.xring_second))
						tmp_compo.add(glyco_seq->getTreeComposition(right_ref.second));
				}
			}
		} else { // Anti-clockwise
			Composition mono_compo = ms.getComposition();
			tmp_compo.add(mono_compo);
			tmp_compo.deduct(cp1);
			if(fp.getMonosaccharideID() != 0){
				
				if(bc.getLinkages().at(fp.getMonosaccharideID()-1).end > ms.getCarbonID(fp.xring_second)) {
					tmp_compo.add(glyco_seq->getSubTreeComposition(fp.getBranchID()));
					tmp_compo.add(bc.getSubComposition(0,fp.getMonosaccharideID()-1));
				}
			} else {
				BranchMap::right_const_iterator right_iter = glyco_seq->getBranchLinks().right.find(fp.getBranchID());
				BOOST_FOREACH(BranchMap::right_reference right_ref, glyco_seq->getBranchLinks().right.equal_range(fp.getBranchID()))
				{
					size_t hinge_id = glyco_seq->getBranchByID(right_ref.second).getLinkages().back().end;

					if(hinge_id > ms.getCarbonID(fp.xring_second))
						tmp_compo.add(glyco_seq->getTreeComposition(right_ref.second));
				}
			}
		}

		return tmp_compo;
	}

	Fragment& Fragment::operator=(const Fragment& rhs)
	{
		if(this != &rhs)
		{
			cleavage_sites = rhs.cleavage_sites;
			glyco_seq = rhs.glyco_seq;
			compo = rhs.compo;
			mod_pool = rhs.mod_pool;
		}
		return *this;
	}

	Fragment::Fragment( const Fragment& fg )
		: param(Param::Instance())
	{
		cleavage_sites = fg.cleavage_sites;
		glyco_seq = fg.glyco_seq;
		compo = fg.compo;
		mod_pool = fg.mod_pool;
	}

	void Fragment::cleanModificationStatus( const FragmentPosition& fp, bool cw /*= true*/ )
	{
		Branch& bc = glyco_seq->getBranchByID(fp.getBranchID());
		Monosaccharide& ms = bc.getUnitByID(fp.getMonosaccharideID());
		
		if(cw == true) { // Clean the NRE part: X, Y, Z.
			if(fp.xring_first == 0 && fp.xring_second == 0) {
				this->cleanChainModificationStatus(fp.getBranchID(), 0, fp.getMonosaccharideID());
				return;
			}

			if(fp.getMonosaccharideID() != 0)
				this->cleanChainModificationStatus(fp.getBranchID(), 0, fp.getMonosaccharideID()-1);

			if(ms.getCarbonID(fp.xring_first) < ms.getRingStart()) { //0-X cleavage
				this->cleanSubUnitModificationStatusByRingID(fp.mac_pos, fp.xring_second, fp.xring_first);
			} else {
				this->cleanSubUnitModificationStatusByRingID(fp.mac_pos, fp.xring_first, fp.xring_second);
			}

		} else { // Clean the RE part: A, B, C
			if(fp.getMonosaccharideID() != bc.getUnitNum()-1)
				this->cleanChainModificationStatus(fp.getBranchID(), fp.getMonosaccharideID()+1, bc.getUnitNum()-1);
			
			if(fp.xring_first == 0 && fp.xring_second == 0)
				return;
			
			if(ms.getCarbonID(fp.xring_first) < ms.getRingStart()) { //0-X cleavage
				this->cleanSubUnitModificationStatusByRingID(fp.mac_pos, fp.xring_first, fp.xring_second);
			} else {
				this->cleanSubUnitModificationStatusByRingID(fp.mac_pos, fp.xring_second, fp.xring_first);
			}
		}
	}

	void Fragment::cleanChainModificationStatus(const size_t& branch_id, const size_t& m1, const size_t& m2 )
	{
		// clean each branch.
		for(size_t i = m1; i < m2+1; i++) {
			std::pair<ModificationByPosition::iterator, ModificationByPosition::iterator> p = this->getModificationSiteIDRange(MacroPosition(branch_id, i));
			while(p.first != p.second) {
				this->modifyModificationStatus(p.first->mod_symbol, p.first->mod_pos, -1);
				p.first++;
			}
		}
	}

	void Fragment::cleanTreeModificationStatus(const size_t& branch_id)
	{
		std::pair<size_t, size_t>& range = glyco_seq->getDescendantBranchIDs(branch_id);
		this->cleanChainModificationStatus(branch_id, range.first, range.second+1);
	}
	void Fragment::cleanSubTreeModificationStatus(const size_t& branch_id)
	{
		std::pair<size_t, size_t>& range = glyco_seq->getDescendantBranchIDs(branch_id);
		this->cleanChainModificationStatus(branch_id, range.first, range.second);
	}

	void Fragment::cleanSubUnitModificationStatusByRingID( const MacroPosition& mac_pos, const size_t& s1, const size_t& s2 )
	{
		// Special notation. No clean
		if(s1 == 0 && s2 == 0) {
			return;
		}

		Monosaccharide& mono = glyco_seq->getBranchByID(mac_pos.branch_id).getUnitByID(mac_pos.mono_id);
		const size_t first = std::min(s1, s2); const size_t second = std::max(s1, s2);

		if(mono.ring_start + second -1 > mono.ring_end || (mono.ring_start + second -1 == mono.ring_end && first == 0)) {
			throw std::runtime_error("Unqualified ring ID");
			return;
		}

		// Convert ring id to carbon id.
		const size_t c0 = first == 0 ? 1 : mono.ring_start + first;
		const size_t c1 = mono.ring_start+second-1 == mono.ring_end ? mono._sites.size()-1 : mono.ring_start+second-1;

		if(s1 < s2) {
			this->cleanSubUnitModificationStatusByCarbonID(mac_pos, c0, c1);
			return;
		} else if (s1 > s2){
			// clean the rest of c0 to c1.
			if(mono._sites.size() - c1 > 0)
				this->cleanSubUnitModificationStatusByCarbonID(mac_pos, c1+1, mono._sites.size()-1);
			if(c0 > 1)
				this->cleanSubUnitModificationStatusByCarbonID(mac_pos, 1, c0-1);
		} else {
			throw std::runtime_error("Strange cleavage!!!");
		}
		
	}

	void Fragment::cleanSubUnitModificationStatusByCarbonID( const MacroPosition& mac_pos, const size_t& c1, const size_t& c2 )
	{
		Monosaccharide& mono = glyco_seq->getBranchByID(mac_pos.branch_id).getUnitByID(mac_pos.mono_id);
		
		ModificationByPosition& mod_pos_index = this->getModificationContainer().get<mod_position>();
		for(size_t i = c1; i <= c2; i++)
		{
			std::pair<ModificationByPosition::iterator, ModificationByPosition::iterator> p = mod_pos_index.equal_range(ModificationPosition(mac_pos.branch_id, mac_pos.mono_id, i));
			while(p.first != p.second)
			{
				this->modifyModificationStatus(p.first, -1);
				++p.first;
			}
		}
	}

	std::string Fragment::getCleavageType() const
	{
		std::string cleavage_str;
		if(cleavage_sites.size() == 0)
			return "SEQ";

		for(CleavageCollection::const_iterator iter = cleavage_sites.begin(); iter!=cleavage_sites.end();iter++)
		{
			cleavage_str.append(iter->first);
			const FragmentPosition& pos = iter->second;
			/*BOOST_FOREACH(const FragmentPosition& pos, iter->second)
			{*/
				size_t id = pos.getMonosaccharideID()+1;
				if(iter->first == "A" || iter->first == "B" || iter->first == "C") {
				} else 
					id = glyco_seq->getBranchByID(0).getUnitNum() - id;

				std::ostringstream ostr;

				// Xring-cleavage.
				if(iter->first == "A" || iter->first == "X") {
					ostr << id << ":" << pos.xring_first << "-" << pos.xring_second;
					cleavage_str.append(ostr.str());
				} else {
					ostr << id;
					cleavage_str.append(ostr.str());
				}
			//}
			
		}
		return cleavage_str;
	}

	void Fragment::printFragment() const
	{
		std::cout << "General type: " << getGeneralType() << std::endl;
		std::cout << "Cleavage type: " << getCleavageType() << std::endl;
		std::cout << "Composition: " << this->getCompositionString() << std::endl;
		std::cout << "Mass: " << this->getMass() << std::endl;
		std::cout << "Modification: ";

		std::set<std::string> mod_types = glyco_seq->getModificationTypes();
		for(std::set<std::string>::iterator iter = mod_types.begin(); iter != mod_types.end(); iter++)
		{
			std::cout << *iter << ": " << this->getModificationSiteNum(*iter, 1) << " ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}

	void Fragment::updateCleavage( const CleavageCollection& cleavage_col )
	{
		// Use the cleavage information to update the fragment.
		CleavageCollection::const_iterator iter = cleavage_col.begin();
		for(; iter != cleavage_col.end(); iter++)
		{
				this->setFragmentation(iter->first, iter->second);
		}

	}

	bool Fragment::isCrossringCleavage() const
	{
		if(cleavage_sites.size() != 1)
			return false;

		if(cleavage_sites.begin()->first == "A" || cleavage_sites.begin()->first == "X")
			return true;
		else 
			return false;
	}

	bool Fragment::isGlycosidicBondCleavage() const
	{
		if(cleavage_sites.size() != 1)
			return false;

		if(cleavage_sites.begin()->first == "B" || cleavage_sites.begin()->first == "C" || cleavage_sites.begin()->first == "Y" || cleavage_sites.begin()->first == "Z")
			return true;
		else 
			return false;
	}

	std::string Fragment::getGeneralType() const
	{
		if(isCrossringCleavage()) return "C";

		if(isGlycosidicBondCleavage()) return "G";

		// Intact.
		if(cleavage_sites.size() == 0) return "I";

		std::multiset<std::string> general_type;
		for(CleavageCollection::const_iterator iter = cleavage_sites.begin(); iter != cleavage_sites.end(); iter++)
		{
			if(iter->first == "B" || iter->first == "C" || iter->first == "Y" || iter->first == "Z")
				general_type.insert("G");
			else if(iter->first == "A" || iter->first == "X")
				general_type.insert("C");
		}

		std::string final_type;
		for(std::multiset<std::string>::iterator iter = general_type.begin();
			iter != general_type.end(); iter++)
			final_type.append(*iter);

		return final_type;
	}

  int Fragment::getMonosaccharidePosition() const
  {
    if(!this->isTerminalCleavage())
      return -1;

    return (int)cleavage_sites.begin()->second.getMonosaccharideID();
  }

  size_t Fragment::getCleavageIndex() const
  {
    size_t id = this->getMonosaccharidePosition()+1;

    if(id == 0) return -1;
    std::string clv_type = this->getFragmentType();
    if(clv_type == "A" || clv_type == "B" || clv_type == "C") {
    } else 
      id = glyco_seq->getBranchByID(0).getUnitNum() - id;

    return id;

  }

  ModificationSites Fragment::getExpandedModificationSites( std::string mod_symbol, int flag /*= 1*/ )
  {
    size_t cleavage_index = this->getCleavageIndex();
    std::string clv_type = this->getFragmentType();
    ModificationSites mod_sites = glyco_seq->getModificationSitesBySymbol(mod_symbol, flag);
    
    if(clv_type == "A") {
      if(cleavage_index != glyco_seq->getBranchByID(0).getUnitNum()) {
        Fragment virtual_frag(glyco_seq);
        virtual_frag.setFragmentation("B", cleavage_index, 0,0);
        mod_sites = virtual_frag.getModificationSitesBySymbol(mod_symbol, flag);
      }

    } else if(clv_type == "X") {
      if(cleavage_index != glyco_seq->getBranchByID(0).getUnitNum()-1) {
        Fragment virtual_frag(glyco_seq);
        virtual_frag.setFragmentation("Y", cleavage_index+1, 0,0);
        mod_sites = virtual_frag.getModificationSitesBySymbol(mod_symbol, flag);
      }
        
    } else {
      mod_sites = this->getModificationSitesBySymbol(mod_symbol, flag);
    }

    return mod_sites;

  }

  std::string Fragment::getFragmentType() const
  {
    std::string cleavage_str;
		if(cleavage_sites.size() == 0)
			return "SEQ";

		for(CleavageCollection::const_iterator iter = cleavage_sites.begin(); iter!=cleavage_sites.end();iter++)
		{
			cleavage_str.append(iter->first);
			
		}
		return cleavage_str;
  }

  gag::ModificationSites Fragment::getReducedModificationSites( std::string mod_symbol, int flag /*= 1*/ )
  {
    size_t cleavage_index = this->getCleavageIndex();
    std::string clv_type = this->getFragmentType();
    ModificationSites mod_sites = ModificationSites();

    if(clv_type == "A") {
      if(cleavage_index > 1) {
        Fragment virtual_frag(glyco_seq);
        virtual_frag.setFragmentation("B", cleavage_index-1, 0,0);
        mod_sites = virtual_frag.getModificationSitesBySymbol(mod_symbol, flag);
      }

    } else if(clv_type == "X") {
      if(cleavage_index > 0) {
        Fragment virtual_frag(glyco_seq);
        virtual_frag.setFragmentation("Y", cleavage_index, 0,0);
        mod_sites = virtual_frag.getModificationSitesBySymbol(mod_symbol, flag);
      }

    } else {
      mod_sites = this->getModificationSitesBySymbol(mod_symbol, flag);
    }

    return mod_sites;
  }

  ModificationSites Fragment::getCompleteModificationSites( std::string mod_symbol, int flag /*= 1*/ )
  {
    return glyco_seq->getModificationSitesBySymbol(mod_symbol, flag);
  }

}
