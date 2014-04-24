/*
 * =====================================================================================
 *
 *       Filename:  GlycanSequence.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/ 5/2012 10:15:09 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/GLYCAN/GlycanSequence.h>
//#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
//#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <boost/lexical_cast.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <iostream>
#include <fstream>

namespace gag
{
	void GlycanSequence::addBranch(Branch& bc)
	{
		glycan_sequence.push_back(bc);
		compo.add(bc.getComposition());
		// Initially, the originmap is set to be the mapping between the branch id and itself.
		originmap.insert(std::make_pair(bc.getBranchID(), bc.getBranchID()));
	}
	void GlycanSequence::addBranchLink(const size_t& re_id, const size_t& nre_id)
	{
		branch_links.insert(BranchMap::value_type(re_id, nre_id));

		// Update composition.
		Branch& bc = (*this).getBranchByID(nre_id);
		Linkage& lk = (*this).getBranchByID(re_id).getLinkages().back();
		bc.getUnitByID(0).removeFunctionalGroup("OH", lk.end);
		//bc.removeModification(0,lk.end,"OH","OH");
		compo.deduct("OH");
	}

	void GlycanSequence::addBranchLinks(const std::set<size_t>& re_ids, const size_t& nre_id)
	{
		BOOST_FOREACH(size_t re_id, re_ids)
		{
			this->addBranchLink(re_id, nre_id);
		}
		// update originmap.
		BranchOrigin::iterator iter1, iter2;
		iter1 = originmap.find(*re_ids.begin());
		iter2 = originmap.find(nre_id);
		iter2->second = iter1->second;


	}
	void GlycanSequence::updateChildrenIDs(const size_t branch_id)
	{
		BranchMap::right_const_iterator right_lower_iter = branch_links.right.lower_bound(branch_id);

		size_t min = branch_id;
		BranchDescendants::const_iterator d_iter = children.find(right_lower_iter->second);
		
		if(d_iter == children.end())
			min = right_lower_iter->second;
		else
			min = d_iter->second;

		children.insert(std::make_pair(branch_id, min));
	}

	void GlycanSequence::update()
	{
		compo.clear();
		
		for(std::vector<Branch>::iterator iter = glycan_sequence.begin(); 
				iter!=glycan_sequence.end(); iter++)
		{
			compo.add(iter->getComposition());
		}
		
	}

	std::pair<size_t, size_t> GlycanSequence::getDescendantBranchIDs(const size_t& branch_id)
	{	
		BranchDescendants::const_iterator d_iter = children.find(branch_id);
		if(d_iter != children.end())
			return std::pair<size_t, size_t>(d_iter->second, branch_id-1);
		else // Leaf branch.
			return std::pair<size_t, size_t>(branch_id, branch_id-1);
	}

	size_t GlycanSequence::getOriginBranchID(const size_t& branch_id)
	{
		BranchOrigin::iterator iter = originmap.find(branch_id);
		if(iter != originmap.end())
			return iter->second;
		else 
			return -1;
	}

	Composition GlycanSequence::getSubComposition(const size_t& id1, const size_t& id2)
	{
		Composition sub_compo;
		for(size_t i = id1; i < id2+1; i++)
		{
			sub_compo.add((*this).getBranchByID(i).getComposition());
		}
		return sub_compo;
	}

	Composition GlycanSequence::getSubTreeComposition(const size_t& branch_id)
	{
		std::pair<size_t, size_t>& range = (*this).getDescendantBranchIDs(branch_id);

			return (*this).getSubComposition(range.first, range.second);
	}
	
	Composition GlycanSequence::getTreeComposition(const size_t& branch_id)
	{
		std::pair<size_t, size_t>& range = (*this).getDescendantBranchIDs(branch_id);

			return (*this).getSubComposition(range.first, range.second+1);
	}

	SequenceCode GlycanSequence::getGAGSequenceCode()
	{
		Branch& bc = glycan_sequence.at(0);

		SequenceCode seq_code;

		for(size_t i = 0; i < bc.getUnitNum(); i++)
		{
			std::string str;
			size_t num = 0;
			std::string symbol = bc.getUnitByID(i).getSymbol();

			FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
			//FunctionalGroup& fg_ac = fgt.getFunctionalGroupBySymbol("Ac");
			//FunctionalGroup& fg_sulf = fgt.getFunctionalGroupBySymbol("HSO4");
			std::string fg_ac = "Ac"; std::string fg_sulf = "HSO3";

			// Check modification.
			if(symbol == "DeltaGlcA"){
				str = "D";
			} else if(symbol == "GlcA"){
				str = "G";
			} else if(symbol == "GlcN"){
				// Check acetation situation.
				if(bc.getUnitByID(i).containFunctionalGroup(fg_ac, 2))
					str = "A";
				else if(bc.getUnitByID(i).containFunctionalGroup(fg_sulf, 2))
					str = "S";
				else
					str = "H";
			} else if(symbol == "GlcNAc"){
				str = "A";
			}
			if(str == "D" || str == "G") {
				for(size_t j = 0; j < bc.getUnitByID(i).getInternalSites().size(); j++)
				{
					if(bc.getUnitByID(i).containFunctionalGroup(fg_sulf, j))
						num = j; 
				}
			} else if(str == "H" || str == "S" || str == "A") {
				for(size_t j = 3; j < bc.getUnitByID(i).getInternalSites().size(); j++)
				{
					if(bc.getUnitByID(i).containFunctionalGroup(fg_sulf, j))
						num += j; 
				}
			}
			//cout << "Sequence: " << str << num;
			seq_code.push_back(std::make_pair(str, num));

		}
		return seq_code;
		//cout << endl;
	}

	std::string GlycanSequence::getGAGSequenceCodeString()
	{
		SequenceCode seq_code = this->getGAGSequenceCode();
		std::string seq_str;
		SequenceCode::iterator iter = seq_code.begin();
		for(; iter!=seq_code.end(); iter++)
		{
			seq_str.append(iter->first);
			seq_str.append(boost::lexical_cast<std::string>(iter->second));
		}
		return seq_str;
	}

	void GlycanSequence::initializeModificationSites( const std::string& mod_symbol )
	{
		std::vector<Branch>& branches = this->getBranches();
		for(size_t bc_id=0; bc_id < branches.size(); bc_id++ )
		{
			std::vector<Monosaccharide>& mono_units = branches[bc_id].getGlycanChainUnits();
			for(size_t mono_id = 0; mono_id < mono_units.size(); mono_id++)
			{
				// Notice: just omit the atom information.
				std::vector<size_t> site_ids = mod_assistant.getModificationSetByPosition(mono_units[mono_id], mod_symbol);
				BOOST_FOREACH(size_t& site_id, site_ids)
				{
					mod_assistant.addModification(ModificationPosition(bc_id, mono_id, site_id), mod_symbol, true);
				}
			}
		}
	}

	void GlycanSequence::modifyModificationStatus( const std::string& mod_symbol, const ModificationPosition& pos, int status)
	{
		// Set the status of the rest of modification type on the same position.
		//int new_status;
		//if(status == 0)
		//	new_status = -1;
		//else if(status == 1)
		//	new_status = 1;
		//else if(status == -1)
		//	new_status = -1;

		//// 1. Locate the position of the modification.
		//ModificationByPosition& mod_by_pos = mod_pos_map.getContainer().get<mod_position>();
		//std::pair<ModificationByPosition::iterator, ModificationByPosition::iterator> p = mod_by_pos.equal_range(pos);
		//
		//if(p.first == p.second) {
		//	std::cout << "The specified position is not found for the modification type" << std::endl;
		//}

		//for(ModificationByPosition::iterator iter = p.first; iter!=p.second; iter++)
		//{
		//	// 2.a Update the information.
		//	if(iter->mod_symbol == mod_symbol) {
		//		mod_pos_map.modifyModificationStatus(iter, status); // Do not change the modification type.
		//	} else {
		//		mod_pos_map.modifyModificationStatus(iter, new_status);
		//	}
		//}
		mod_assistant.modifyModificationStatus(mod_symbol, pos, status);
	}

	void GlycanSequence::updateModification( const std::string& mod_symbol, ModificationPosition& pos )
	{
		Monosaccharide& mono = this->getBranchByID(pos.getBranchID()).getUnitByID(pos.getMonosaccharideID());
		mod_assistant.modifyFunctionalGroup(mono, mod_symbol);
		// Update the modification status to unavailable.
		this->modifyModificationStatus(mod_symbol, pos);
		this->getBranchByID(pos.getBranchID()).update();
		this->update();
	}

	void GlycanSequence::updateModification( const std::string& mod_symbol, ModificationPosition& pos, size_t site_id )
	{
		Monosaccharide& mono = this->getBranchByID(pos.getBranchID()).getUnitByID(pos.getMonosaccharideID());
		mod_assistant.modifyFunctionalGroup(mono, mod_symbol, site_id);
		// Update site id.
		pos.site_id = site_id;
		// Update the modification status to unavailable.
		this->modifyModificationStatus(mod_symbol, pos);
		this->getBranchByID(pos.getBranchID()).update();
		this->update();
	}

	void GlycanSequence::buildByGAGComposition( const GlycanComposition& glycan_compo, std::string init_mono)
	{
		Branch bc(0);

		FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
		std::map<std::string, int>::const_iterator d_iter = glycan_compo.mono_compo.find("DeltaGlcA");
		std::map<std::string, int>::const_iterator glca_iter = glycan_compo.mono_compo.find("GlcA");
		std::map<std::string, int>::const_iterator glcn_iter = glycan_compo.mono_compo.find("GlcN");
		std::map<std::string, int>::const_iterator s_iter = glycan_compo.mod_compo.find("SO3");
		std::map<std::string, int>::const_iterator ac_iter = glycan_compo.mod_compo.find("Ac");

		int dglca_num = d_iter == glycan_compo.mono_compo.end() ? 0 : d_iter->second;
		int glca_num = glca_iter == glycan_compo.mono_compo.end() ? 0 : glca_iter->second;
		int glcn_num = glcn_iter == glycan_compo.mono_compo.end() ? 0 : glcn_iter->second;
		int s_num = s_iter == glycan_compo.mod_compo.end() ? 0 : s_iter->second;
		int ac_num = ac_iter == glycan_compo.mod_compo.end() ? 0 : ac_iter->second;

		std::string init_error = "Please check the configuration of initial monosaccharide unit.";
		std::string num_error = "The number of monosaccharides for each type does not fit the model.";
		std::string sulfur_error = "The number of modifications doesn't fit";

		/* Check if the composition meet GAG requirement. */
		if(dglca_num > 1 || abs(glca_num - glcn_num) > 1) {
			std::cout << num_error << std::endl;
			return;
		}
		if(dglca_num == 1 && init_mono != "GlcA") {
			std::cout << init_error << std::endl;
			return;
		}
		if(ac_num > glcn_num || (ac_num + s_num) > (dglca_num + glca_num + 3 * glcn_num) )
		{
			std::cout << sulfur_error << std::endl;
			return;
		}

		// Decide the sequence.
		int cur_index = 0;

		if(dglca_num == 1) {
			bc.addUnit(fgt.getFunctionalGroupBySymbol("DeltaGlcA"));
			cur_index++;
		} 

		if(glca_num == 0 && glcn_num == 0) {
			this->addBranch(bc);
			return;
		}

		std::vector<std::string> units;
		units.push_back("GlcN"); units.push_back("GlcA");

		std::vector<std::string> links;
		links.push_back("beta"); links.push_back("alpha");

		if(init_mono == "GlcA") {
			if(cur_index == 1){
				bc.addLinkage(Linkage(cur_index-1, 1, 4, "alpha"));
			} else {
				std::swap(units[0], units[1]);
				std::swap(links[0], links[1]);
			}
		} else if(init_mono == "GlcN") {

		}

		while(1) {
			for(int i = 0; i < 2; i++) {
				bc.addUnit(fgt.getFunctionalGroupBySymbol(units[i]));
				cur_index++;
				if(cur_index != dglca_num + glca_num + glcn_num) {
					bc.addLinkage(Linkage(cur_index-1, 1, 4, links[i]));
				} else {
					break;
				}
			}
			if(cur_index == dglca_num + glca_num + glcn_num)
				break;
		}		

		this->addBranch(bc);

		// add modification information.
		for(std::map<std::string, int>::const_iterator iter = glycan_compo.mod_compo.begin(); iter != glycan_compo.mod_compo.end(); iter++)
			this->initializeModificationSites(iter->first);

		mod_constraint = glycan_compo.mod_compo;

	}

	std::set<std::string> GlycanSequence::getModificationTypes()
	{
		std::set<std::string> mod_types;
		for(std::map<std::string, int>::iterator iter = mod_constraint.begin();
			iter != mod_constraint.end(); iter++)
		{
			mod_types.insert(iter->first);
		}
		return mod_types;
	}

	int GlycanSequence::getModificationConstraint( const std::string& mod_symbol )
	{
		std::map<std::string, int>::iterator iter = mod_constraint.find(mod_symbol);
		if(iter != mod_constraint.end())
			return iter->second;
		else
			return 0;
	}

	ModificationSites GlycanSequence::getModificationSitesBySymbol(const std::string& mod_symbol)
	{
		std::set<ModificationPosition> mod_set;
		// Partial key extract.
		std::pair<ModificationByCompositeKey::iterator, ModificationByCompositeKey::iterator> p = mod_assistant.getModificationContainer().get<mod_compo_key>().equal_range(boost::make_tuple(mod_symbol));

		while(p.first != p.second) {
			mod_set.insert(p.first->mod_pos);
			++p.first;
		}

		return mod_set;
	}

	ModificationSites GlycanSequence::getModificationSitesBySymbol(const std::string& mod_symbol, int status)
	{
		std::set<ModificationPosition> mod_set;
		// Partial key extract.
		std::pair<ModificationByCompositeKey::iterator, ModificationByCompositeKey::iterator> p = mod_assistant.getModificationContainer().get<mod_compo_key>().equal_range(boost::make_tuple(mod_symbol, status));

		while(p.first != p.second) {
			mod_set.insert(p.first->mod_pos);
			++p.first;
		}

		return mod_set;
	}

	void GlycanSequence::printStructure()
	{
		// Iterate over all branches and explore the structure.
		for(std::vector<Branch>::iterator iter = this->getBranches().begin(); iter != this->getBranches().end(); iter++)
		{
			iter->printStructure();
			std::cout << "Children: " << this->getDescendantBranchIDs(iter->getBranchID()).first << "-" << this->getDescendantBranchIDs(iter->getBranchID()).second << std::endl;
		}

		// Print modification information.
		ModificationByPosition& mod_pos_index = this->getModificationSetByPosition();
		ModificationByPosition::iterator mod_iter = mod_pos_index.begin();
		for(; mod_iter != mod_pos_index.end(); mod_iter++)
		{
			std::cout << "Position: ";
			std::cout << mod_iter->mod_pos.getBranchID() << " " << mod_iter->mod_pos.getMonosaccharideID() << " " << mod_iter->mod_pos.site_id << " " << mod_iter->mod_pos.atom << std::endl;

			std::string status;
			if(mod_iter->mod_status == 1)
				status = "Available";
			else if(mod_iter->mod_status == 0)
				status = "Occupied";
			else if(mod_iter->mod_status == -1)
				status = "Unavailable";

			std::cout << "Type: " << mod_iter->mod_symbol << std::endl;
			std::cout << "Status: " << status << std::endl;

		}

	}

	void GlycanSequence::addModification( const std::string& mod_symbol, ModificationPosition& pos )
	{
		Monosaccharide& mono = this->getBranchByID(pos.getBranchID()).getUnitByID(pos.getMonosaccharideID());
		mod_assistant.modifyFunctionalGroup(mono, mod_symbol, pos.site_id);
		this->getBranchByID(pos.getBranchID()).update();
		this->update();
	}

	void GlycanSequence::addModificationConstraint( const std::string& mod_symbol, int number )
	{

		// Insert/Update the records.
		std::map<std::string, int>::iterator mod_iter = mod_constraint.find(mod_symbol);
		if(mod_iter == mod_constraint.end()) // Not found.
			mod_constraint.insert(std::make_pair(mod_symbol, number));
		else // Found: add corresponding composition number.
			mod_iter->second += number;

		// Initialize modification site.
		this->initializeModificationSites(mod_symbol);
	}

  ModificationSites GlycanSequence::getComplementaryModificationSites( const ModificationSites& ms, std::string mod_symbol )
  {
    ModificationSites complete_sites = this->getModificationSitesBySymbol(mod_symbol, 1);
    ModificationSites diff_sites;

    std::set_difference(complete_sites.begin(), complete_sites.end(), ms.begin(), ms.end(), std::inserter(diff_sites, diff_sites.end()));

    return diff_sites;	
  }



}
