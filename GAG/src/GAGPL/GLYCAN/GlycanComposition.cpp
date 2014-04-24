/********************************************************************
	created:	2013/04/04
	created:	4:4:2013   9:31
	filename: 	GlycanComposition.cpp
	file path:	GAGPL/GLYCAN
	file base:	GlycanComposition
	file ext:	cpp
	author:		Han Hu (HH), hh.earlydays@gmail.com
	
	purpose:	
*********************************************************************/

#include <GAGPL/GLYCAN/GlycanComposition.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <iostream>

namespace gag
{
	//GlycanComposition::GlycanComposition(const std::map<std::string, int>& compo)
	//{
	//	std::map<std::string, int>::const_iterator map_iter = compo.begin();
	//	for(; map_iter != compo.end(); map_iter++)
	//	{
	//		// Check if the functional group is defined.
	//		std::map<std::string, int>::iterator glycan_iter = mono_compo.find(map_iter->first);
	//		if(glycan_iter == mono_compo.end()) // Not found.
	//			mono_compo.insert(std::make_pair(map_iter->first, map_iter->second));
	//		else // Found: update corresponding composition number.
	//			glycan_iter->second = map_iter->second;
	//	}
	//}

	void GlycanComposition::addGlycanComposition( const std::string& mono, int number )
	{
		// Check if the functional group is defined.
		
		if(!fgt.containFunctionalGroup(mono)) {
			// Skip the mono.
			std::cout << "Functional group " << mono << " was not found!" << std::endl;
			return;
		}

		// Insert/Update the records.
		std::map<std::string, int>::iterator glycan_iter = mono_compo.find(mono);
		if(glycan_iter == mono_compo.end()) // Not found.
			mono_compo.insert(std::make_pair(mono, number));
		else // Found: update corresponding composition number.
			glycan_iter->second = number;
	}

	//GlycanSequence GlycanComposition::getHSBackboneSequence(const std::string init_mono)
	//{
	//	GlycanSequence gs;
	//	Branch bc(0);

	//	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	//	std::map<std::string, int>::iterator d_iter = glycan_compo.find("DeltaGlcA");
	//	std::map<std::string, int>::iterator glca_iter = glycan_compo.find("GlcA");
	//	std::map<std::string, int>::iterator glcn_iter = glycan_compo.find("GlcN");
	//	std::map<std::string, int>::iterator s_iter = glycan_compo.find("HSO3");
	//	std::map<std::string, int>::iterator ac_iter = glycan_compo.find("Ac");

	//	int dglca_num = d_iter == glycan_compo.end() ? 0 : d_iter->second;
	//	int glca_num = glca_iter == glycan_compo.end() ? 0 : glca_iter->second;
	//	int glcn_num = glcn_iter == glycan_compo.end() ? 0 : glcn_iter->second;
	//	int s_num = s_iter == glycan_compo.end() ? 0 : s_iter->second;
	//	int ac_num = ac_iter == glycan_compo.end() ? 0 : ac_iter->second;

	//	std::string init_error = "Please check the configuration of initial monosaccharide unit.";
	//	std::string num_error = "The number of monosaccharides for each type does not fit the model.";
	//	std::string sulfur_error = "The number of modifications doesn't fit";

	//	/* Check if the composition meet GAG requirement. */
	//	if(dglca_num > 1 || abs(glca_num - glcn_num) > 1) {
	//		std::cout << num_error << std::endl;
	//		return gs;
	//	}
	//	if(dglca_num == 1 && init_mono != "GlcA") {
	//		std::cout << init_error << std::endl;
	//		return gs;
	//	}
	//	if(ac_num > glcn_num || (ac_num + s_num) > (dglca_num + glca_num + 3 * glcn_num) )
	//	{
	//		std::cout << sulfur_error << std::endl;
	//		return gs;
	//	}

	//	// Decide the sequence.
	//	int cur_index = 0;
	//	
	//	if(dglca_num == 1) {
	//		bc.addUnit(fgt.getFunctionalGroupBySymbol("DeltaGlcA"));
	//		cur_index++;
	//	} 

	//	if(glca_num == 0 && glcn_num == 0) {
	//		gs.addBranch(bc);
	//		return gs;
	//	}

	//	std::vector<std::string> units;
	//	units.push_back("GlcN"); units.push_back("GlcA");
	//	
	//	std::vector<std::string> links;
	//	links.push_back("beta"); links.push_back("alpha");
	//	
	//	if(init_mono == "GlcA") {
	//		if(cur_index == 1){
	//			bc.addLinkage(Linkage(cur_index-1, 1, 4, "alpha"));
	//		} else {
	//			std::swap(units[0], units[1]);
	//			std::swap(links[0], links[1]);
	//		}
	//	} else if(init_mono == "GlcN") {
	//	
	//	}

	//	while(1) {
	//		for(int i = 0; i < 2; i++) {
	//			bc.addUnit(fgt.getFunctionalGroupBySymbol(units[i]));
	//			cur_index++;
	//			if(cur_index != dglca_num + glca_num + glcn_num) {
	//				bc.addLinkage(Linkage(cur_index-1, 1, 4, links[i]));
	//			} else {
	//				break;
	//			}
	//		}
	//		if(cur_index == dglca_num + glca_num + glcn_num)
	//			break;
	//	}		

	//	gs.addBranch(bc);

	//	// add modification information.
	//	for(std::map<std::string, int>::iterator iter = mod_compo.begin(); 
	//		iter != mod_compo.end(); iter++)
	//		gs.initializeModificationSites(iter->first);

	//	gs.mod_constraint = mod_compo;
	//	
	//	return gs;

	//}

	void GlycanComposition::addModificationComposition( const std::string& mod_symbol, int number )
	{

		// Insert/Update the records.
		std::map<std::string, int>::iterator mod_iter = mod_compo.find(mod_symbol);
		if(mod_iter == mod_compo.end()) // Not found.
			mod_compo.insert(std::make_pair(mod_symbol, number));
		else // Found: add corresponding composition number.
			mod_iter->second += number;
	}


}