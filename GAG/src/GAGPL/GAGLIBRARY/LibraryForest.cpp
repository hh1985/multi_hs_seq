#include <GAGPL/GAGLIBRARY/LibraryForest.h>
#include <boost/foreach.hpp>
#include <GAGPL/MATH/combination.h>
#include <GAGPL/CHEMISTRY/Modifier.h>

namespace gag
{

	ModificationSiteList LibraryForest::getCandidateAcetationSites()
	{
		ModificationSiteList sitelist;
		size_t branch_num = gag_seq.getBranchSize();
		for(size_t bc_id = 0; bc_id < branch_num; bc_id++)
		{
			size_t mono_num = gag_seq.getBranchByID(bc_id).getUnitNum();
			for(size_t mono_id = 0; mono_id < mono_num; mono_id++)
			{
				// Identify the position of the acetation.
				if(gag_seq.getBranchByID(bc_id).getUnitByID(mono_id).getSymbol() == "GlcN")
				{
					// TBD. More intelligent way of getting the acetate number.
					ModificationPosition site(bc_id, mono_id, 2, "N");
					sitelist.push_back(site);
				}
			}
		}
		return sitelist;
	}

	ModificationSiteList LibraryForest::getCandidateSulfationSites()
	{
		ModificationSiteList sitelist;
		size_t branch_num = gag_seq.getBranchSize();
		for(size_t bc_id = 0; bc_id < branch_num; bc_id++)
		{
			size_t mono_num = gag_seq.getBranchByID(bc_id).getUnitNum();
			for(size_t mono_id = 0; mono_id < mono_num; mono_id++)
			{
				// Identify the position of the acetation.
				std::string symbol = gag_seq.getBranchByID(bc_id).getUnitByID(mono_id).getSymbol();
				if(symbol == "GlcN") {
					// TBD. More intelligent way of getting the acetate number.
					ModificationPosition site1(bc_id, mono_id, 2, "N");
					ModificationPosition site2(bc_id, mono_id, 3, "O");
					ModificationPosition site3(bc_id, mono_id, 6, "O");
					sitelist.push_back(site1);
					sitelist.push_back(site2);
					sitelist.push_back(site3);
				} else if(symbol == "DeltaGlcA" || symbol == "GlcA") {
					ModificationPosition site(bc_id, mono_id, 2, "O");
					sitelist.push_back(site);
				}
			}
		}
		return sitelist;
	}

	void LibraryForest::generateTreeByComposition(const size_t& num_s, const size_t& num_ac)
	{
		// Based on the gag sequence backbone, guess all the combinations of the sulfation sites.
		ModificationSiteList mod_ac = this->getCandidateAcetationSites();
		ModificationSiteList mod_sulf = this->getCandidateSulfationSites();
		if(mod_ac.size() < num_ac || mod_sulf.size() < num_s)
			throw std::runtime_error("Illegal number of sulfation or acetation.");

		if(mod_ac.size() == 0) {
			ModificationPatternList pattern_sulf = this->guessModificationPositions(mod_sulf, num_s);
			ModificationSiteList selection_sulf, selection_ac;
			BOOST_FOREACH(selection_sulf, pattern_sulf)
			{
				//this->applyModification(selection_sulf, selection_ac);
				mod_pattern_list.push_back(std::make_pair(selection_sulf, selection_ac));
				
			}
		} else {
			ModificationPatternList pattern_ac = this->guessModificationPositions(mod_ac, num_ac);
			ModificationSiteList selection_ac;
			BOOST_FOREACH(selection_ac, pattern_ac)
			{
				ModificationPosition site_ac;
				ModificationPosition site_sulf;
				ModificationSiteList temp_mod_sulf = mod_sulf;
				// remove the overlapping site from sulfation sites.
				BOOST_FOREACH(site_ac, selection_ac)
				{
					ModificationSiteList::iterator iter_sulf = temp_mod_sulf.begin();
					//for(; iter_sulf != temp_mod_sulf.end(); iter_sulf++)
					while(iter_sulf != temp_mod_sulf.end())
					{
						if(site_ac == *iter_sulf) {
							iter_sulf = temp_mod_sulf.erase(iter_sulf);
							break;
						}
						else
							++iter_sulf;
					}
	
				}

				ModificationPatternList pattern_sulf = this->guessModificationPositions(temp_mod_sulf, num_s);
				ModificationSiteList selection_sulf;
				BOOST_FOREACH(selection_sulf, pattern_sulf)
				{
					//this->applyModification(selection_sulf, selection_ac);
					mod_pattern_list.push_back(std::make_pair(selection_sulf, selection_ac));
				}
				
			}
		}
		
	}


	// Consider about overloading the function in the future after defining the modification class.
	LibraryTree LibraryForest::applyModification(GlobalModPattern& global_pattern)
	{
		ModificationSiteList& sulf_sites = global_pattern.first;
		ModificationSiteList& ac_sites = global_pattern.second;
		// Copy constructor.
		GlycanSequence gs(gag_seq);

		ModificationSiteList::iterator iter = sulf_sites.begin();
		//std::string mod_s_compo = "H2SO4";
		Modifier modifier;
		for(; iter != sulf_sites.end(); iter++)
		{
			//if(iter->atom == "N") {
			//	gs.getBranchByID(iter->branch_id).addModification(iter->mono_id, iter->site_id, "NH2", mod_s_compo, "H2O");		
			//	gs.getBranchByID(iter->branch_id).getUnitByID(iter->mono_id).getInternalSiteByCarbonID(iter->site_id).addFunctionalGroup("HSO4");
			//} else if(iter->atom == "O") {
			//	gs.getBranchByID(iter->branch_id).addModification(iter->mono_id, iter->site_id, "OH", mod_s_compo, "H2O");
			//	gs.getBranchByID(iter->branch_id).getUnitByID(iter->mono_id).getInternalSiteByCarbonID(iter->site_id).addFunctionalGroup("HSO4");
			//}
			modifier.modifyGlycanSequence(gs, "S", *iter, iter->site_id);
			// Add the modification information into the InternalSite records. This should not modify the composition.
			
		}

		ModificationSiteList::iterator iter1 = ac_sites.begin();
		std::string mod_a_compo = "C2H4O2";
		for(; iter1 != ac_sites.end(); iter1++)
		{
			//if(iter1->atom == "N") {
			//	gs.getBranchByID(iter1->branch_id).addModification(iter1->mono_id, iter1->site_id, "NH2", mod_a_compo, "H2O");
			//	gs.getBranchByID(iter1->branch_id).getUnitByID(iter1->mono_id).getInternalSiteByCarbonID(iter1->site_id).addFunctionalGroup("Ac");
			//}
			modifier.modifyGlycanSequence(gs, "Ac", *iter1, iter1->site_id);
		}
		gs.update();

		size_t max_cleavage_num = 2;
		LibraryTree tree(gs, max_cleavage_num);
		
		return tree;
	}

	ModificationPatternList LibraryForest::guessModificationPositions( ModificationSiteList& sites, size_t num )
	{
		ModificationPatternList mod_patterns;
		// Initialization.
		ModificationSiteList mod_sites(sites.begin(), sites.begin() + num);
		
		// Generate all the combinations of the pattern.
		do {
			mod_patterns.push_back(mod_sites);
		} 
		while(stdcomb::next_combination(sites.begin(), sites.end(), mod_sites.begin(), mod_sites.end()));

		return mod_patterns;
	}

}