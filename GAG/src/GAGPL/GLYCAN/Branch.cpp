/*
 * =====================================================================================
 *
 *       Filename:  Branch.cpp
 *
 *    Description:  Implementation file of class Branch.
 *
 *        Version:  1.0
 *        Created:  05/ 4/2012  9:07:08 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#include <iostream>
#include <GAGPL/GLYCAN/Branch.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>

namespace gag
{

	void Branch::addUnit(Monosaccharide& mono_unit)
	{	
		
		// Remove OH at link.end.
		if(mono_chain.size() != 0)
		{
			Linkage& lk = links.back();
			//Composition& cp = is.getFunctionalGroups().at(0).getComposition();
			FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
			//:) TBD The functional group should be specified by name.
			//FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol("OH");
			//mono_unit.removeFunctionalGroup(lk.end, fg);
			mono_unit.removeFunctionalGroup("OH", lk.end);
		}
		mono_chain.push_back(mono_unit);
		compo.add(mono_unit.getComposition());

	}

	void Branch::addLinkage(const Linkage& link)
	{
		// Modify the unit.
		//std::vector<Monosaccaride>::iterator iter = mono_chain.back();
		//Monosaccharide& mono = mono_chain.back();
		Monosaccharide& mono = mono_chain.at(link.nre_id);
		
		//InternalSite& is = mono.getInternalSites().at(link.start); 
		FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
		//:) TBD The functional group should be specified by name.
		FunctionalGroup& fg_h = fgt.getFunctionalGroupBySymbol("H");
		FunctionalGroup& fg_oh = fgt.getFunctionalGroupBySymbol("OH");
		FunctionalGroupChain chain;

		chain.push_back(std::make_pair(link.start, fg_oh));
		chain.push_back(std::make_pair(0, fg_h));

		mono.removeFunctionalGroupByChain(chain);
		// Retrieve the composition on specified position.
		//Composition& cp = is.getFunctionalGroups().at(0).site_gps.at(1);
		//FunctionalGroup fg("OH");
		//is.remove(fg, "H");

		links.push_back(link);
		compo.deduct("H");
	}

	void Branch::addExtension(Composition& cp)
	{
		re_extension.push_back(cp);
		compo.add(cp);
	}

	std::vector<Linkage> Branch::getNeighborLinks(const size_t mono_id)
	{
		// The structure is guaranteed to be linear.
		std::vector<Linkage> lk;
		if(links.size() == 1){
		}
		else if(mono_id == links.size()-1) {
			lk.push_back(links.at(mono_id-1));
		} else if(mono_id == 0) {
			lk.push_back(links.at(mono_id));
		} else {
			lk.push_back(links.at(mono_id - 1));
			lk.push_back(links.at(mono_id));
		}
		return lk;
	}

	Composition Branch::getSubComposition(const size_t start, const size_t end)
	{
		Composition sub_compo;

		for(std::vector<Monosaccharide>::iterator iter = mono_chain.begin() + start; 
				iter != mono_chain.begin()+end+1; iter++)
		{
			sub_compo.add(iter->getComposition());
		}

		return sub_compo;
	}
	
	// Recalculate the overall mass.
	void Branch::update()
	{
		compo.clear();

		for(std::vector<Monosaccharide>::iterator iter = mono_chain.begin(); 
				iter!=mono_chain.end(); iter++)
		{
			compo.add(iter->getComposition());
		}
		for(std::vector<Composition>::iterator iter = re_extension.begin(); 
			iter != re_extension.end(); iter++)
		{
			compo.add(*iter);
		}

	}
	//void Branch::addModification(const size_t mono_id, const size_t site_id, const std::string& plus)
	//{
	//	Monosaccharide& mono = mono_chain.at(mono_id);
	//	mono.addComposition(site_id, plus);
	//	compo.add(plus);
	//}
	//void Branch::addModification(const size_t mono_id, const size_t site_id, FunctionalGroup& ori, const Composition& plus, const Composition& minus)
	//{
	//	Monosaccharide& mono = mono_chain.at(mono_id);
	//	mono.addFunctionalGroup(site_id, ori, plus, minus);
	//	compo.add(plus);
	//	compo.deduct(minus);
	//}

	//void Branch::addModification(const size_t mono_id, const size_t site_id, const std::string& ori, const std::string& plus, const std::string& minus)
	//{
	//	Monosaccharide& mono = mono_chain.at(mono_id);
	//	mono.addFunctionalGroup(site_id, ori, plus, minus);
	//	compo.add(plus);
	//	compo.deduct(minus);
	//}

	//void Branch::replace(const size_t mono_id, const size_t site_id, FunctionalGroup& ori, FunctionalGroup& re)
	//{
	//	Monosaccharide& mono = mono_chain.at(mono_id);
	//	mono.replace(site_id, ori, re);
	//	compo.deduct(ori.getComposition());
	//	compo.add(re.getComposition());
	//}

	//void Branch::removeModification(const size_t mono_id, const size_t site_id, FunctionalGroup& ori, const Composition& minus)
	//{
	//	Monosaccharide& mono = mono_chain.at(mono_id);
	//	mono.removeFunctionalGroup(site_id, ori, minus);
	//	compo.deduct(minus);
	//}
	//void Branch::removeModification(const size_t mono_id, const size_t site_id, const std::string& ori, const std::string& minus)
	//{
	//	Monosaccharide& mono = mono_chain.at(mono_id);
	//	mono.removeFunctionalGroup(site_id, ori, minus);
	//	compo.deduct(minus);
	//}
	Composition Branch::getExtensionComposition()
	{
		Composition temp_compo;
		for(size_t i=0; i < re_extension.size(); i++)
			temp_compo.add(re_extension.at(i));

		return temp_compo;
	}

	void Branch::printStructure()
	{
		std::cout << "Branch: " << this->getBranchID()  << std::endl;
		std::cout << "Composition: " << this->getCompositionString()<< std::endl;
		std::vector<Monosaccharide>& mono_chain = this->getGlycanChainUnits();

		for(std::vector<Monosaccharide>::iterator iter = mono_chain.begin(); iter != mono_chain.end(); iter++)
		{
			iter->printStructure();
		}
	}

}





