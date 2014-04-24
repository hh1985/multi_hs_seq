/*
 * =====================================================================================
 *
 *       Filename:  Monosaccharide.cpp
 *
 *    Description:  Implementation file for class Monosaccharide.
 *
 *        Version:  1.0
 *        Created:  05/ 3/2012  7:47:07 AM
 *       Revision:  none
 *       Compiler: 	msvc 
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/GLYCAN/Monosaccharide.h>
//#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>

namespace gag
{
	//void Monosaccharide::add(size_t pos, const FunctionalGroup& ori, const Composition& plus, const Composition& minus)
	//{
	//	std::multiset<FunctionalGroup>& fg = internalsites.at(pos).getFunctionalGroups();

	//	std::multiset<FunctionalGroup>::iterator iter = fg.find(ori);
	//	
	//	if(iter != fg.end())
	//	{
	//		// Modify the composition of the functional group.
	//		(*iter).add(plus);
	//		(*iter).deduct(minus);

	//		// TBD: Modify the name. This is hard!!!
	//		// Update the composition of the site.
	//		internalsites.at(pos).add(plus);
	//		internalsites.at(pos).deduct(minus);

	//		compo.add(plus);
	//		compo.deduct(minus);
	//	}

	//}

	//void Monosaccharide::add(size_t pos, const std::string& ori, const std::string& plus, const std::string& minus)
	//{
	//	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	//	FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol(ori);
	//	Composition cp1(plus);
	//	Composition cp2(minus);
	//	(*this).add(pos, fg, cp1, cp2);
	//}
	//
	//void Monosaccharide::addComposition(size_t pos, const std::string& plus)
	//{
	//	internalsites.at(pos).add(plus);
	//	compo.add(plus);
	//}

	//void Monosaccharide::remove(size_t pos, FunctionalGroup& fg, const Composition& lost)
	//{
	//	std::multiset<FunctionalGroup>& fgs = internalsites.at(pos).getFunctionalGroups();

	//	std::multiset<FunctionalGroup>::iterator iter = fgs.find(fg);

	//	if(iter != fgs.end())
	//	{
	//		if(!(lost == (*iter).getComposition())){
	//			(*iter).deduct(lost);
	//			internalsites.at(pos).deduct(lost);
	//			compo.deduct(lost);
	//		}
	//		else {
	//			fgs.erase(iter);
	//			internalsites.at(pos).deduct(fg.getComposition());
	//			compo.deduct(fg.getComposition());
	//		}
	//	}
	//}

	//void Monosaccharide::remove(size_t pos, const std::string& ori, const std::string& lost)
	//{
	//	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	//	FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol(ori);
	//	Composition cp(lost);
	//	(*this).remove(pos, fg, cp);
	//}

	//void Monosaccharide::remove(size_t pos, FunctionalGroup& fg)
	//{
	//	std::multiset<FunctionalGroup>& fgs = internalsites.at(pos).getFunctionalGroups();

	//	std::multiset<FunctionalGroup>::iterator iter = fgs.find(fg);

	//	if(iter != fgs.end()) {
	//		fgs.erase(iter);
	//		internalsites.at(pos).deduct(fg.getComposition());
	//		compo.deduct(fg.getComposition());
	//	}
	//}

	//void Monosaccharide::remove(std::map<size_t, FunctionalGroup>& substract)
	//{
	//	std::map<size_t, FunctionalGroup>::iterator iter = substract.begin();
	//	for(;iter!=substract.end(); iter++)
	//	{
	//		(*this).remove(iter->first, iter->second);
	//	}
	//}

	//void Monosaccharide::replace(size_t pos, FunctionalGroup& ori, FunctionalGroup& re)
	//{
	//	internalsites.at(pos).replace(ori, re);
	//	compo.deduct(ori.getComposition());
	//	compo.add(re.getComposition());
	//}
	//InternalSite& Monosaccharide::getInternalSiteByRingID(const size_t r_id)
	//{
	//	if(r_id == 0)
	//		return internalsites.at(0);
	//	else if(r_id > 0 && (ring_start + r_id -1 <= ring_end))
	//		return internalsites.at(ring_start + r_id -1);
	//	else 
	//		throw std::runtime_error("Unqualified ring id.");
	//}
	//
	//Composition Monosaccharide::getCleavageShift()
	//{
	//	//Composition temp_compo;
	//	return internalsites.at(ring_start).getCleavageShift(0);
	//}

	//size_t Monosaccharide::getRingID(const size_t c_id)
	//{
	//	if(c_id == 0)
	//		return 0;
	//	else if(c_id <= ring_start)
	//		return 1;
	//	else if(c_id <= ring_end)
	//		return (c_id - ring_start + 1);
	//	else if(c_id > ring_end)
	//		return (ring_end - ring_start + 1);
	//	else 
	//		throw std::runtime_error("Unqualified carbon ID");
	//}

	//size_t Monosaccharide::getCarbonID(const size_t r_id)
	//{
	//	if(r_id == 0)
	//		return 0;
	//	else if(r_id >= 1 && (ring_start + r_id -1 <= ring_end))
	//		return (ring_start + r_id -1);
	//	else if(ring_start + r_id -1 > ring_end)
	//		throw std::runtime_error("Unqualified ring ID");

	//	return 9999;
	//}

	//Composition Monosaccharide::getSubCompositionByCarbonID(const size_t& s1, const size_t& s2)
	//{
	//	Composition chain_compo;

	//	if(s1 > s2)
	//		throw std::runtime_error("Unqualified carbon ID pair!");

	//	for(std::vector<InternalSite>::iterator iter = internalsites.begin() + s1; iter!= internalsites.begin() + s2 + 1; iter++)
	//	{
	//		chain_compo.add(iter->getComposition());
	//	}

	//	return chain_compo;
	//}
	//
	//Composition Monosaccharide::getSubCompositionByRingID(const size_t& s1, const size_t& s2)
	//{
	//	if(s1 > s2 || ring_start + s2 -1 > ring_end)
	//		throw std::runtime_error("Unqualified ring ID");

	//	Composition ring_compo;
	//	if(s1 == 0 && s2 == 0) // The whole ring. 
	//		ring_compo = (*this).getSubCompositionByCarbonID(0, internalsites.size()-1);
	//	else if(s1 == 0 && ring_start+s2-1 == ring_end) 
	//		ring_compo = (*this).getSubCompositionByCarbonID(1, internalsites.size()-1);
	//	else if(s1 == 0 && ring_start+s2-1 < ring_end) // Start from Carbon ID 1.
	//		ring_compo = (*this).getSubCompositionByCarbonID(1, ring_start + s2 -1);
	//	else if(ring_start+s2-1 < ring_end)
	//		ring_compo = (*this).getSubCompositionByCarbonID(ring_start+s1, ring_start+s2-1);
	//	else if(ring_start+s2-1 == ring_end)
	//		ring_compo = (*this).getSubCompositionByCarbonID(ring_start+s1, internalsites.size()-1);
	//	else 
	//		throw std::runtime_error("Strange things happened!");
	//	//ring_compo.updateString();
	//	return ring_compo;
	//}

	//bool Monosaccharide::hasFunctionalGroup( const size_t& pos, FunctionalGroup& fg )
	//{
	//	std::vector<InternalSite>::iterator iter = internalsites.begin()+pos;
	//	if(iter->hasFunctionalGroup(fg))
	//		return true;
	//	else 
	//		return false;
	//}

	//void Monosaccharide::addFunctionalGroupByChain( FunctionalGroupChain& chain )
	//{
	//	
	//}

	//void Monosaccharide::removeFunctionalGroupByChain( FunctionalGroupChain& chain )
	//{
	//	
	//}





}
