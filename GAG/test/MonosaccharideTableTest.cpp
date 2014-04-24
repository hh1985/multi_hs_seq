/*
 * =====================================================================================
 *
 *       Filename:  MonosaccharideUnitTableTest.cpp
 *
 *    Description:  Test file.
 *
 *        Version:  1.0
 *        Created:  05/ 3/2012 10:51:07 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#include <GAGPL/GLYCAN/Monosaccharide.h>
#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/GLYCAN/InternalSite.h>
#include <iostream>

//#include <GAGPL/CHEMISTRY/PeriodicTable.h>
using namespace std;
using namespace gag;

void exploreInternalSite(InternalSite& site)
{
	cout << "Composition: " << site.getCompositionString() << endl;
	std::multiset<FunctionalGroup>& fg = site.getFunctionalGroups();
	for(std::multiset<FunctionalGroup>::iterator iter = fg.begin(); iter != fg.end(); iter++)
	{
		cout << (*iter).getCompositionString() << endl;
	}
}

int main ( int argc, char *argv[] )
{


	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	fgt.load();

	MonosaccharideUnitTable& mut = MonosaccharideUnitTable::Instance();
	mut.load();

	Monosaccharide& mono = mut.getMonosaccharideBySymbol("NeuAc");

	cout << "Composition: " << mono.getCompositionString() << endl;

	cout << "Name: " << mono.getName() << endl;

	cout << "Ring: " << mono.getRingPosByRingID().first << " " << mono.getRingPosByRingID().second << endl;

	std::vector<InternalSite>& sites = mono.getInternalSites();

	for(std::vector<InternalSite>::iterator iter = sites.begin(); iter != sites.end(); iter++)
	{
		cout << (*iter).getCompositionString() << endl;
	}
	
	cout << "Function: getInternalSiteByCarbonID 3" << endl;
	InternalSite& site = mono.getInternalSiteByCarbonID(3);
	exploreInternalSite(site);

	cout << "Function: getInternalSiteByRingID 5" << endl;
	InternalSite& site1 = mono.getInternalSiteByRingID(5);
	exploreInternalSite(site1);

	cout << "Function: getSubCompositionByCarbonID 1-5" << endl;
	Composition cp1 = mono.getSubCompositionByCarbonID(1,5);
	cout << "CP1: " << cp1.getCompositionString() << endl;
	
	cout << "Function: getSubCompositionByRingID 1-5" << endl;
	Composition cp2 = mono.getSubCompositionByRingID(1,5);
	cout << "CP2: " << cp2.getCompositionString() << endl;
	Composition cp3 = mono.getSubCompositionByRingID(2,4);
	cout << "CP3: " << cp3.getCompositionString() << endl;

	cout << "Function: add() " << endl;
	size_t test_pos = 4;
	cout << "Before: " << endl;
	exploreInternalSite(mono.getInternalSiteByCarbonID(test_pos));
	cout << "After " << endl;
	
	FunctionalGroup& fg1 = fgt.getFunctionalGroupBySymbol("OH");
	Composition cp11("HSO4");
	Composition cp21("H2O");
	mono.add(test_pos,fg1,cp11,cp21);
	InternalSite& site2 = mono.getInternalSiteByCarbonID(test_pos);
	exploreInternalSite(site2);

	cout << "Function: remove() " << endl;
	cout << "Before: " << endl;
	exploreInternalSite(mono.getInternalSiteByCarbonID(test_pos));
	cout << "After " << endl;
	Composition cp31("SO3");
	mono.remove(test_pos, fg1, cp31);
	exploreInternalSite(site2);

	cout << "Cleavage shift: ";
	Composition temp_compo = mono.getCleavageShift();
	cout << temp_compo.getCompositionString() << endl;
	FunctionalGroup& fg_nh = fgt.getFunctionalGroupBySymbol("2-AB");
	FunctionalGroup& fg_oh = fgt.getFunctionalGroupBySymbol("OH");
	mono.replace(2,fg_oh,fg_nh);
	
	cout << "After replace: ";
	cout << mono.getCleavageShift().getCompositionString() << endl;
	
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
