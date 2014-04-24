/********************************************************************
	created:	2013/01/23
	created:	23:1:2013   20:04
	filename: 	ModifierTest.cpp
	file path:	GAG\test
	file base:	ModifierTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include <GAGPL/CHEMISTRY/Modifier.h>
#include <GAGPL/CHEMISTRY/FunctionalGroup.h>
#include <GAGPL/GLYCAN/Branch.h>
#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <iostream>

using namespace std;
using namespace gag;

void exploreInternalSite(Site& site)
{
	cout << "Site Composition: " << site.getSiteCompositionString() << endl;
	std::multimap<std::string, FunctionalGroup>& fg = site.fg_map;
	for(std::multimap<std::string, FunctionalGroup>::iterator iter = fg.begin(); iter != fg.end(); iter++)
	{
		cout << (*iter).second.getCompositionString() << endl;
	}
}

void exploreMonosaccharide(Monosaccharide& mono)
{
	cout << "Name: " << mono.getSymbol() << endl;
	cout << "Monosaccharide Composition: " << mono.getCompositionString() << endl;
	std::vector<Site>& sites = mono.getInternalSites();

	for(std::vector<Site>::iterator iter = sites.begin(); iter != sites.end(); iter++)
	{
		exploreInternalSite(*iter);
	}
}

void exploreBranch(Branch& bc)
{
	cout << "Branch: " << bc.getBranchID()  << endl;
	cout << "Composition: " << bc.getCompositionString()<< endl;
	std::vector<Monosaccharide>& mono_chain = bc.getGlycanChainUnits();

	for(std::vector<Monosaccharide>::iterator iter = mono_chain.begin(); iter != mono_chain.end(); iter++)
	{
		exploreMonosaccharide(*iter);
	}

}

int main( int argc, char *argv[]) 
{
	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();

	// Test 1. 2-AB labeling sequence.
	GlycanSequence gs;
	Branch bc(0);
	//Monosaccharide ms = fgt.getFunctionalGroupBySymbol("GlcA");
	bc.addUnit(fgt.getFunctionalGroupBySymbol("GlcA"));
	bc.addLinkage(Linkage(0,1,4,"alpha"));

	bc.addUnit(fgt.getFunctionalGroupBySymbol("GlcN"));
	bc.addLinkage(Linkage(1,1,4,"belta"));

	bc.addUnit(fgt.getFunctionalGroupBySymbol("GlcA"));
	bc.addLinkage(Linkage(2,1,4,"alpha"));

	bc.addUnit(fgt.getFunctionalGroupBySymbol("GlcNAc"));
	bc.update();

	gs.addBranch(bc);

	// Append modification.
	cout << "Test 2-AB labeling" << endl;
	Modifier modi;
	GlycanSequence gs1 = gs;
	FunctionalGroup& mono1 = gs.getBranchByID(0).getUnitByID(3);
	modi.modifyFunctionalGroup(mono1, "2-AB");
	gs1.update();
	std::cout << "Composition: " << gs1.getCompositionString() << std::endl;

	for(std::vector<Branch>::iterator iter = gs1.getBranches().begin(); iter != gs1.getBranches().end(); iter++)
	{
		exploreBranch(*iter);
	}
	// Test 2. Boons
	cout << "Test boons." << endl;
	GlycanSequence gs2 = gs;
	FunctionalGroup& mono2 = gs.getBranchByID(0).getUnitByID(3);
	modi.modifyFunctionalGroup(mono2, "PENTA", 1);
	gs2.update();
	std::cout << "Composition: " << gs2.getCompositionString() << std::endl;

	for(std::vector<Branch>::iterator iter = gs2.getBranches().begin(); iter != gs2.getBranches().end(); iter++)
	{
		exploreBranch(*iter);
	}

	// Test 3: containFunctionalGroup
	cout << "Contain test: " << endl;
	Monosaccharide mono = fgt.getFunctionalGroupBySymbol("GlcA");
	modi.modifyFunctionalGroup(mono, "S", 2);
	
	if(mono.containFunctionalGroup("HSO3", 2))
		cout << "contain HSO4" << endl;
	else 
		cout << "Test again." << endl;



}