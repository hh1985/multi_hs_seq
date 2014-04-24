#include <GAGPL/GLYCAN/GlycanSequence.h>
//#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/GLYCAN/GlycanComposition.h>
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

int main(int argc, char* argv[])
{
	// Load Functional Group Table.
	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();

	GlycanSequence gs;

	// Build up a Branch of only two units.
	Branch bc1(0);
	Monosaccharide ms1 = fgt.getFunctionalGroupBySymbol("Man");
	// Unsaturated unit. This part should be realized by a class called
	// Reaction in the future.
	bc1.addUnit(ms1);	
	Linkage lk1(0,1,6,"alpha");
	bc1.addLinkage(lk1);
	gs.addBranch(bc1);

	Branch bc2(1);
	Monosaccharide ms2 = fgt.getFunctionalGroupBySymbol("Man");
	bc2.addUnit(ms2);
	Linkage lk2(0,1,3,"alpha");
	bc2.addLinkage(lk2);
	gs.addBranch(bc2);

	Branch bc3(2);
	Monosaccharide ms3 = fgt.getFunctionalGroupBySymbol("Man");
	bc3.addUnit(ms3);
	Linkage lk3(0,1,6,"alpha");
	bc3.addLinkage(lk3);
	gs.addBranch(bc3);
	gs.addBranchLink(0,2);
	gs.addBranchLink(1,2);
	gs.updateChildrenIDs(2);

	Branch bc4(3);
	Monosaccharide ms4 = fgt.getFunctionalGroupBySymbol("GlcNAc");
	bc4.addUnit(ms4);
	Linkage lk4(0,1,4,"belta");
	bc4.addLinkage(lk4);
	gs.addBranch(bc4);

	Branch bc5(4);
	Monosaccharide ms5 = fgt.getFunctionalGroupBySymbol("Man");
	bc5.addUnit(ms5);
	Linkage lk5(0,1,3,"alpha");
	bc5.addLinkage(lk5);
	gs.addBranch(bc5);

	Branch bc6(5);
	Monosaccharide ms6 = fgt.getFunctionalGroupBySymbol("Man");
	bc6.addUnit(ms6);
	Linkage lk6(0,1,4,"belta");
	bc6.addLinkage(lk6);
	Monosaccharide ms7 = fgt.getFunctionalGroupBySymbol("GlcNAc");
	bc6.addUnit(ms7);
	Linkage lk7(1,1,4,"belta");
	bc6.addLinkage(lk7);
	Monosaccharide ms8 = fgt.getFunctionalGroupBySymbol("GlcNAc");
	bc6.addUnit(ms8);
	gs.addBranch(bc6);
	gs.addBranchLink(2,5);
	gs.addBranchLink(3,5);
	gs.addBranchLink(4,5);
	gs.updateChildrenIDs(5);

	// Print results.
	for(std::vector<Branch>::iterator iter = gs.getBranches().begin(); iter != gs.getBranches().end(); iter++)
	{
		exploreBranch(*iter);
		cout << "Children: " << gs.getDescendantBranchIDs(iter->getBranchID()).first << "-" << gs.getDescendantBranchIDs(iter->getBranchID()).second << endl;
	}

	std::cout << "************************************" << std::endl;
	std::cout << "  Test 2: glycan composition        " << std::endl;
	std::cout << "************************************" << std::endl;
	GlycanComposition g_compo;
	g_compo.addGlycanComposition("DeltaGlcA", 1);
	g_compo.addGlycanComposition("GlcA", 2);
	g_compo.addGlycanComposition("GlcN", 3);
	g_compo.addModificationComposition("Ac", 1);
	g_compo.addModificationComposition("SO3", 5);

	GlycanSequence gs1;
	gs1.buildByGAGComposition(g_compo);

	for(std::vector<Branch>::iterator iter = gs1.getBranches().begin(); iter != gs1.getBranches().end(); iter++)
	{
		exploreBranch(*iter);
		cout << "Children: " << gs1.getDescendantBranchIDs(iter->getBranchID()).first << "-" << gs1.getDescendantBranchIDs(iter->getBranchID()).second << endl;
	}

	std::cout << "************************************" << std::endl;
	std::cout << "Test 3: Modification Site" << std::endl;
	std::cout << "************************************" << std::endl;
	ModificationByPosition& mod_pos_index = gs1.getModificationSetByPosition();
	gs1.modifyModificationStatus("SO3", ModificationPosition(0, 0, 2));

	gs1.modifyModificationStatus("Ac", ModificationPosition(0, 3, 2));
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

	std::cout << "************************************" << std::endl;
	std::cout << "     Test 3: Derivatization         " << std::endl;
	std::cout << "************************************" << std::endl;
	
	ModificationPosition derv(0, 5, 1);
	gs1.addModification("AnMan", derv);

	for(std::vector<Branch>::iterator iter = gs1.getBranches().begin(); iter != gs1.getBranches().end(); iter++)
	{
		exploreBranch(*iter);
		cout << "Children: " << gs1.getDescendantBranchIDs(iter->getBranchID()).first << "-" << gs1.getDescendantBranchIDs(iter->getBranchID()).second << endl;
	}



	cin.get();

	return EXIT_SUCCESS;
}