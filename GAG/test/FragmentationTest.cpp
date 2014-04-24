#include <GAGPL/GLYCAN/GlycanSequence.h>
//#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/FRAGMENTATION/FragmentationTable.h>
#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <GAGPL/CHEMISTRY/Modification.h>
#include <GAGPL/MISC/Param.h>
#include <iostream>
#include <cmath>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>

using namespace std;
using namespace gag;
using namespace param;

int main()
{
	// The test is only on single fragments.
	std::cout << "**************************************" << std::endl;
	std::cout << "         Test 1: Fragmentation        " << std::endl;
	std::cout << "**************************************" << std::endl;

	// For the overall information. Refer to class library
	GlycanComposition g_compo;
	g_compo.addGlycanComposition("DeltaGlcA", 1);
	g_compo.addGlycanComposition("GlcA", 2);
	g_compo.addGlycanComposition("GlcN", 3);
	g_compo.addModificationComposition("Ac", 1);
	g_compo.addModificationComposition("SO3", 6);

	GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
	gs->buildByGAGComposition(g_compo);

	// Create the fragment.
	Fragment fg(gs);
	CleavageCollection cc;
	// Set the fragment type: Y5 A5(1,5)
	FragmentPosition fp1(0,0,0,0);
	FragmentPosition fp2(0,4,1,5);

	cc.insert(std::make_pair("Y", fp1));
	cc.insert(std::make_pair("A", fp2));
	//fg.setFragmentation("Y", fp1);
	//fg.setFragmentation("A", fp2);
	fg.updateCleavage(cc);
	fg.printFragment();

	Fragment new_fg(gs);
	CleavageCollection cc1;
	cc1.insert(std::make_pair("B", FragmentPosition(0, 0, 0, 0)));
	new_fg.updateCleavage(cc1);
	//new_fg.setFragmentation("X", 5, "", 0, 2);
	new_fg.printFragment();

	Fragment ad_fg(gs);
	CleavageCollection cc2;
	cc2.insert(std::make_pair("A", FragmentPosition(0, 2, 2, 4)));
	cc2.insert(std::make_pair("X", FragmentPosition(0, 1, 2, 4)));
	ad_fg.updateCleavage(cc2);
	ad_fg.printFragment();
	ModificationSites mod_sites = ad_fg.getModificationSitesBySymbol("Ac", 1);

	std::cout << "**************************************" << std::endl;
	std::cout << "         Test 2: Modification         " << std::endl;
	std::cout << "**************************************" << std::endl;

	ModificationByPosition& mod_pos_index = fg.getModificationSetByPosition();
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

	std::cout << std::endl;
	std::cout << "The total available modification sites: " << std::endl;
	
	std::set<std::string> mod_types = gs->getModificationTypes();
	for(std::set<std::string>::iterator iter = mod_types.begin(); iter != mod_types.end(); iter++)
	{
		std::cout << *iter << ": " << fg.getModificationSiteNum(*iter, 1) << std::endl;
	}

	std::cout << "**************************************" << std::endl;
	std::cout << "    Test 3: Multiple fragmentation    " << std::endl;
	std::cout << "**************************************" << std::endl;

	// Terminal cleavage
	std::cout << "Test fragment type B" << std::endl;
	for(size_t i = 1; i <=6; i++)
	{
		//std::cout << "Fragment B" << i << ":" << std::endl;
		Fragment fg1(gs);
		fg1.setFragmentation("B", i, "", 0, 0);
		std::cout << fg1.getCleavageType() << std::endl;
		//printFragment(fg1);
		fg1.printFragment();
	}
	std::cout << "Test fragment type C" << std::endl;
	for(size_t i = 1; i <=6; i++)
	{
		//std::cout << "Fragment C" << i << ":" << std::endl;
		Fragment fg1(gs);
		fg1.setFragmentation("C", i, "", 0, 0);
		std::cout << fg1.getCleavageType() << std::endl;
		//printFragment(fg1);
		fg1.printFragment();
	}
	std::cout << "Test fragment type A" << std::endl;
	for(size_t i = 1; i<=6; i++)
	{
		for(size_t x = 0; x<=3; x++)
			for(size_t y = x+2; y-x<5 && y<=5; y++)
			{
				//std::cout << "Fragment A" << i << ": " << x << "-" << y << std::endl;
				if(x == 0 && y == 4) continue;
				if(x == 1 && y == 3) continue;
				Fragment fg1(gs);
				fg1.setFragmentation("A", i, "", x, y);
				std::cout << fg1.getCleavageType() << std::endl;
				//printFragment(fg1);
				fg1.printFragment();
			}	
	}
	std::cout << "Test fragment type Y" << std::endl;
	for(size_t i = 0; i <=5; i++)
	{
		//std::cout << "Fragment Y" << i << ":" << std::endl;
		Fragment fg1(gs);
		fg1.setFragmentation("Y", i, "", 0, 0);
		std::cout << fg1.getCleavageType() << std::endl;
		//printFragment(fg1);
		fg1.printFragment();
	}
	std::cout << "Test fragment type Z" << std::endl;
	for(size_t i = 0; i <=5; i++)
	{
		//std::cout << "Fragment Z" << i << ":" << std::endl;
		Fragment fg1(gs);
		fg1.setFragmentation("Z", i, "", 0, 0);
		std::cout << fg1.getCleavageType() << std::endl;
		//printFragment(fg1);
		fg1.printFragment();
	}
	std::cout << "Test fragment type X" << std::endl;
	for(size_t i = 0; i<=5; i++)
	{
		for(size_t x = 0; x<=3; x++)
			for(size_t y = x+2; y-x<5 && y<=5; y++)
			{
				//std::cout << "Fragment X" << i << ": " << x << "-" << y << std::endl;
				Fragment fg1(gs);
				fg1.setFragmentation("X", i, "", x, y);
				std::cout << fg1.getCleavageType() << std::endl;
				//printFragment(fg1);
				fg1.printFragment();
			}	
	}


  std::cout << "**************************************" << std::endl;
  std::cout << "    Test 4: Fragment reflection       " << std::endl;
  std::cout << "**************************************" << std::endl;
  
  
  int chain_length = gs->getBranchByID(0).getUnitNum();
  std::cout << "Chain length: " << chain_length << "\n";

  // Normal one.
  Fragment virtual_frag(gs);
  virtual_frag.setFragmentation("A", 3, 0, 2);
  std::cout << "Fragment type: " << virtual_frag.getFragmentType() << "\n";
  std::cout << "Cleavage index: " << virtual_frag.getCleavageIndex() << "\n";
  //printFragment(fg1);
  virtual_frag.printFragment();

  ModificationSites sulfate_sites = virtual_frag.getModificationSitesBySymbol("SO3",1);
  printModificationSites(sulfate_sites);
  ModificationSites expanded_sites = virtual_frag.getExpandedModificationSites("SO3",1);
  printModificationSites(expanded_sites);
  ModificationSites reduced_sites = virtual_frag.getReducedModificationSites("SO3", 1);
  printModificationSites(reduced_sites);

  Fragment virtual_frag1(gs);
  virtual_frag1.setFragmentation("X", 2, 3, 5);
  std::cout << "Fragment type: " << virtual_frag1.getFragmentType() << "\n";
  std::cout << "Cleavage index: " << virtual_frag1.getCleavageIndex() << "\n";
  //printFragment(fg1);
  virtual_frag1.printFragment();

  sulfate_sites = virtual_frag1.getModificationSitesBySymbol("SO3",1);
  printModificationSites(sulfate_sites);
  expanded_sites = virtual_frag1.getExpandedModificationSites("SO3",1);
  printModificationSites(expanded_sites);
  reduced_sites = virtual_frag1.getReducedModificationSites("SO3", 1);
  printModificationSites(reduced_sites);
  

  std::cout << "\nExtreme case:\n";
  // Extreme case.
  Fragment extreme_frag1(gs);
  extreme_frag1.setFragmentation("A", chain_length, 0, 2);
  extreme_frag1.printFragment();
  std::cout << "Modification sites:\n"; 
  printModificationSites(extreme_frag1.getModificationSitesBySymbol("SO3", 1));
  std::cout << "Expanded sites:\n";
  printModificationSites(extreme_frag1.getExpandedModificationSites("SO3", 1));
  std::cout << "Reduced sites:\n";
  printModificationSites(extreme_frag1.getReducedModificationSites("SO3", 1));

  Fragment extreme_frag2(gs);
  extreme_frag2.setFragmentation("X", chain_length-1, 1, 5);
  extreme_frag2.printFragment();
  std::cout << "Modification sites:\n";
  printModificationSites(extreme_frag2.getModificationSitesBySymbol("SO3", 1));
  std::cout << "Expanded sites:\n";
  printModificationSites(extreme_frag2.getExpandedModificationSites("SO3", 1));
  std::cout << "Reduced sites:\n";
  printModificationSites(extreme_frag2.getReducedModificationSites("SO3", 1));
  
  Fragment extreme_frag3(gs);
  extreme_frag3.setFragmentation("A", 1, 1, 5);
  extreme_frag3.printFragment();
  std::cout << "Modification sites:\n";
  printModificationSites(extreme_frag3.getModificationSitesBySymbol("SO3", 1));
  std::cout << "Expanded sites:\n";
  printModificationSites(extreme_frag3.getExpandedModificationSites("SO3", 1));
  std::cout << "Reduced sites:\n";
  printModificationSites(extreme_frag3.getReducedModificationSites("SO3", 1));

  Fragment extreme_frag4(gs);
  extreme_frag4.setFragmentation("X", 0, 0, 2);
  extreme_frag4.printFragment();
  std::cout << "Modification sites:\n";
  printModificationSites(extreme_frag4.getModificationSitesBySymbol("SO3", 1));
  std::cout << "Expanded sites:\n";
  printModificationSites(extreme_frag4.getExpandedModificationSites("SO3", 1));
  std::cout << "Reduced sites:\n";
  printModificationSites(extreme_frag4.getReducedModificationSites("SO3", 1));

	std::cout << std::endl;
	std::cout << "Done!" << std::endl;

	//cin.get();

	// Match against the library.
	return EXIT_SUCCESS;
}