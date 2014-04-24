/*
 * =====================================================================================
 *
 *       Filename:  FunctionalGroupTableTest.cpp
 *
 *    Description:  Test file for class FunctionalGroupTable.
 *
 *        Version:  1.0
 *        Created:  04/30/2012  1:39:20 AM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#include <iostream>

#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>

int main ( int argc, char *argv[] )
{

	using namespace std;
	using namespace gag;

	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();

	// Go through all the functional groups.
	std::vector<std::string> fg_vec = fgt.getAllFunctionalGroupSymbols();
	for(std::vector<std::string>::iterator iter = fg_vec.begin(); iter != fg_vec.end(); iter++) {
		FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol(*iter);
		fg.printFunctionalGroupInformation();
	}

	// Do some common modification.
	// 1. N-S
	FunctionalGroup fg_nh2 = fgt.getFunctionalGroupBySymbol("NH2");
	FunctionalGroup fg_s = fgt.getFunctionalGroupBySymbol("H2SO4");
	FunctionalGroup fg_h = fgt.getFunctionalGroupBySymbol("H");
	FunctionalGroup fg_oh = fgt.getFunctionalGroupBySymbol("OH");
	fg_nh2.removeFunctionalGroup(fg_h);
	fg_s.removeFunctionalGroup(fg_oh);
	fg_nh2.addFunctionalGroup(fg_s);
	cout << "Modified N-S: " << endl;
	fg_nh2.printFunctionalGroupInformation();

	// 2. O-S
	FunctionalGroup fg_oh_1 = fgt.getFunctionalGroupBySymbol("OH");
	FunctionalGroup fg_s_1 = fgt.getFunctionalGroupBySymbol("H2SO4");
	fg_oh_1.removeFunctionalGroup(fg_h);
	fg_s_1.removeFunctionalGroup(fg_oh);
	fg_oh_1.addFunctionalGroup(fg_s_1);
	cout << "Modified 2-OS: " << endl;
	fg_oh_1.printFunctionalGroupInformation();

	// 3. N-Ac
	FunctionalGroup fg_nh2_1 = fgt.getFunctionalGroupBySymbol("NH2");
	FunctionalGroup fg_ac = fgt.getFunctionalGroupBySymbol("AcOH");
	fg_nh2_1.removeFunctionalGroup(fg_h);
	fg_ac.removeFunctionalGroup(fg_oh);
	fg_nh2_1.addFunctionalGroup(fg_ac);
	cout << "Modified NAc: " << endl;
	fg_nh2_1.printFunctionalGroupInformation();

	FunctionalGroup fg_mono = fgt.getFunctionalGroupBySymbol("GlcN");
	fg_mono.addFunctionalGroup(fg_h);
	fg_mono.removeFunctionalGroup(fg_oh, 1);
	fg_mono.addFunctionalGroup(fg_h, 1);

	//FunctionalGroup fg_2ab = fgt.getFunctionalGroupBySymbol("2-AB");
	//
	////FunctionalGroup fg_nh_2 = fg_2ab.reorganizeSites("NH2", 1);
	////fg_nh_2.removeFunctionalGroup(fg_h);
	//fg_2ab.reorganizeSites("NH2",1);
	//fg_2ab.removeFunctionalGroup(fg_h);

	//fg_mono.printFunctionalGroupInformation();
	//fg_2ab.printFunctionalGroupInformation();
	//fg_nh_2.printFunctionalGroupInformation();

	FunctionalGroup fg_mono_1 = fgt.getFunctionalGroupBySymbol("GlcN");

	FunctionalGroupChain chain1;
	chain1.push_back(std::make_pair(2, fg_nh2));
	chain1.push_back(std::make_pair(0, fg_h));
	fg_mono_1.replaceFunctionalGroupByChain(chain1, fg_s_1);

	cout << "Modification by Chain..." << endl;
	fg_mono_1.printFunctionalGroupInformation();

	cout << fg_mono_1.getSubCompositionByRingID(0,2).getCompositionString() << endl;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
