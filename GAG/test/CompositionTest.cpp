/*
 * =====================================================================================
 *
 *       Filename:  CompositionTest.cpp
 *
 *    Description:  Test file for class Composition.
 *
 *        Version:  1.0
 *        Created:  04/28/2012 12:26:57 AM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/CHEMISTRY/Composition.h>
#include <iostream>

int main ( int argc, char *argv[] )
{
	using namespace gag;
	// The Periodic Table should be loaded at least once;
	PeriodicTable& pt = PeriodicTable::Instance();
	//pt.load();

	/* Set composition. */
	Composition cp("C7H20OHCl3");
	
	const std::map<std::string, int>& cp_map = cp.get();

	for(std::map<std::string, int>::const_iterator iter = cp_map.begin(); iter != cp_map.end(); iter++)
	{
		std::cout << iter->first << " " << iter->second << std::endl;
	}

	std::cout << "Mass: " << cp.getMass() << std::endl;
	std::cout << "size is: " << cp.getCompositionString().length() << std::endl;
	
	/* Add composition. */
	cp.add("CCa2");
	std::cout << "String: " << cp.getCompositionString() << std::endl;
	
	//cp_map = cp.get();
	for(std::map<std::string, int>::const_iterator iter = cp_map.begin(); iter != cp_map.end(); iter++)
	{
		std::cout << iter->first << " " << iter->second << std::endl;
	}

	std::cout << "Mass: " << cp.getMass() << std::endl;
	std::cout << "size is: " << cp.getCompositionString().length() << std::endl;
	/* Substract composition */
	cp.deduct("COH");
	std::cout << "String: " << cp.getCompositionString() << std::endl;
	
	//cp_map = cp.get();
	for(std::map<std::string, int>::const_iterator iter = cp_map.begin(); iter != cp_map.end(); iter++)
	{
		std::cout << iter->first << " " << iter->second << std::endl;
	}

	std::cout << "Mass: " << cp.getMass() << std::endl;
	std::cout << "size is: " << cp.getCompositionString().length() << std::endl;
	
	/* Assignment */
	Composition cm = cp;
	std::cout << "CP String: " << cp.getCompositionString() << std::endl;
	const std::map<std::string, int>& cm_map = cm.get();

	for(std::map<std::string, int>::const_iterator iter = cp_map.begin(); iter != cp_map.end(); iter++)
	{
		std::cout << iter->first << " " << iter->second << std::endl;
	}

	std::cout << "CP Mass: " << cp.getMass() << std::endl;
	std::cout << "CP size is: " << cp.getCompositionString().length() << std::endl;
	
	/* clear */
	cp.clear();
	std::cout << "After clear, size is: " << cp.getCompositionString().length() << std::endl;

	std::cout << "**************************************" << std::endl;
	std::cout << "         Test 2: Comparison            " << std::endl;
	std::cout << "**************************************" << std::endl;

	Composition comp0("C3HO3");
	Composition comp1("H2O");
	std::cout << "Compare" << comp0.getCompositionString() << " and " << comp1.getCompositionString() << std::endl;
	std::cout << "Results: ";
	if(comp1 < comp0)
		std::cout << "smaller" << std::endl;
	else
		std::cout << "larger" << std::endl;


	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
