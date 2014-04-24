/*
 * =====================================================================================
 *
 *       Filename:  PeriodicTableTest.cpp
 *
 *    Description:  Test file for class Periodic
 *
 *        Version:  1.0
 *        Created:  04/26/2012  4:17:25 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#include <iostream>

#include <GAGPL/CHEMISTRY/PeriodicTable.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
using namespace std;
using namespace gag;

void printIsotope(const Isotope& iso) {
	cout << "Isotope information" << endl;
	cout << "Mass: " << iso.mass << endl;
	cout << "Abd: " << iso.abundance << endl;
	cout << "Shift: " << iso.neutron_shift << endl;
	cout << "Nominal mass: " << iso.nominal_mass << endl;
}
int main ( int argc, char *argv[] )
{



	PeriodicTable& pt = PeriodicTable::Instance();
	//pt.load();

	cout << "We are going to add Ca" << endl;

	Element e = pt.getElementBySymbol("Ca");
	
	cout << e.name << ' ' << e.isotopes.size() << endl;

	// Iterate over all isotopes for given element.
	IsotopeSetSequential::iterator it = e.isotopes.begin();
	for (; it != e.isotopes.end(); it++)
	{
		//cout << e.isotopes[i].mass << ": " << e.isotopes[i].abundance;
		//cout << " - " << e.isotopes[i].neuron_shift << endl;
		cout << "Mass: " << it->mass << endl;
		cout << "Abd: " << it->abundance << endl;
		cout << "Shift: " << it->neutron_shift << endl;
		cout << "Nominal mass: " << it->nominal_mass << endl;
	}

	// Test if the isotope for give shift or neutron number can be correctly identified.
	const Isotope isotope1 = pt.getIsotopeByRelativeShift("C", 1);
	printIsotope(isotope1);
	const Isotope isotope2 = pt.getIsotopeByNominalMass("O", 18);
	printIsotope(isotope2);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
