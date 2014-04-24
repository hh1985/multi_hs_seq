/*
 * =====================================================================================
 *
 *       Filename:  PeriodicTable.cpp
 *
 *    Description:  The specific implementation for loading periodic table data.  
 *
 *        Version:  1.0
 *        Created:  4/24/2012 3:06:40 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include "GAGPL/CHEMISTRY/PeriodicTable.h"
#include <iostream>

namespace gag
{

	PeriodicTable& PeriodicTable::Instance()
	{
		static PeriodicTable pt;
		return pt;
	}

	void PeriodicTable::load(const std::string& filename)
	{
		std::cout << "Load periodic table !" << std::endl;
		using boost::property_tree::ptree;
		ptree pt;
		
		read_xml(filename, pt);

		//ptree es = pt.get_child("parameters.ElementIsotopes");

		//std::pair<ptree::assoc_iterator, ptree::assoc_iterator> pcld = 
		//	es.equal_range("Element");
		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.ElementIsotopes"))
		{
			//std::cout << v.first << std::endl;
			if(v.first == "Element")
			{
				Element elem;
				elem.symbol = v.second.get<std::string>("Symbol");
				elem.name = v.second.get<std::string>("Name");
				elem.atomicity = v.second.get<int>("Atomicity");
				//std::cout << "Before Isotope" << std::endl;
				std::pair<ptree::assoc_iterator, ptree::assoc_iterator> ei = v.second.equal_range("Isotope");
				for(ptree::assoc_iterator iter = ei.first; iter != ei.second; iter++)
				{
					// Notice that this round method only works for positive number.
					size_t nominal_mass = int(iter->second.get<double>("Mass") + 0.5);
					
					Isotope iso = {iter->second.get<double>("Mass"), iter->second.get<float>("Probability"), iter->second.get<size_t>("Shift"), nominal_mass};

					//std::cout << iso.mass << " " << iso.abundance << std::endl;
					elem.isotopes.push_back(iso);					
				}

				elements.insert(std::make_pair(elem.symbol, elem));
			}
		}	
		// Delete params and keep all information only in elements.
		pt.erase("parameters.ElementIsotopes");
	
	}


	const Element PeriodicTable::getElementBySymbol(const std::string& symbol) const
	{
		std::map<std::string, Element>::const_iterator i = elements.find(symbol);
		return i != elements.end() ? i->second : Element();
	}

	const Isotope PeriodicTable::getIsotopeByRelativeShift( const std::string& symbol, const size_t shift/*=0*/ ) const
	{
		const Element ele = this->getElementBySymbol(symbol);
		// Nominal shift.
		const IsotopeSet::nth_index<1>::type& isotope_index = ele.isotopes.get<1>();
		IsotopeSet::nth_index<1>::type::const_iterator shift_iter = isotope_index.find(shift);
		return shift_iter != isotope_index.end() ? *shift_iter : Isotope();
	}

	const Isotope PeriodicTable::getIsotopeByNominalMass( const std::string& symbol, const size_t mass ) const
	{
		const Element ele = this->getElementBySymbol(symbol);
		// Nominal shift.
		const IsotopeSet::nth_index<2>::type& isotope_index = ele.isotopes.get<2>();
		IsotopeSet::nth_index<2>::type::const_iterator shift_iter = isotope_index.find(mass);
		return shift_iter != isotope_index.end() ? *shift_iter : Isotope();

	}

	size_t PeriodicTable::getMaxRelativeShift( const std::string& symbol ) const
	{
		const Element ele = this->getElementBySymbol(symbol);

		IsotopeSetSequential::reverse_iterator it = ele.isotopes.rbegin();

		return it->neutron_shift;
	}

}

