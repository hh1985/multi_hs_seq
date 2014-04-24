/*
 * =====================================================================================
 *
 *       Filename:  MonosaccharideUnitTable.cpp
 *
 *    Description:  The implementation file for MonosaccharideUnitTable loaded from 
 *    							xml file.
 *
 *        Version:  1.0
 *        Created:  04/30/2012 11:46:16 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <iostream>

namespace gag
{
	MonosaccharideUnitTable& MonosaccharideUnitTable::Instance()
	{
		static MonosaccharideUnitTable mut;
		return mut;
	}

	void MonosaccharideUnitTable::load(const std::string& filename)
	{
		std::cout << "Load monosaccharide unit table once!" << std::endl;
		using boost::property_tree::ptree;
		ptree pt;

		read_xml(filename, pt);

		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.MonosaccharideUnitSet"))
		{
			if(v.first == "Monosaccharide")
			{
				size_t start = v.second.get<size_t>("Ring.Start");
				size_t end = v.second.get<size_t>("Ring.End");
				Monosaccharide ms(v.second.get<std::string>("Composition"),v.second.get<std::string>("Symbol"),v.second.get<std::string>("Name"), start, end);

				//ms(v.second.get<std::string>("Composition"));

				BOOST_FOREACH(ptree::value_type &s, v.second.get_child("Sites"))
				{
					std::string& compo = s.second.get<std::string>("Composition");
					
					Site is(compo);

					if(s.second.count("Subset")>0) {
						BOOST_FOREACH(ptree::value_type &m, s.second.get_child("Subset"))
						{
							FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
							//fgt.load();
							FunctionalGroup& ff = fgt.getFunctionalGroupBySymbol(m.second.data());
							is.addFunctionalGroup(ff);
						}
					}	

					ms.addInternalSite(is);
					
				}
				monos.insert(std::make_pair(ms.getSymbol(), ms));
			}
		}

		pt.erase("parameters.MonosaccharideUnitSet");
	}

	Monosaccharide MonosaccharideUnitTable::getMonosaccharideBySymbol(const std::string& symbol)
	{
		std::map<std::string, Monosaccharide>::iterator iter = 
			monos.find(symbol);

		if(iter != monos.end())
			return iter->second;
		else 
			return Monosaccharide();
	}
}
