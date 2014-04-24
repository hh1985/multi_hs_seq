/********************************************************************
	created:	2012/06/27
	created:	27:6:2012   15:55
	filename: 	ModificationTable.cpp
	file path:	GAG\src\GAGPL\CHEMISTRY
	file base:	ModificationTable
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include <GAGPL/CHEMISTRY/ModificationTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <boost/tuple/tuple.hpp>
#include <boost/make_shared.hpp>
#include <iostream>

namespace gag
{
	ModificationTable& ModificationTable::Instance()
	{
		static ModificationTable mt;
		return mt;
	}

	void ModificationTable::load(const std::string& filename)
	{
		std::cout << "Load modification once!" << std::endl;
		FunctionalGroupTable& ftable = FunctionalGroupTable::Instance();
		//ftable.load();

		using boost::property_tree::ptree;
		ptree pt;

		read_xml(filename, pt);

		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.ModificationSet"))
		{
			if(v.first == "Modification") {
				// Only one name is allowed.
				std::string name = v.second.get<std::string>("Name");
				// Multiple symbols are allowed.
				std::set<std::string> symbols;
				std::vector<Reaction> reactions;
				
				std::pair<ptree::assoc_iterator, ptree::assoc_iterator> ei = v.second.equal_range("Symbol");
				for(ptree::assoc_iterator iter = ei.first; iter != ei.second; iter++)
				{
					symbols.insert(iter->second.data());
				}

				// Iterate over all the reactions.
				std::pair<ptree::assoc_iterator, ptree::assoc_iterator> reaction_pair = v.second.equal_range("Reaction");
				for(ptree::assoc_iterator iter = reaction_pair.first; iter != reaction_pair.second; iter++)
				{
					Reaction rt;
					
					if(iter->second.count("Position") > 0)
						rt.position = iter->second.get<size_t>("Position");

					BOOST_FOREACH(ptree::value_type &s, iter->second.get_child("Target"))
					{
						if(s.first == "Core") {

							if(s.second.count("Minus") > 0) {
								
								std::pair<ptree::assoc_iterator, ptree::assoc_iterator> minus_pair = s.second.equal_range("Minus");
								for(ptree::assoc_iterator assoc_iter = minus_pair.first; assoc_iter != minus_pair.second; assoc_iter++)
								{
									rt.core_operation.push_back(std::make_pair(Removement, assoc_iter->second.data()));
								}
							}
							if(s.second.count("Plus") > 0) {
								//size_t position = s.second.get("Plus.<xmlattr>.Site", 0);
								std::pair<ptree::assoc_iterator, ptree::assoc_iterator> plus_pair = s.second.equal_range("Plus");
								for(ptree::assoc_iterator assoc_iter = plus_pair.first; assoc_iter != plus_pair.second; assoc_iter++)
								{
									rt.core_operation.push_back(std::make_pair(Addition, assoc_iter->second.data()));
								}
							}
						} else if(s.first == "Subgroup") {
							
							OperationSet oper_set;

							if(s.second.count("Minus") > 0) {
								size_t fg_position = s.second.get("Minus.<xmlattr>.Site", 0);
								std::pair<ptree::assoc_iterator, ptree::assoc_iterator> minus_pair = s.second.equal_range("Minus");
								for(ptree::assoc_iterator assoc_iter = minus_pair.first; assoc_iter != minus_pair.second; assoc_iter++)
								{
									oper_set.push_back(boost::make_tuple(fg_position, Removement, assoc_iter->second.data()));
								}
							} 
							if(s.second.count("Plus") > 0) {
								size_t fg_position = s.second.get("Plus.<xmlattr>.Site",0);
								std::pair<ptree::assoc_iterator, ptree::assoc_iterator> plus_pair = s.second.equal_range("Plus");
								for(ptree::assoc_iterator assoc_iter = plus_pair.first; assoc_iter != plus_pair.second; assoc_iter++)
								{
									oper_set.push_back(boost::make_tuple(fg_position, Addition, assoc_iter->second.data()));
								}
							}
							rt.sub_fg_operation.insert(std::make_pair(s.second.get<std::string>("FunctionalGroup"), oper_set));
						}
						

					}

					OperationSet oper_set;
					if(iter->second.count("Reactant") > 0) {
						BOOST_FOREACH(ptree::value_type &s, iter->second.get_child("Reactant"))
						{
							if(s.first == "FunctionalGroup")
								rt.reactant_operation.first = s.second.data();
							else if(s.first ==  "Minus") {
								size_t fg_position = s.second.get("<xmlattr>.Site", 0);
								oper_set.push_back(boost::make_tuple(fg_position, Removement, s.second.data()));
							} else if(s.first == "Plus") {
								size_t fg_position = s.second.get("<xmlattr>.Site", 0);
								oper_set.push_back(boost::make_tuple(fg_position, Addition, s.second.data()));
							}
						}
					}
					
					rt.reactant_operation.second = oper_set;

					reactions.push_back(rt);
				}

				ModificationRule mr;
				// Iterate over all the rules.
				std::pair<ptree::assoc_iterator, ptree::assoc_iterator> rule_pair = v.second.equal_range("Rule");
				for(ptree::assoc_iterator rule_iter = rule_pair.first; rule_iter != rule_pair.second; rule_iter++)
				{	
					std::string fg_symbol = rule_iter->second.get<std::string>("FunctionalGroup");

					std::vector<size_t> sites;
					std::pair<ptree::assoc_iterator, ptree::assoc_iterator> site_pair = rule_iter->second.equal_range("Position");
					for(ptree::assoc_iterator site_iter = site_pair.first; site_iter != site_pair.second; site_iter++)
					{
						sites.push_back((size_t)atoi(site_iter->second.data().c_str()));
					}

					mr.insert(std::make_pair(fg_symbol, sites));
				}

				// Store the information into mod_table;
				ModificationPtr mod_ptr = boost::make_shared<Modification>(name, symbols, reactions, mr);
				// mod_table.insert(mod_ptr);
				// Update the reference table: mod_by_name and mod_by_shortcut.
				mod_by_name.insert(std::make_pair(name, mod_ptr));
				for(std::set<std::string>::iterator iter = mod_ptr->getSymbols().begin(); iter!=mod_ptr->getSymbols().end(); iter++)
				{
					mod_by_symbols.insert(std::make_pair(*iter, mod_ptr));
				}




			}
		}

	}

	Modification ModificationTable::getModificationBySymbol( const std::string& symbol )
	{
		std::map<std::string, ModificationPtr>::iterator iter = mod_by_symbols.find(symbol);
		return iter != mod_by_symbols.end() ? *(iter->second) : Modification();
	}

	Modification ModificationTable::getModificationByName( const std::string& name)
	{
		std::map<std::string, ModificationPtr>::iterator iter = mod_by_name.find(name);
		return iter != mod_by_name.end() ? *(iter->second) : Modification();
	}




}