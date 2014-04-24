/*
 * =====================================================================================
 *
 *       Filename:  FunctionalGroupTable.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/29/2012 10:16:24 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <boost/lexical_cast.hpp>
#include <iostream>

namespace gag
{
	FunctionalGroupTable& FunctionalGroupTable::Instance()
	{
		static FunctionalGroupTable cgt;
		return cgt;
	}

	void FunctionalGroupTable::load(const std::string& filename)
	{
		std::cout << "Load functional group table once!\n" << std::endl;
		PeriodicTable& ptable = PeriodicTable::Instance();
		//ptable.load();

		using boost::property_tree::ptree;
		ptree pt;

		read_xml(filename, pt);

		FunctionalGroupReference fg_ref;

		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.FunctionalGroupSets"))
		{
			if(v.first == "FunctionalGroup")
			{
				FunctionalGroupProtoType fg_pt;

				fg_pt._symbol = v.second.get<std::string>("Symbol");			
				fg_pt._name = v.second.get<std::string>("Name");

				//const std::string compo_string = v.second.get<std::string>("Composition");
				//fg.getComposition().update(compo_string );
				fg_pt._composition_string = v.second.get<std::string>("Composition");
				
				if(v.second.count("Ring") > 0) {
					BOOST_FOREACH(ptree::value_type &m, v.second.get_child("Ring"))
					{
						if(m.first == "Start")
							fg_pt._start = boost::lexical_cast<size_t>(m.second.data());
						else if(m.first == "End")
							fg_pt._end = boost::lexical_cast<size_t>(m.second.data());
						else {
							// Do nothing.
						}
					}
				} else {
					fg_pt._start = 0;
					fg_pt._end = 0;
				}
				
				if(v.second.count("Sites")>0) // Dependent functional group.
				{
					BOOST_FOREACH(ptree::value_type &s, v.second.get_child("Sites"))
					{
						
						// Composition nodes.
						if(s.first == "Site")
						{
							// The repeating number of Site.
							size_t loop1 = s.second.get("<xmlattr>.Count", 1);
							
							// Check the contents of Site and decide if conversion is needed.
							// 1. A Site has to have a core or functional group as its children.
							// 2. If a core, store the information into _sites.
							// 3. If a functional group, covert it to Site which will be stored in _sites.
							std::string core;
							// tag Subset.FunctionalGroup or FunctionalGroup.
							std::vector<std::string> satellites;
							
							if(s.second.count("Core") > 0) {
								core = s.second.get<std::string>("Core");
								
								// Subset is not necessarily present.
								if(s.second.count("Subset") > 0) {
									BOOST_FOREACH(ptree::value_type &t, s.second.get_child("Subset"))
									{
										if(t.first == "FunctionalGroup") {
											size_t loop2 = t.second.get("<xmlattr>.Count", 1);	
											for(size_t i = 0; i<loop2; i++)
												satellites.push_back(t.second.data());
										}
									}
								}
							} else if(s.second.count("FunctionalGroup") == 1) {
								// Only one functional group is allowed for replacement.
								// core will be kept empty.
								satellites.push_back(s.second.get<std::string>("FunctionalGroup")); 	
							} else {
								std::cout << fg_pt._name << std::endl;
								throw std::runtime_error("Functional group for replacing is missing!");
							}

							for(size_t i = 0; i < loop1; i++) 
								fg_pt._sites.push_back(std::make_pair(core, satellites));

						} 
					}
				} else { // Independent functional group.
					// Currently, do nothing.
				}
				
				fg_ref.insert(std::make_pair(fg_pt._symbol, fg_pt));
			}
		}
		
		pt.erase("parameters.FunctionalGroupSets");
		this->buildtree(fg_ref);

	}

	// Convert information from fg_ref to functionalgroups.
	void FunctionalGroupTable::buildtree(const FunctionalGroupReference& fg_ref)
	{
		FunctionalGroupReference::const_iterator iter = fg_ref.begin();
		// Iterate over the reference table.
		for(; iter != fg_ref.end(); iter++) {
			// Create a FunctionalGroup object using symbol information from fg_ref.
			this->createFunctionalGroup(fg_ref, iter->first);
			//functionalgroups.insert(std::make_pair(iter->first, ff));
		}
	}

	FunctionalGroup FunctionalGroupTable::createFunctionalGroup(const FunctionalGroupReference& fg_ref, const std::string& symbol)
	{
		FunctionalGroupReference::const_iterator iter = fg_ref.find(symbol);
		if(iter == fg_ref.end())
			throw std::runtime_error("Symbol cannot be found in the reference table!");

		// If this functional group has been created previously, return it directly.
		if(this->containFunctionalGroup(symbol))
			return this->getFunctionalGroupBySymbol(symbol);

		// If not defined, create it from the scratch.
		FunctionalGroup ff(symbol, iter->second._name, iter->second._start, iter->second._end);

		// Set the composition of the functional group.
		ff.add(iter->second._composition_string);

		// Check if it is an independent functional group.
		if(!iter->second._sites.empty()) { // Dependent functional group.
			// Iterate over all the sites.
			std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator site_iter = iter->second._sites.begin();
			for(; site_iter != iter->second._sites.end(); site_iter++)
			{			
				if(site_iter->first.empty()){ // If core is empty. Convert the functional group to Site.

					if(site_iter->second.empty() || site_iter->second.size() > 1)
						throw std::runtime_error("There should be only one functional group for replacing site.");

					std::string child_symbol = site_iter->second.front();
					FunctionalGroup child_fg = this->createFunctionalGroup(fg_ref, child_symbol);
					Site temp_st = child_fg.getConvertedSite();

					// Update _sites.
					ff._sites.push_back(temp_st);				

				} else {
					// Create Site
					PeriodicTable& ptable = PeriodicTable::Instance();
					Site temp_st;
					temp_st.core = ptable.getElementBySymbol(site_iter->first);
					std::vector<std::string>::const_iterator fg_iter = site_iter->second.begin();
					for (; fg_iter != site_iter->second.end(); fg_iter++)
					{
						FunctionalGroup temp_fg = this->createFunctionalGroup(fg_ref, *fg_iter);
						std::string temp_symbol(*fg_iter);
						temp_st.fg_map.insert(std::make_pair(temp_symbol, temp_fg));
					}
					ff._sites.push_back(temp_st);
				}
			}

		} else { // Independent functional group.
				
		}

		functionalgroups.insert(std::make_pair(symbol, ff));		
		// Update functionalgroups.
		return ff;
	}

	FunctionalGroup FunctionalGroupTable::getFunctionalGroupBySymbol(const std::string& symbol) const
	{
		std::map<std::string, FunctionalGroup>::const_iterator i = functionalgroups.find(symbol);

		return i != functionalgroups.end() ? i->second : FunctionalGroup();
	}

	bool FunctionalGroupTable::containFunctionalGroup( const std::string& symbol )
	{
		std::map<std::string, FunctionalGroup>::iterator iter = functionalgroups.find(symbol);
		return iter != functionalgroups.end();
	}

	std::vector<std::string> FunctionalGroupTable::getAllFunctionalGroupSymbols() const
	{
		std::vector<std::string> keys;
		std::pair<std::string, FunctionalGroup> fg;
		BOOST_FOREACH(fg, functionalgroups)
		{
			keys.push_back(fg.first);
		}
		return keys;
	}

	//void FunctionalGroupTable::appendModificationToFunctionalGroup( FunctionalGroup& fg, Modification& mod)
	//{
	//	Reaction temp_rt;
	//	BOOST_FOREACH(temp_rt, mod.getModificationReactions())
	//	{
	//		// apply the corresponding reactions to the core.
	//		std::pair<AtomOperation, std::string> single_core_operation;
	//		BOOST_FOREACH(single_core_operation, temp_rt.core_operation)
	//		{
	//			if(single_core_operation.first == Addition) {
	//				FunctionalGroup fg_operation = this->getFunctionalGroupBySymbol(single_core_operation.second);
	//				fg.addFunctionalGroup(fg_operation, temp_rt.position);
	//			} else if(single_core_operation.first == Removement) {
	//				fg.removeFunctionalGroup(single_core_operation.second, temp_rt.position);
	//			}
	//		}
	//		
	//		// apply the corresponding reactions to the sub-functionalgroup.
	//		boost::tuple<size_t, AtomOperation, std::string> single_operation;
	//		std::pair<std::string, OperationSet> fg_operation_pair;
	//		BOOST_FOREACH(fg_operation_pair, temp_rt.sub_fg_operation)
	//		{				
	//			if(!fg.containFunctionalGroup(fg_operation_pair.first, temp_rt.position)) {
	//				continue;
	//			}
	//			FunctionalGroupChain chain;
	//			chain.push_back(std::make_pair(temp_rt.position, this->getFunctionalGroupBySymbol(fg_operation_pair.first)));
	//			BOOST_FOREACH(single_operation, fg_operation_pair.second)
	//			{
	//				FunctionalGroupChain chain_temp(chain);
	//				chain_temp.push_back(std::make_pair(single_operation.get<0>(), this->getFunctionalGroupBySymbol(single_operation.get<2>())));
	//				if(single_operation.get<1>() == Addition) {						
	//					//fg.getSubFunctionalGroup(fg_operation_pair.first, temp_rt.position).addFunctionalGroup(fg_operation, single_operation.get<0>());
	//					fg.addFunctionalGroupByChain(chain_temp);
	//				} else if(single_operation.get<1>() == Removement) {
	//					//fg.getSubFunctionalGroup(fg_operation_pair.first, temp_rt.position).removeFunctionalGroup(single_operation.get<2>(), single_operation.get<0>());
	//					fg.removeFunctionalGroupByChain(chain_temp);
	//				}
	//			}
	//		}

	//		// Get the reactant functional group and modify it.
	//		FunctionalGroup fg_ext = this->getFunctionalGroupBySymbol(temp_rt.reactant_operation.first);
	//		BOOST_FOREACH(single_operation, temp_rt.reactant_operation.second)
	//		{
	//			if(single_operation.get<1>() == Addition) {
	//				FunctionalGroup fg_operation = this->getFunctionalGroupBySymbol(single_operation.get<2>());
	//				fg_ext.addFunctionalGroup(fg_operation, single_operation.get<0>());
	//			} else if(single_operation.get<1>() == Removement) {
	//				fg_ext.removeFunctionalGroup(single_operation.get<2>(), single_operation.get<0>());
	//			}
	//		}

	//		if(!temp_rt.core_operation.empty()) { // Add to the core.
	//			fg.addFunctionalGroup(fg_ext, temp_rt.position);
	//		} else if(!temp_rt.sub_fg_operation.empty()) { // Add to the sub functional group.
	//			// Check and see if the functional group contain such a sub functional group.
	//			std::multimap<std::string, OperationSet>::iterator oper_iter = temp_rt.sub_fg_operation.begin();
	//			for(; oper_iter != temp_rt.sub_fg_operation.end(); oper_iter++)
	//			{
	//				if(fg.containFunctionalGroup(oper_iter->first, temp_rt.position)) {
	//					FunctionalGroupChain chain;
	//					chain.push_back(std::make_pair(temp_rt.position, this->getFunctionalGroupBySymbol(oper_iter->first)));
	//					chain.push_back(std::make_pair(0, fg_ext));
	//					//fg.getSubFunctionalGroup(oper_iter->first, temp_rt.position).addFunctionalGroup(fg_ext);
	//					fg.addFunctionalGroupByChain(chain);
	//				}
	//			}
	//		}
	//	}
	//}

	//void FunctionalGroupTable::appendModificationToFunctionalGroup( FunctionalGroup& fg, Modification mod, const size_t idx )
	//{
	//	// Use this function at risk.
	//	for(std::vector<Reaction>::iterator iter = mod.getModificationReactions().begin(); iter != mod.getModificationReactions().end(); iter++)
	//	{
	//		iter->position = idx;
	//	}
	//		
	//	this->appendModificationToFunctionalGroup(fg, mod);
	//		
	//}
}
