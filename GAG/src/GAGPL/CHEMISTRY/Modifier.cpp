#include <GAGPL/CHEMISTRY/Modification.h>
#include <GAGPL/CHEMISTRY/Modifier.h>
#include <iostream>

namespace gag
{
	void Modifier::modifyFunctionalGroup( FunctionalGroup& fg, const Modification& mod)
	{
		//Reaction temp_rt;
		BOOST_FOREACH(Reaction temp_rt, mod.getModificationReactions())
		{
			// apply the corresponding reactions to the core.
			std::pair<AtomOperation, std::string> single_core_operation;
			BOOST_FOREACH(single_core_operation, temp_rt.core_operation)
			{
				if(single_core_operation.first == Addition) {
					FunctionalGroup fg_operation = fun_table.getFunctionalGroupBySymbol(single_core_operation.second);
					fg.addFunctionalGroup(fg_operation, temp_rt.position);
				} else if(single_core_operation.first == Removement) {
					fg.removeFunctionalGroup(single_core_operation.second, temp_rt.position);
				}
			}

			// apply the corresponding reactions to the sub-functionalgroup.
			boost::tuple<size_t, AtomOperation, std::string> single_operation;
			std::pair<std::string, OperationSet> fg_operation_pair;
			BOOST_FOREACH(fg_operation_pair, temp_rt.sub_fg_operation)
			{				
				if(!fg.containFunctionalGroup(fg_operation_pair.first, temp_rt.position)) {
					continue;
				}

				FunctionalGroupChain chain;
				chain.push_back(std::make_pair(temp_rt.position, fun_table.getFunctionalGroupBySymbol(fg_operation_pair.first)));
				BOOST_FOREACH(single_operation, fg_operation_pair.second)
				{
					FunctionalGroupChain chain_temp(chain);
					chain_temp.push_back(std::make_pair(single_operation.get<0>(), fun_table.getFunctionalGroupBySymbol(single_operation.get<2>())));
					if(single_operation.get<1>() == Addition) {						
						//fg.getSubFunctionalGroup(fg_operation_pair.first, temp_rt.position).addFunctionalGroup(fg_operation, single_operation.get<0>());
						fg.addFunctionalGroupByChain(chain_temp);
					} else if(single_operation.get<1>() == Removement) {
						//fg.getSubFunctionalGroup(fg_operation_pair.first, temp_rt.position).removeFunctionalGroup(single_operation.get<2>(), single_operation.get<0>());
						fg.removeFunctionalGroupByChain(chain_temp);
					}
				}
			}

			// Get the reactant functional group and modify it.
			FunctionalGroup fg_ext = fun_table.getFunctionalGroupBySymbol(temp_rt.reactant_operation.first);
			BOOST_FOREACH(single_operation, temp_rt.reactant_operation.second)
			{
				if(single_operation.get<1>() == Addition) {
					FunctionalGroup fg_operation = fun_table.getFunctionalGroupBySymbol(single_operation.get<2>());
					fg_ext.addFunctionalGroup(fg_operation, single_operation.get<0>());
				} else if(single_operation.get<1>() == Removement) {
					fg_ext.removeFunctionalGroup(single_operation.get<2>(), single_operation.get<0>());
				}
			}

			// Append the reactant functional group to the core or sub-functional group.
			if(!temp_rt.core_operation.empty()) { // Add to the core.
				fg.addFunctionalGroup(fg_ext, temp_rt.position);
			} else if(!temp_rt.sub_fg_operation.empty()) { 
				// Add to the sub functional group.
				// Check and see if the functional group contain such a sub functional group.
				std::multimap<std::string, OperationSet>::iterator oper_iter = temp_rt.sub_fg_operation.begin();
				for(; oper_iter != temp_rt.sub_fg_operation.end(); oper_iter++)
				{
					if(fg.containFunctionalGroup(oper_iter->first, temp_rt.position)) {
						FunctionalGroupChain chain;
						chain.push_back(std::make_pair(temp_rt.position, fun_table.getFunctionalGroupBySymbol(oper_iter->first)));
						chain.push_back(std::make_pair(0, fg_ext));

						fg.addFunctionalGroupByChain(chain);
					}
				}
			}
		}
	}

	void Modifier::modifyFunctionalGroup( FunctionalGroup& fg, const std::string& mod_symbol)
	{
		Modification mod = mod_table.getModificationBySymbol(mod_symbol);

		this->modifyFunctionalGroup(fg, mod);
	}

	void Modifier::modifyFunctionalGroup( FunctionalGroup& fg, const std::string& mod_symbol, const size_t index )
	{
		Modification mod = mod_table.getModificationBySymbol(mod_symbol);
		// Use this function at risk.
		// Modify the position of modification.
		for(std::vector<Reaction>::iterator iter = mod.getModificationReactions().begin(); iter != mod.getModificationReactions().end(); iter++)
		{
			iter->position = index;
		}

		this->modifyFunctionalGroup(fg, mod);
	}

	std::vector<size_t> Modifier::getModificationSetByPosition( const std::string& fg_symbol, const std::string& mod_symbol)
	{
		// Identify if the functional group exists.
		if(!fun_table.containFunctionalGroup(fg_symbol))
			return std::vector<size_t>();

		Modification mod = mod_table.getModificationBySymbol(mod_symbol);
		return mod.getModificationSites(fg_symbol);
	}

	std::vector<size_t> Modifier::getModificationSetByPosition( const FunctionalGroup& fg, const std::string& mod_symbol)
	{
		Modification mod = mod_table.getModificationBySymbol(mod_symbol);
		std::vector<size_t> all_sites = mod.getModificationSites(fg._symbol);
		std::vector<size_t> new_sites;
		BOOST_FOREACH(size_t& site, all_sites)
		{
			if(isFunctionalGroupModifiable(fg, mod_symbol, site))
				new_sites.push_back(site);
		}

		return new_sites;
	}

	bool Modifier::isFunctionalGroupModifiable(const FunctionalGroup& fg, Modification& mod )
	{
		BOOST_FOREACH(Reaction& temp_rt, mod.getModificationReactions())
		{
			// apply the corresponding reactions to the core.
			std::pair<AtomOperation, std::string> single_core_operation;
			BOOST_FOREACH(single_core_operation, temp_rt.core_operation)
			{
				if(single_core_operation.first == Removement) {
					if(!fg.attachFunctionalGroup(single_core_operation.second, temp_rt.position))
						return false;
				}
			}

			// apply the corresponding reactions to the sub-functionalgroup.
			boost::tuple<size_t, AtomOperation, std::string> single_operation;
			std::pair<std::string, OperationSet> fg_operation_pair;
			BOOST_FOREACH(fg_operation_pair, temp_rt.sub_fg_operation)
			{				
				if(!fg.containFunctionalGroup(fg_operation_pair.first, temp_rt.position)) {
					continue;
				}

				FunctionalGroupChain chain;
				chain.push_back(std::make_pair(temp_rt.position, fun_table.getFunctionalGroupBySymbol(fg_operation_pair.first)));
				BOOST_FOREACH(single_operation, fg_operation_pair.second)
				{
					FunctionalGroupChain chain_temp(chain);
					chain_temp.push_back(std::make_pair(single_operation.get<0>(), fun_table.getFunctionalGroupBySymbol(single_operation.get<2>())));
					if(single_operation.get<1>() == Removement) {
						//fg.getSubFunctionalGroup(fg_operation_pair.first, temp_rt.position).removeFunctionalGroup(single_operation.get<2>(), single_operation.get<0>());
						if(!fg.containFunctionalGroupByChain(chain_temp))
							return false;
					}
				}
			}
		}

		return true;
	}


	bool Modifier::isFunctionalGroupModifiable(const FunctionalGroup& fg, const std::string& mod_symbol, const size_t index)
	{
		Modification mod = mod_table.getModificationBySymbol(mod_symbol);
		// Use this function at risk.
		// Modify the position of modification.
		for(std::vector<Reaction>::iterator iter = mod.getModificationReactions().begin(); iter != mod.getModificationReactions().end(); iter++)
		{
			iter->position = index;
		}

		return this->isFunctionalGroupModifiable(fg, mod);
	}

	void Modifier::modifyModificationStatus( const std::string& mod_symbol, const ModificationPosition& pos, int status )
	{
		// Set the status of the rest of modification type on the same position.
		int new_status;
		if(status == 0)
			new_status = -1;
		else if(status == 1)
			new_status = 1;
		else if(status == -1)
			new_status = -1;

		// 1. Locate the position of the modification.
		ModificationByPosition& mod_by_pos = mod_pool.get<mod_position>();
		std::pair<ModificationByPosition::iterator, ModificationByPosition::iterator> p = mod_by_pos.equal_range(pos);

		if(p.first == p.second) {
			std::cout << "The specified position is not found for the modification type" << std::endl;
		}

		for(ModificationByPosition::iterator iter = p.first; iter!=p.second; iter++)
		{
			// 2.a Update the information.
			if(iter->mod_symbol == mod_symbol) {
				this->modifyModificationStatus(iter, status); // Do not change the modification type.
			} else {
				this->modifyModificationStatus(iter, new_status);
			}
		}
	}

	void Modifier::modifyModificationStatus( ModificationByPosition::iterator pos_iter, int new_status )
	{
		ModificationByPosition& mod_by_pos = mod_pool.get<mod_position>();
		ModificationInfo mod_info = *pos_iter;
		mod_info.mod_status = new_status;
		mod_by_pos.replace(pos_iter, mod_info);
	}

	size_t Modifier::getModificationSiteNum( const std::string& mod_symbol, int flag) const
	{
		std::pair<ModificationByCompositeKey::const_iterator, ModificationByCompositeKey::const_iterator> p = mod_pool.get<mod_compo_key>().equal_range(boost::make_tuple(mod_symbol, flag));

		return std::distance(p.first, p.second);
	}

	size_t Modifier::getModificationSiteNum( const std::string& mod_symbol) const
	{
		std::pair<ModificationByCompositeKey::const_iterator, ModificationByCompositeKey::const_iterator> p = mod_pool.get<mod_compo_key>().equal_range(boost::make_tuple(mod_symbol));

		return std::distance(p.first, p.second);
	}

	std::pair<ModificationByPosition::iterator, ModificationByPosition::iterator> Modifier::getModificationSiteIDRange( const MacroPosition& mac_pos )
	{
		ModificationByPosition& mod_pos_index = mod_pool.get<mod_position>();

		ModificationPosition mod_pos_0(mac_pos.branch_id, mac_pos.mono_id, 0);
		ModificationPosition mod_pos_1(mac_pos.branch_id, mac_pos.mono_id+1, 0);

		ModificationByPosition::iterator it0 = mod_pos_index.lower_bound(mod_pos_0);
		ModificationByPosition::iterator it1 = mod_pos_index.upper_bound(mod_pos_1);

		return std::make_pair(it0, it1);
	}

	ModificationSites Modifier::getModificationSitesBySymbol( const std::string& mod_symbol, int flag )
	{
		std::pair<ModificationByCompositeKey::const_iterator, ModificationByCompositeKey::const_iterator> p = mod_pool.get<mod_compo_key>().equal_range(boost::make_tuple(mod_symbol, flag));

		std::set<ModificationPosition> mod_set;
		while(p.first != p.second)
		{
			mod_set.insert(p.first->mod_pos);
			++p.first;
		}

		return mod_set;
	}

	ModificationSites Modifier::getModificationSitesBySymbol( const std::string& mod_symbol)
	{
		std::pair<ModificationByCompositeKey::const_iterator, ModificationByCompositeKey::const_iterator> p = mod_pool.get<mod_compo_key>().equal_range(boost::make_tuple(mod_symbol));

		std::set<ModificationPosition> mod_set;
		while(p.first != p.second)
		{
			mod_set.insert(p.first->mod_pos);
			++p.first;
		}

		return mod_set;
	}

	int Modifier::getModificationStatusByPosition( const ModificationPosition& mod_pos, const std::string& mod_symbol )
	{
		std::pair<ModificationByPosition::const_iterator, ModificationByPosition::const_iterator> p = mod_pool.get<mod_position>().equal_range(mod_pos);

		int status = -100;
		while(p.first != p.second)
		{
			if(p.first->mod_symbol == mod_symbol) {
				status = p.first->mod_status;
				break;
			}
			++p.first;
		}

		return status;
	}

}