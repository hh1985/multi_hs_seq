#include <GAGPL/GAGLIBRARY/SequenceSpace.h>
#include <GAGPL/CHEMISTRY/Modification.h>
#include <GAGPL/MATH/combination.h>
#include <iostream>
#include <algorithm>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace gag
{
	
	SequenceList SequenceSpace::generateModificationCombination()
	{
		SequenceList seq_list;

		// 0. Sequentially list the modification symbol.
		set<string> mod_sym_set = gs->getModificationTypes();

		set<string>::iterator iter = mod_sym_set.begin();

		// 0.1 Get available modification sites.
		
		this->nextModificationSymbol(mod_sym_set, iter, seq_list, ModificationSequence(), gs->getModificationMap());

		return seq_list;
	}

	void SequenceSpace::nextModificationSymbol(set<string>& mod_set, set<string>::iterator iter, SequenceList& seq_list, ModificationSequence mod_seq, Modifier mod_manager)
	{
		// Get currently available modification sites.
		ModificationSites current_sites = mod_manager.getModificationSitesBySymbol(*iter, 1);
		size_t mod_num = (size_t)gs->getModificationConstraint(*iter);

		//ModificationSites mod_sites; 
    size_t i = 0;

    // Non-const version of the container.
    std::vector<ModificationPosition> current_vec(current_sites.begin(), current_sites.end());
    std::vector<ModificationPosition> mod_vec;

		for(ModificationSites::iterator mod_iter = current_sites.begin(); mod_iter != current_sites.end(); mod_iter++)
		{
			mod_vec.push_back(*mod_iter);
			if(++i >= mod_num)
				break;
		}

    //ModificationSites::iterator CurIter(current_sites.begin());
    
		do 
		{
			// Set the status of mod_sites to "occupied".
			Modifier new_manager(mod_manager);
			for(auto mod_iter = mod_vec.begin(); mod_iter != mod_vec.end(); mod_iter++)
			{
				new_manager.modifyModificationStatus(*iter, *mod_iter, 0);
			}
			// Select the current combination, set the status and move to the next symbol.
			ModificationSequence current_seq(mod_seq);
      ModificationSites mod_sites(mod_vec.begin(), mod_vec.end());
			current_seq.insert(std::make_pair(*iter, mod_sites));
			
			set<std::string>::iterator it_copy = iter;
			it_copy++;
			if(it_copy == mod_set.end()) {
				seq_list.push_back(current_seq);
				continue;
			}
			this->nextModificationSymbol(mod_set, it_copy, seq_list, current_seq, new_manager);
		} while(stdcomb::next_combination(current_vec.begin(), current_vec.end(), mod_vec.begin(), mod_vec.end()));

	}

  bool SequenceSpace::nextModificationSymbol( set<string>& mod_set, std::set<string>::iterator mod_iter, ModificationSites prev_occupied_sites )
  {
    if(mod_iter == mod_set.end()) return false;

    // Generate available candidate sites.
    std::vector<ModificationPosition>& candidate_sites = container_seq.at(*mod_iter);
    // Current occupied sites.
    std::vector<ModificationPosition>& occupied_sites = prototype_seq.at(*mod_iter);

    // Not end.
    bool status = nextModificationSymbol(mod_set, std::next(mod_iter, 1), prev_occupied_sites);
    if(!status) { // Can not move deeper.
      if(stdcomb::next_combination(candidate_sites.begin(), candidate_sites.end(),occupied_sites.begin(), occupied_sites.end()))
      {
        // Update prev_occupied_sites.
        std::copy(occupied_sites.begin(), occupied_sites.end(), std::inserter(prev_occupied_sites, prev_occupied_sites.end()));

        // New sequence. Update container_seq and prototype_seq.
        this->updateSites(mod_set, std::next(mod_iter, 1), prev_occupied_sites);
          
        // Update occupied sites.
        //std::copy(occupied_sites.begin(), occupied_sites.end(), std::inserter(prev_occupied_sites, prev_occupied_sites.end()));

        //bool status = nextModificationSymbol(mod_set, std::next(mod_iter, 1), prev_occupied_sites);
        return true;
      } else { // Can not move further.
        // Return to last level.
        return false;
      }
    } else {
      return true;
    }
  }

	std::string SequenceSpace::getSequenceString( const ModificationSequence& mod_seq )
	{
		SequenceCode seq_code = getSequenceCode(mod_seq);
		
		std::string seq_str;
		SequenceCode::iterator iter = seq_code.begin();
		for(; iter!=seq_code.end(); iter++)
		{
			seq_str.append(iter->first);
			seq_str.append(boost::lexical_cast<std::string>(iter->second));
		}
		return seq_str;

	}

	SequenceCode SequenceSpace::getSequenceCode( const ModificationSequence& mod_seq )
	{
		SequenceCode seq_code;

		std::string str;

		Branch& bc = gs->getBranchByID(0);

		ModificationSequence::const_iterator iter = mod_seq.find("SO3");
		const ModificationSites& s_sites = iter->second;

		for(size_t i = 0; i < bc.getUnitNum(); i++)
		{
			std::string symbol = bc.getUnitByID(i).getSymbol();
			
			if(symbol == "DeltaGlcA") {
				str = "D";
			} else if(symbol == "GlcA") {
				str = "G";
			} else if(symbol == "GlcN") { // A, S or H.
				// Check the 2-N position.
				ModificationPosition mod_pos(0, i, 2);
				ModificationSequence::const_iterator ac_iter = mod_seq.find("Ac");
				ModificationSites::const_iterator s_pos_iter = s_sites.find(mod_pos);

				if(ac_iter != mod_seq.end()) {

					const ModificationSites& ac_sites = ac_iter->second;

					ModificationSites::const_iterator ac_pos_iter = ac_sites.find(mod_pos);

					if(ac_pos_iter != ac_sites.end())
						str = "A";
					else if(s_pos_iter != s_sites.end())
						str = "S";
					else {
						str = "H";
					}
				} else {
					if(s_pos_iter != s_sites.end())
						str = "S";
					else {
						str = "H";
					}
				}
	
			} else if(symbol == "GlcNAc") {
				str = "A";
			} else {
				str = "X"; // Unknown type.
			}
			if(iter == mod_seq.end()) {
				seq_code.push_back(std::make_pair(str,0));
				continue;
			}

			size_t num = 0;

			// Get the number of SO3 on the ring.
			if(str == "D" || str == "G") {
				ModificationSites::const_iterator pos_iter = s_sites.find(ModificationPosition(0, i, 2));
				if(pos_iter != s_sites.end()) {
					num = 2;
				}
			} else if(str == "H" || str == "S" || str == "A") {
				for(size_t j = 3; j < bc.getUnitByID(i).getInternalSites().size(); j++)
				{
					ModificationSites::const_iterator pos_iter = s_sites.find(ModificationPosition(0, i, j));
					if(pos_iter != s_sites.end())
						num += pos_iter->site_id;
				}
			} else {
				// Do nothing.
			}

			seq_code.push_back(std::make_pair(str,num));

		}

		return seq_code;
	}

	MergedCode SequenceSpace::getMergedCode( const ModificationSequence& mod_seq )
	{
		MergedCode merged_code;

		std::string str;

		Branch& bc = gs->getBranchByID(0);

		ModificationSequence::const_iterator iter = mod_seq.find("SO3");
		const ModificationSites& s_sites = iter->second;

		for(size_t i = 0; i < bc.getUnitNum(); i++)
		{
			std::string symbol = bc.getUnitByID(i).getSymbol();

			if(symbol == "DeltaGlcA") {
				str = "D";
			} else if(symbol == "GlcA") {
				str = "G";
			} else if(symbol == "GlcN") { // A, S or H.
				// Check the 2-N position.
				ModificationPosition mod_pos(0, i, 2);
				ModificationSequence::const_iterator ac_iter = mod_seq.find("Ac");
				ModificationSites::const_iterator s_pos_iter = s_sites.find(mod_pos);

				if(ac_iter != mod_seq.end()) {

					const ModificationSites& ac_sites = ac_iter->second;

					ModificationSites::const_iterator ac_pos_iter = ac_sites.find(mod_pos);

					if(ac_pos_iter != ac_sites.end())
						str = "A";
					/*else if(s_pos_iter != s_sites.end())
						str = "S";*/
					else {
						str = "H";
					}
				} else {
					//if(s_pos_iter != s_sites.end())
					//	str = "S";
					//else {
					//	str = "H";
					//}
					str = "H";
				}

			} else if(symbol == "GlcNAc") {
				str = "A";
			} else {
				str = "X"; // Unknown type.
			}
			if(iter == mod_seq.end()) {
				merged_code.push_back(std::make_pair(str,0));
				continue;
			}

			size_t num = 0;

			// Get the number of SO3 on the ring.
			if(str == "D" || str == "G") {
				ModificationSites::const_iterator pos_iter = s_sites.find(ModificationPosition(0, i, 2));
				if(pos_iter != s_sites.end()) {
					num = 1;
				}
			} else if(str == "H" || str == "A") {
				for(size_t j = 1; j < bc.getUnitByID(i).getInternalSites().size(); j++)
				{
					ModificationSites::const_iterator pos_iter = s_sites.find(ModificationPosition(0, i, j));
					if(pos_iter != s_sites.end())
						num += 1;
				}
			} else {
				// Do nothing.
			}

			merged_code.push_back(std::make_pair(str,num));

		}

		return merged_code;
	}

	std::string SequenceSpace::getMergedCodeString( const ModificationSequence& mod_seq )
	{
		SequenceCode seq_code = getMergedCode(mod_seq);

		std::string seq_str;
		SequenceCode::iterator iter = seq_code.begin();
		for(; iter!=seq_code.end(); iter++)
		{
			seq_str.append(iter->first);
			seq_str.append(boost::lexical_cast<std::string>(iter->second));
		}
		return seq_str;
	}

  ModificationSequence SequenceSpace::nextModificationSequence()
  {
    // For the initial situation.
    if(container_seq.size() == 0 && prototype_seq.size() == 0) {
      this->initializeModificationSequence();
      return this->getCurrentSequence();
    } else {
      std::set<std::string> mod_set = gs->getModificationTypes();

      bool flag = this->nextModificationSymbol(mod_set, mod_set.begin(), ModificationSites());

      if(flag) // Get new sequence.
        return this->getCurrentSequence();
      else // The end of enumeration.
        return ModificationSequence();
    }
  }

  void SequenceSpace::updateSites( set<string>& mod_set, std::set<string>::iterator mod_iter, ModificationSites prev_occupied_sites )
  {
    if(mod_iter == mod_set.end()) return;

    // Total candidate sites.
    ModificationSites total_sites = gs->getModificationSitesBySymbol(*mod_iter, 1);
    ModificationSites available_sites = getSiteDifference(total_sites, prev_occupied_sites);

    std::vector<ModificationPosition> available_vec(available_sites.begin(), available_sites.end());
    container_seq[*mod_iter] = available_vec;

    int mod_num = gs->getModificationConstraint(*mod_iter);
    
    // Select the first n sites for the prototype_seq.
    std::vector<ModificationPosition> initial_vec(available_vec.begin(), std::next(available_vec.begin(), std::min((size_t)mod_num, available_vec.size())));

#ifdef _DEBUG
  if(((size_t)mod_num > available_vec.size()))
  {
    std::cout << *mod_iter << " Not correct!\n";
    std::cout << "Mod number: " << mod_num;
    printModificationSites(available_sites);
  }
#endif // _DEBUG

    prototype_seq[*mod_iter] = initial_vec;

    // Add the current choice into pre_occupied_sites.
    std::copy(initial_vec.begin(), initial_vec.end(), std::inserter(prev_occupied_sites, prev_occupied_sites.end()));
    // Go to a deeper level.
    this->updateSites(mod_set, std::next(mod_iter, 1), prev_occupied_sites);

    return;
  }

  void SequenceSpace::initializeModificationSequence()
  {
    std::set<std::string> mod_set = gs->getModificationTypes();
    this->updateSites(mod_set, mod_set.begin(), ModificationSites());
  }

  ModificationSequence SequenceSpace::getCurrentSequence() const
  {
    ModificationSequence mod_seq;
    for(auto iter = prototype_seq.begin(); iter != prototype_seq.end(); iter++)
    {
      ModificationSites mod_sites(iter->second.begin(), iter->second.end());
      mod_seq.insert(std::make_pair(iter->first, mod_sites));
    }
    return mod_seq;
  }

  //bool SequenceSpace::isQualifiedSequence(const ModificationSequence& mod_seq)
  //{
  //  auto sulfate_iter = mod_seq.find("SO3");
  //  auto ac_iter = mod_seq.find("Ac");

  //  bool decision = true;

  //  // No sulfate.
  //  if(sulfate_iter == mod_seq.end()) return decision;

  //  // MacroPosition =>
  //  //    set of site_id
  //  std::map<MacroPosition, std::map<size_t, ModificationDistribution::iterator>> residue_sulfate;

  //  // Iterate over all sulfated sites, storing the information into residue_sulfate.
  //  for(auto site_iter = sulfate_iter->second.begin(); site_iter != sulfate_iter->second.end(); site_iter++)
  //  {
  //    residue_sulfate[site_iter->macro_pos].insert(std::make_pair(site_iter->site_id, site_iter));
  //  }

  //  
  //  for(auto residue_iter = residue_sulfate.begin(); residue_iter != residue_sulfate.end(); residue_iter++)
  //  {
  //    // Check the type of the residue.
  //    std::string residue_type = gs->getBranchByID(residue_iter->first.branch_id).getUnitByID(residue_iter->first.mono_id).getSymbol();

  //    if(residue_type == "DeltaGlcA" || residue_type == "GlcA") continue;

  //    if(residue_type == "GlcN" || residue_type == "GlcNAc") {
  //      auto third_iter = residue_iter->second.find(3);
  //      auto sixth_iter = residue_iter->second.find(6);
  //      auto second_iter = residue_iter->second.find(2);
  //    } else  {
  //      cout << residue_type << "\n";
  //      throw std::runtime_error("Undefined residue type!");
  //    }
  //  }
  //  
  //  return decision;
  //}


}