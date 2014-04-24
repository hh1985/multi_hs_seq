/********************************************************************
	created:	2013/04/24
	created:	24:4:2013   10:27
	filename: 	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\src\GAGPL\GAGLIBRARY\SequencePrediction.cpp
	file path:	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\src\GAGPL\GAGLIBRARY
	file base:	SequencePrediction
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include <GAGPL/GAGLIBRARY/SequencePrediction.h>
#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <GAGPL/MATH/MassConversion.h>
#include <algorithm>

namespace gag
{
	ModificationSites SequencePrediction::getComplementaryModificationSites( const ModificationSites& ms)
	{
		ModificationSites complete_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);
		ModificationSites diff_sites;

		std::set_difference(complete_sites.begin(), complete_sites.end(), ms.begin(), ms.end(), std::inserter(diff_sites, diff_sites.end()));

		return diff_sites;		
	}
	void SequencePrediction::build(std::multimap<MonoPeakPtr, NodeItem>& matched_results )
	{
    //std::vector<std::string> type_vec;
    type_vec.push_back("G");
    type_vec.push_back("C");
    type_vec.push_back("GG");
    type_vec.push_back("GC");
    type_vec.push_back("CG");
    type_vec.push_back("CC");

    // Setting the parameter map.    
    //std::map<std::string, double> type_score;
    type_score.insert(std::make_pair("G", param.getParameter<double>("gb_cleavage").first));
    type_score.insert(std::make_pair("C", param.getParameter<double>("cr_cleavage").first));
    type_score.insert(std::make_pair("CG",param.getParameter<double>("int_cleavage").first));
    type_score.insert(std::make_pair("GC", param.getParameter<double>("int_cleavage").first));
    type_score.insert(std::make_pair("GG", param.getParameter<double>("int_cleavage").first));
    type_score.insert(std::make_pair("CC", param.getParameter<double>("int_cleavage").first));

    // 0. Initialization of modification distribution.
    //this->initialize();

    // 1. Group assignments and remove redundant information.
    //AssignmentElite assignment_rep = this->groupAssignments(matched_results);
       
    // 2. Estimate Ac distribution. Notice that the distribution of Ac has 
    // to be estimated before the distribution of SO3.
    this->generateModificationDistribution("Ac", /*assignment_rep, */ matched_results, false);
    
    // 2.5 Clean the incorrect assignment based on Ac distribution.
    this->filterAssignments("Ac", matched_results);
    //this->filterAssignments("Ac", assignment_rep);

    // 3. Estimate SO3 distribution.
    this->generateModificationDistribution("SO3", /*assignment_rep , */ matched_results, false);

	}

  bool SequencePrediction::isNRECleavage( const ModificationSites& ms) const
  {
		ModificationSites complete_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);

		ModificationSites head_sites;
		
		// 1. Get the modification on the first ring.
		MacroPosition mac_pos(0,0);
		for(ModificationSites::iterator iter = complete_sites.begin(); iter != complete_sites.end(); iter++)
		{
			if(mac_pos == MacroPosition(0,0)) {
				mac_pos = iter->macro_pos;
				head_sites.insert(*iter);
			} else if(iter->macro_pos == mac_pos){
				head_sites.insert(*iter);
			} else {
				break;
			}
		}
		
		// 2. Decide if the current modification sites include the first one.
		return std::includes(ms.begin(), ms.end(), head_sites.begin(), head_sites.end());

	}

  bool SequencePrediction::isRECleavage( const ModificationSites& ms) const
  {
		ModificationSites complete_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);

		ModificationSites tail_sites;

		// 1. Get the modification on the last ring.
		MacroPosition mac_pos(0,0);
		for(ModificationSites::reverse_iterator iter = complete_sites.rbegin(); iter != complete_sites.rend(); iter++)
		{
			if(mac_pos == MacroPosition(0,0)) {
				mac_pos = iter->macro_pos;
				tail_sites.insert(*iter);
			} else if(iter->macro_pos == mac_pos){
				tail_sites.insert(*iter);
			} else {
				break;
			}
		}

		// 2. Decide if the current modification sites include the first one.
		return std::includes(ms.begin(), ms.end(), tail_sites.begin(), tail_sites.end());
	}

  void SequencePrediction::printModificationDistribution( std::string mod_symbol )
  {
		// Print title.
		std::cout << "Position\tModification\tIntensity" << std::endl;

		for(ModificationDistribution::const_iterator iter = mod_pattern[mod_symbol].begin(); iter != mod_pattern[mod_symbol].end(); iter++)
		{
			std::cout << iter->first.printString() << "\t" << mod_symbol << "\t" << iter->second << std::endl;
		}
	}

  void SequencePrediction::printModificationDistribution( const ModificationDistribution& dist) const
  {
    // Print title.
    std::cout << "Position\tModification\tIntensity" << std::endl;

    for(ModificationDistribution::const_iterator iter = dist.begin(); iter != dist.end(); iter++)
    {
      std::cout << iter->first.printString() << "\t" << mod_symbol << "\t" << iter->second << std::endl;
    }
  }

  AssignmentGroup SequencePrediction::groupAssignments(const AssignmentSpace& matched_results, bool loss /* = TRUE */)
  {
    std::cout << "Grouping assignments based on " << mod_symbol << "\n";

    AssignmentGroup assignment_rep;

    /* Step 1: organize different groups by the modification sites and number */
    // Key: the pair of modification sites and mod number.
    // Value: the rep structure information.
    // For sulfation, the number of SO3 cannot be used for grouping, simply assign the number as 0.
    typedef std::map<std::pair<ModificationSites, int>, std::map<MonoPeakPtr, NodeItem>> ClusterAssignmentMap;
    // Separate glycosidic-bond cleavage and cross-ring cleavage.
    //ClusterAssignmentMap gb_frag_register;
    //ClusterAssignmentMap cr_frag_register;
    std::vector<ClusterAssignmentMap> frag_cluster(2, ClusterAssignmentMap());
    ClusterAssignmentMap* frag_reg_ptr = nullptr;

    for(AssignmentSpace::const_iterator assign_iter = matched_results.begin();assign_iter!= matched_results.end(); assign_iter++) {
      const FragmentPtr temp_frag = assign_iter->second.getFragment();
           
      // Remove unqualified assignments.
      if(!this->isQualifiedAssignment(*temp_frag)) continue;
         
      if(temp_frag->isGlycosidicBondCleavage()) {
        frag_reg_ptr = &frag_cluster[0];
      } else if(temp_frag->isCrossringCleavage()) {
        frag_reg_ptr = &frag_cluster[1];
      } else {
        ;
      }

      ModificationSites frag_sites = temp_frag->getModificationSitesBySymbol(mod_symbol, 1);

      int mod_num = assign_iter->second.getModificationNum(mod_symbol);

      int group_num = (loss ? -1 : mod_num);
      std::pair<ModificationSites, int> cluster_key = std::make_pair(frag_sites, group_num);

      temp_frag->printFragment();
      printModificationSites(frag_sites);
      std::cout << frag_sites.size() << "\t" << mod_num << "\n";

      //AssignmentMember member;
      //member.mono_pk_ptr = assign_iter->first;
      //member.node = assign_iter->second;
      ClusterAssignmentMap::iterator frag_iter = (*frag_reg_ptr).find(cluster_key);
      if(frag_iter == (*frag_reg_ptr).end()) {
        std::cout << "Not found!\n";
      } else {
        std::cout << "Found!\n";
      }
      (*frag_reg_ptr)[cluster_key].insert(*assign_iter);
      std::cout << "Updated size: " << (*frag_reg_ptr)[cluster_key].size() << "\n";
    }

    /* Step 2: Take the group information */
    for(size_t i = 0; i < 2; i++) {
      for(ClusterAssignmentMap::iterator cluster_iter = frag_cluster[i].begin(); cluster_iter != frag_cluster[i].end(); cluster_iter++) {
        
        std::map<MonoPeakPtr, NodeItem>& cluster_member = cluster_iter->second;
        
        std::map<MonoPeakPtr, NodeItem>::iterator best_iter = cluster_member.begin();
        for(std::map<MonoPeakPtr, NodeItem>::iterator int_iter = cluster_member.begin(); int_iter != cluster_member.end(); int_iter++) {

          if(loss) {
            // Get the member with maximum mod number.
            if(int_iter->second.getModificationNum(mod_symbol) > best_iter->second.getModificationNum(mod_symbol)) {
              best_iter = int_iter;
            } else if(int_iter->second.getModificationNum(mod_symbol) == best_iter->second.getModificationNum(mod_symbol)) {
              if(int_iter->second.getCompositionShiftList().size() < best_iter->second.getCompositionShiftList().size()) {
                best_iter = int_iter;
              }
            }
            
            
          } 
        }
        double weight = (loss ? 1.0 : (double)cluster_member.size());

        // Do something with assginment_rep.
        assignment_rep.insert(AssignmentMember(best_iter->first, best_iter->second, weight));
        
        
      }
    }

    std::cout << "\n";	

    return assignment_rep;
  }

	void SequencePrediction::estimateUniquenessValue( CleavageElite& clv_elite, CleavageSpace& clv_pool )
	{
    std::cout << "\nEstimating uniqueness value:\n"; 
		// 1. Collecting all unique mass values from elite.
    CleavageByExperimentalMass& clv_by_mass = clv_elite.get<exp_mass>();
		for(CleavageByExperimentalMass::iterator mass_iter = clv_by_mass.begin(); mass_iter != clv_by_mass.end(); mass_iter++)
		{
      std::cout << "\n";
		// 2. Checking the abstract assignments in clv_pool.
			std::set<CleavageItem> temp_gb_vec; // Glycosidic-bond.
			std::set<CleavageItem> temp_cr_vec; // Cross-ring.
			std::set<CleavageItem> temp_int_vec; // Internal.
	
			// Search the alternative explanations in the pool.
			std::pair<CleavageByExperimentalMass::iterator, CleavageByExperimentalMass::iterator> mass_pair = clv_pool.get<exp_mass>().equal_range(mass_iter->experimental_mass); 

			while(mass_pair.first != mass_pair.second)
			{
        //std::cout << mass_pair.first->experimental_mass << "\t" << mass_pair.first->general_type << "\n";
        if(mass_pair.first->isGlycosidicBondCleavage()) {
          temp_gb_vec.insert(*(mass_pair.first));
        } else if(mass_pair.first->isCrossringCleavage()) { 
          temp_cr_vec.insert(*(mass_pair.first));
        } else if(mass_pair.first->isInternalCleavage()) {
          temp_int_vec.insert(*(mass_pair.first));
        } else {
          throw std::runtime_error("Undefined cleavage type!");
        }

				mass_pair.first++;
			}

			// 3. Calculate the uniqueness value of terminal assignment.
			std::cout << "Mass: " << mass_iter->experimental_mass << " GB assignment(s): " << temp_gb_vec.size() << " CR assignment(s): " << temp_cr_vec.size() << " Internal cleavage: " << temp_int_vec.size() << "\n";

			// At least one assignment for terminal glycosidic-bond fragment is required. Should be more intelligent.
			if(temp_gb_vec.size() == 0 && temp_cr_vec.size() == 0) continue;

			// Storing the modification sites for reducing assignments.
      // Consider treating the mod number difference as punishment.
			std::set<std::pair<ModificationSites, int>> total_mod_set;

			double total_score = 0.0;
      double weight0 = param.getParameter<double>("gb_cleavage").first;
			for(std::set<CleavageItem>::iterator iter = temp_gb_vec.begin(); iter != temp_gb_vec.end();  iter++) {
				total_mod_set.insert(std::make_pair(iter->mod_sites, iter->mod_num));
				total_score += 1.0;
			}

			// For the rest of the assignments, check and see if they can be reduced to the total_mod_set.
			double weight1 = param.getParameter<double>("cr_cleavage").first; // Weight for cross-ring cleavage.
			for(std::set<CleavageItem>::iterator iter = temp_cr_vec.begin(); iter != temp_cr_vec.end();  iter++) {
				std::set<std::pair<ModificationSites, int>>::iterator mod_iter = total_mod_set.find(std::make_pair(iter->mod_sites, iter->mod_num));
				// If found, do nothing. Otherwise, add to the total score.
				if(mod_iter == total_mod_set.end()) {
					total_mod_set.insert(std::make_pair(iter->mod_sites, iter->mod_num));
					total_score += weight1;
				} else {
          std::cout << "Redundant and reduced!\n";
        }
			}

      double weight2 = param.getParameter<double>("int_cleavage").first; // Weight for internal cleavage.
			for(std::set<CleavageItem>::iterator iter = temp_int_vec.begin(); iter != temp_int_vec.end();  iter++) {
				std::set<std::pair<ModificationSites, int>>::iterator mod_iter = total_mod_set.find(std::make_pair(iter->mod_sites, iter->mod_num));
				if(mod_iter == total_mod_set.end()) {
					total_mod_set.insert(std::make_pair(iter->mod_sites, iter->mod_num));
					total_score += weight2;
				} else {
          std::cout << "Redundant and reduced!\n";
        }
			}

			// Updating the score for the current abstract fragment.
      // The replacing might be costly.
      CleavageItem clv_item = *mass_iter;
      if(clv_item.isGlycosidicBondCleavage())
        clv_item.setUniquenessValue(weight0/total_score);
      else if(clv_item.isCrossringCleavage())
        clv_item.setUniquenessValue(weight1/total_score);
      else
        throw std::runtime_error("Strange cleavage type.");

      //clv_item.estimated_mod_num = clv_item.mod_num * clv_item.uniqueness_confidence;
      clv_by_mass.replace(mass_iter, clv_item);
			//mass_iter->uniqueness_confidence = 1.0/total_score;
			std::cout << "The total score is: " << clv_item.uniqueness_confidence << "\n";
		}

		
	}

	void SequencePrediction::filterAssignments( const std::string& mod_symbol, AssignmentSpace& assignment_pool)
	{
		const ModificationDistribution& mod_dist = mod_pattern[mod_symbol];
		int total_number = gs->getModificationConstraint(mod_symbol);

    // 1. Sort modification sites by likelihood values. Take the top ones.
    std::multimap<double, ModificationPosition> score_map;
    for(ModificationDistribution::const_iterator iter = mod_dist.begin(); iter != mod_dist.end(); iter++)
      score_map.insert(std::make_pair(iter->second, iter->first));

    std::multimap<double, ModificationPosition>::reverse_iterator rev_iter = score_map.rbegin();
    ModificationSites mod_sites;
    for(; rev_iter != score_map.rend(); rev_iter++)
    {
        total_number--;
        if(total_number < 0) break;

        // Set the corresponding position as unavailable.
        mod_sites.insert(rev_iter->second);
    }

    // 2. Keep only the data qualifying for the modification.
    AssignmentSpace::iterator assign_iter = assignment_pool.begin();
    while(assign_iter != assignment_pool.end())
    //for(AssignmentSpace::iterator iter = assignment_pool.begin(); iter != assignment_pool.end(); iter++)
    {
        // 2.1 For each node item, check if it meets the requirement of modification status.
        ModificationSites available_sites = assign_iter->second.getFragment()->getModificationSitesBySymbol(mod_symbol, 1);

        // Decide the supposed number of modifications.
        ModificationSites modified_sites;
        std::set_intersection(available_sites.begin(), available_sites.end(), mod_sites.begin(), mod_sites.end(), std::inserter(modified_sites, modified_sites.end()));

        std::cout << "MZ: " << assign_iter->first->mz << "\tZ: " << assign_iter->first->z << "\tMass: " << msmath::calculateMass(assign_iter->first->mz, assign_iter->first->z) << "\tType: " << assign_iter->second.getFragment()->getCleavageType() << "\tSite size: " << modified_sites.size() << "\tAc: " << assign_iter->second.getModificationNum("Ac") << "\tSupposed number: " << modified_sites.size() << "\tSO3: " << assign_iter->second.getModificationNum("SO3") << "\n";

        // If the deduced occupied site number != specified site number: unreasonable assignment.
        if((int)modified_sites.size() != assign_iter->second.getModificationNum(mod_symbol)) {
            std::cout << "\tDeleted";
            // Notice the usage of post-increment!!!
            assignment_pool.erase(assign_iter++);
        } else {
            // Pass!
            ++assign_iter;
        }

        std::cout << "\n";
    }

    // 3. Update the modification status.
    for(ModificationSites::const_iterator iter = mod_sites.begin(); iter != mod_sites.end(); iter++)
        gs->modifyModificationStatus(mod_symbol, *iter);

    // Update the status for all the rest of the fragments. Needs improvement in the future.
    std::set<FragmentPtr> frag_set;
    for(AssignmentSpace::const_iterator iter = assignment_pool.begin(); iter != assignment_pool.end(); iter++)
    {
        frag_set.insert(iter->second.getFragment());
    }

    for(std::set<FragmentPtr>::iterator frag_iter = frag_set.begin(); frag_iter != frag_set.end(); frag_iter++)
        for(ModificationSites::const_iterator pos_iter = mod_sites.begin(); pos_iter != mod_sites.end(); pos_iter++)
            (*frag_iter)->modifyModificationStatus(mod_symbol, *pos_iter);
	
  }

  bool SequencePrediction::generateModificationDistribution( std::string symbol, const AssignmentSpace& assignment_pool, bool separate /* = true */ )
  {
    // 0. Check and see if the mod_type exists.
    mod_symbol = symbol;
    mod_status = "G";
        
    std::cout << "Generate modification distribution -- " << mod_symbol << "\n";
        
    int total_mod_num = gs->getModificationConstraint(mod_symbol);
    if(total_mod_num == 0) return false;

    ModificationSites complete_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);

    bool flag = (mod_symbol == "Ac" ? false : true);
    
    //AssignmentGroup assignment_rep = this->groupAssignments(assignment_pool, flag);
    // 1. Group assignments by fragment type to create independent observations.
    AssignmentSpace assignment_rep = this->groupAssignmentsByFragmentType(assignment_pool, flag);

    // 2. Sort the assignments by their uniqueness values.
    typedef std::multimap<double, AssignmentSpace::iterator> AssignmentCluster;
    AssignmentCluster uni_assignments;
    std::set<double> uni_keys;
    for(AssignmentSpace::iterator assign_iter = assignment_rep.begin(); assign_iter != assignment_rep.end(); assign_iter++) {
      double uni_score = estimateUniquenessValue(assign_iter, assignment_pool);
      uni_assignments.insert(std::make_pair(uni_score, assign_iter));
      uni_keys.insert(uni_score);
    }
    
    CleavageContainer clv_rep, clv_cr_rep;
    for(std::set<double>::reverse_iterator rev_iter = uni_keys.rbegin(); rev_iter != uni_keys.rend(); rev_iter++)
    {
      double key = *rev_iter;
      std::pair<AssignmentCluster::iterator, AssignmentCluster::iterator> p = uni_assignments.equal_range(key);

      std::cout << "Uniqueness score: " << key << "\n";
      AssignmentSpace ideo_uni;
      while(p.first != p.second) {
        
        std::cout << "MZ: " << p.first->second->first->mz << "\tCharge: " << p.first->second->first->z << "\tType: " << p.first->second->second.getCleavageType() << "\n";
        if(!flag && p.first->second->second.getGeneralType() != "G") {
        } else {
          ideo_uni.insert(*(p.first->second));
        }
        
        p.first++;
      }

      std::cout << "Size before grouping: " << ideo_uni.size() << "\n";
      AssignmentSpace clv_prior = this->groupAssignmentsByModificationSites(ideo_uni, flag);

      std::cout << "Size after grouping: " << clv_prior.size() << "\n";

      // 2.1 Convert the assignments into cleavage item.
      for(AssignmentSpace::iterator ideo_iter = clv_prior.begin(); ideo_iter != clv_prior.end(); ideo_iter++)
      {

        CleavageItem clv_item = this->convertNodeItem(*(ideo_iter->first), ideo_iter->second);

        if(clv_item.mod_sites.size() == 0 || clv_item.mod_sites.size() == complete_sites.size()) continue;

        // The missing site number smaller than the missing mod number.
        ModificationSites comp_sites = this->getComplementaryModificationSites(clv_item.mod_sites);
        if((int)comp_sites.size() < (total_mod_num - clv_item.mod_num))
          continue;

        // Update the uniqueness of the cleavage item.
        clv_item.setUniquenessValue(key);

        if(separate && clv_item.isCrossringCleavage())
           clv_cr_rep.insert(clv_item);
        else
          clv_rep.insert(clv_item);
      }
    }

    // Initialize modification distribution.
    this->initialize();

    // The tree storing background information.
    CleavageTree last_tree(gs, mod_symbol);

    // Interaction between the distribution and graph.
    this->updateGlobalDistribution(clv_rep, last_tree);

    //// 1. Convert the assignments to abstract fragments.
    //CleavageContainer clv_pool, clv_rep, clv_cr_rep;

    //// 1.1 Conversion of assignment representative.
    //for(AssignmentGroup::const_iterator rep_iter = assignment_rep.begin(); rep_iter != assignment_rep.end(); rep_iter++) 
    //{
    //  CleavageItem clv_item = this->convertNodeItem(*(rep_iter->mono_pk_ptr), rep_iter->node);
    //      
    //  // Ignore NULL site and full site.
    //  if(clv_item.mod_sites.size() == 0 || clv_item.mod_sites.size() == complete_sites.size())
    //    continue;

    //  // The missing site number smaller than the missing mod number.
    //  ModificationSites comp_sites = this->getComplementaryModificationSites(clv_item.mod_sites);
    //  if((int)comp_sites.size() < (total_mod_num - clv_item.mod_num))
    //    continue;

    //  clv_item.addRepeatNumber((int)rep_iter->weight);

    //  if(separate) {
    //    if(clv_item.isCrossringCleavage()) {
    //      clv_cr_rep.insert(clv_item);
    //    } else {
    //      clv_rep.insert(clv_item);
    //    }
    //  } else {
    //    clv_rep.insert(clv_item);
    //  }
    //  
    //}

    //  // 1.2 Conversion of assignment pool.
    //  for(AssignmentSpace::const_iterator pool_iter = assignment_pool.begin(); pool_iter != assignment_pool.end(); pool_iter++)
    //    clv_pool.insert(this->convertNodeItem(*(pool_iter->first), pool_iter->second));

    //  std::cout << "Elite size: " << clv_rep.size() << "\n";
    //  std::cout << "Crossing Elite size: " << clv_cr_rep.size() << "\n";
    //  std::cout << "Pool size: " << clv_pool.size() << "\n";

    //  // 2. Estimate the uniqueness value.
    //  this->estimateUniquenessValue(clv_rep, clv_pool);

    //  // Initialize modification distribution.
    //  this->initialize();

    //  // The tree storing background information.
    //  CleavageTree last_tree(gs, mod_symbol);

    //  // Interaction between the distribution and graph.
    //  this->updateGlobalDistribution(clv_rep, last_tree);

    //  // 3. Doing the same for crossing-ring cleavage.  Except the modification distribution allocation policy.
      if(separate && mod_symbol == "SO3") {
        std::cout << "Add cross-ring cleavage into consideration\n";

        mod_status = "C";
        //this->estimateUniquenessValue(clv_cr_rep, clv_pool);

        this->updateGlobalDistribution(clv_cr_rep, last_tree);
      }
      
    return true;
  }

  void SequencePrediction::initialize()
  {
    std::cout << "\nInitialize: " << std::endl;

    ModificationDistribution mod_dist;
    //this->initializeModificationDistribution(mod_dist);
    this->equalInitialization(mod_dist);

    mod_pattern[mod_symbol] = mod_dist;
    printModificationDistribution(mod_symbol);

  }

  CleavageItem SequencePrediction::convertNodeItem(const MonoPeak& mono_pk, const NodeItem& node)
  {
    double mass = msmath::calculateMass(mono_pk.mz, -1*mono_pk.z);

    ModificationSites mod_sites = node.getFragment()->getModificationSitesBySymbol(mod_symbol, 1);

    std::cout << "MZ: " << mono_pk.mz << "\tCharge: " << -1 * mono_pk.z << "\tMass: " << mass << std::endl;

    return CleavageItem(mass, mod_sites, node.getModificationNum(mod_symbol), node.getCleavageNum(), node.getGeneralType());

    
  }

  void SequencePrediction::updateLocalDistribution( ModificationDistribution& mod_dist, CleavageItem& clv, CleavageNodePtr child, CleavageNodePtr parent )
  {
    std::cout << "Update local distribution!\n";

    const ModificationSites& child_sites = child->getModificationSites();
    const ModificationSites& parent_sites = parent->getModificationSites();
    const ModificationSites& current_sites = clv.mod_sites;

    // The background number for the cleavage.
    double bg_num = this->getAccumulatedDensity(mod_dist, child_sites, current_sites);

    ModificationSites begin_sites;
    double ch_num = this->getAccumulatedDensity(mod_dist, begin_sites, child_sites);

    double lambda = clv.getConfidenceValue();
      
    // Update the modification number.
    double updated_mod_num = lambda * (clv.mod_num - ch_num) + (1-lambda) * bg_num;
    double total_num = this->getAccumulatedDensity(mod_dist, child_sites, parent_sites);

    // Not negative values. In the new version, there might be a piercing mechanism.
    ModificationSites bg_sites = this->getDiffSites(child_sites, current_sites);
    ModificationSites comp_sites = this->getDiffSites(current_sites, parent_sites);

    if(updated_mod_num >= total_num) {
      updated_mod_num = total_num * (1 - 1e-2);
    } else if(updated_mod_num < 0.0) {
      updated_mod_num = 1e-2; // Prevent negative value.
    } else if(updated_mod_num > (double)bg_sites.size()) {
      updated_mod_num = (double)bg_sites.size(); // Prevent value larger than 1.0.
    }

    std::cout << "Mass: " << clv.experimental_mass << "\t" << clv.getGeneralType() << " ";
    std::cout << mod_symbol << ": " << clv.mod_num << "\t" << updated_mod_num << "\t" << bg_num << "\t" << lambda << "\n";
    printModificationSites(clv.mod_sites);
    std::cout << "Confidence: " << clv.getConfidenceValue() << "\n";
    std::cout << clv.uniqueness_confidence << "\t" << clv.background_confidence << "\t" << clv.intactness_confidence << "\n";

    // Count the number of sites on each residue.
    std::map<size_t, int> residue_map;

    for(ModificationSites::iterator site_iter = bg_sites.begin(); site_iter != bg_sites.end(); site_iter++) {
      size_t mono_id = site_iter->getMonosaccharideID();
      if(residue_map.find(mono_id)== residue_map.end()) {
        residue_map[mono_id] = 1;
      } else {
        residue_map[mono_id] += 1;
      }
    }

    double mono_avg1 = updated_mod_num / (double)residue_map.size();
    
    // Reallocate the mod number.
    //if(mono_avg1 <= 1.0) { // The averaging is on the residue level.
    //  /*std::cout << "Scale the bg number!\n";
    //  double scale_factor1 = updated_mod_num / bg_num;
    //  for(ModificationSites::iterator mod_iter = bg_sites.begin(); mod_iter != bg_sites.end(); mod_iter++)
    //    mod_dist[*mod_iter] *= scale_factor1;

    //  ModificationSites comp_sites = this->getDiffSites(current_sites, parent_sites);
    //  double scale_factor2 = (total_num - updated_mod_num) / (total_num - bg_num);

    //  for(ModificationSites::iterator mod_iter = comp_sites.begin(); mod_iter != comp_sites.end(); mod_iter++)
    //    mod_dist[*mod_iter] *= scale_factor2;*/

    //  std::cout << "Averaging over the monosaccharide residue." << std::endl;
    //  
    //  for(ModificationSites::iterator site_iter = bg_sites.begin(); site_iter != bg_sites.end(); site_iter++) {
    //    size_t mono_id = site_iter->getMonosaccharideID();
    //    mod_dist[*site_iter] = mono_avg1 / (double)residue_map[mono_id];

    //    std::cout << "ID: " << mono_id << "\t" << mod_dist[*site_iter] << std::endl;
    //  }

    //  std::map<size_t, int> residue_map2;

    //  for(ModificationSites::iterator site_iter = comp_sites.begin(); site_iter != comp_sites.end(); site_iter++) {
    //    size_t mono_id = site_iter->getMonosaccharideID();
    //    if(residue_map2.find(mono_id)== residue_map2.end()) {
    //      residue_map2[mono_id] = 1;
    //    } else {
    //      residue_map2[mono_id] += 1;
    //    }
    //  }

    //  double mono_avg2 = (total_num - updated_mod_num) / (double) residue_map2.size();
    //  for(ModificationSites::iterator site_iter = comp_sites.begin(); site_iter != comp_sites.end(); site_iter++) {
    //    size_t mono_id = site_iter->getMonosaccharideID();
    //    mod_dist[*site_iter] = mono_avg2 / (double)residue_map2[mono_id];

    //    std::cout << "ID: " << mono_id << "\t" << mod_dist[*site_iter] << std::endl;
    //  }

    //} else  { // The averaging is on the site level.
      std::cout << "Average mod number over the sites!\n";
      this->averageBySite(mod_dist, bg_sites, updated_mod_num);
      /*for(ModificationSites::iterator mod_iter = bg_sites.begin(); mod_iter != bg_sites.end(); mod_iter++) {
        double avg1 = updated_mod_num / (double)bg_sites.size();
        mod_dist[*mod_iter] = avg1;
        std::cout << "Average 1: " << avg1 << "\n";
      }*/

      this->averageBySite(mod_dist, comp_sites, total_num - updated_mod_num);
      /*for(ModificationSites::iterator mod_iter = comp_sites.begin(); mod_iter != comp_sites.end(); mod_iter++) {
        double avg2 = (total_num - updated_mod_num) / (double)comp_sites.size();
        mod_dist[*mod_iter] = avg2;
        std::cout << "Average 2: " << avg2 << "\n";
      }*/

    //} 

    //std::cout << "Total num: " << total_num << "\tBg num: " << bg_num << "\tScale factor 1: " << scale_factor1 << "\tScale factor 2: " << scale_factor2 << "\n\n";

    printModificationDistribution(mod_symbol);
    std::cout << "End updating local distribution!\n\n";
  }

  void SequencePrediction::updateLocalDistribution( CleavageItem& clv, CleavageTree& clv_tree )
  {
    

    CleavageNodePtr child = *(clv_tree.getSubsets(clv).begin());
    CleavageNodePtr parent = *(clv_tree.getSupersets(clv).begin());

    const ModificationSites& child_sites = child->getModificationSites();
    const ModificationSites& parent_sites = parent->getModificationSites();
    const ModificationSites& current_sites = clv.mod_sites;

    const std::vector<CleavageItem>& child_items = child->getCleavageItems();

    // Find the one with maximum confidence value.
    int child_index = 0; double max_child_conf(0.0);
    for(size_t i = 0; i < child_items.size(); i++)
    {
      double conf = child_items[i].getConfidenceValue();
      if(conf > max_child_conf) {
        child_index = (int)i;
        max_child_conf = conf;
      }
    }
    // incompatible.
    if(clv.mod_num - child_items[child_index].mod_num > (int)clv.mod_sites.size() - (int)child_items[child_index].mod_sites.size()) {
      std::cout << "Incompatible child!\n";
      child = *(clv_tree.getSubsets(child).begin());
    }

    const std::vector<CleavageItem>& parent_items = parent->getCleavageItems();

    int parent_index = 0; double max_parent_conf(0.0);
    for(size_t i = 0 ; i < parent_items.size(); i++)
    {
      double conf = parent_items[i].getConfidenceValue();
      if(conf > max_parent_conf) {
        parent_index = (int)i;
        max_parent_conf = conf;
      }
    }
    if(parent_items[parent_index].mod_num - clv.mod_num > (int)parent_items[parent_index].mod_sites.size() - (int)clv.mod_sites.size())
    {
      std::cout << "Incompatible parent!\n";
      parent = *(clv_tree.getSupersets(parent).begin());
    }

    this->updateLocalDistribution(mod_pattern[mod_symbol], clv, child, parent);
  }

  bool SequencePrediction::isQualifiedAssignment(const Fragment& frag) const
  {
    if(frag.isTerminalCleavage()) {
      return true;
    } else {
      return false;
    }
  }

  void SequencePrediction::scalingModificationDistribution()
  {
    ModificationDistribution& mod_dist = mod_pattern[mod_symbol];

    std::vector<double> density;
    for(ModificationDistribution::iterator iter = mod_dist.begin(); iter != mod_dist.end(); iter++)
    {
      density.push_back(iter->second);
    }

    double min_value = *std::min_element(density.begin(), density.end());
    double max_value = *std::max_element(density.begin(), density.end());
      
    if(min_value == max_value) return;

    double sum = 0.0;
    for(ModificationDistribution::iterator iter = mod_dist.begin(); iter != mod_dist.end(); iter++)
    {
      iter->second = (iter->second - min_value + 1e-6) / (max_value - min_value + 1e-6);
      sum += iter->second;
    }

    int total_num = gs->getModificationConstraint(mod_symbol);
    // Normalization.

    for(ModificationDistribution::iterator iter = mod_dist.begin(); iter != mod_dist.end(); iter++)
    {
      iter->second = total_num * iter->second / sum;
    }
      
  }

  double SequencePrediction::getAccumulatedDensity( const ModificationDistribution& mod_dist, const ModificationSites& child_sites, const ModificationSites& parent_sites )
  {
    ModificationSites diff_sites;
    std::set_difference(parent_sites.begin(), parent_sites.end(), child_sites.begin(), child_sites.end(),  std::inserter(diff_sites, diff_sites.end()));

    double sum = 0.0;
    for(ModificationSites::iterator iter = diff_sites.begin(); iter != diff_sites.end(); iter++)
    {
      ModificationDistribution::const_iterator const_iter = mod_dist.find(*iter);
      if(const_iter != mod_dist.end())
        sum += const_iter->second;
      else
        std::cout << "Site missing!\n";
    }
    return sum;
  }

  ModificationSites SequencePrediction::getDiffSites( const ModificationSites& child_sites, const ModificationSites& parent_sites )
  {
    ModificationSites diff_sites;
    std::set_difference(parent_sites.begin(), parent_sites.end(), child_sites.begin(), child_sites.end(),  std::inserter(diff_sites, diff_sites.end()));

    return diff_sites;
  }

  void SequencePrediction::updateGlobalDistribution( CleavageElite& clv_rep, CleavageTree& last_tree )
  {

    std::cout << "\n****Update global distribution!****\n" << std::endl;

    // 1. Sort assignments by their uniqueness value.
    std::set<double> conf_keys;
    for(CleavageByUniquenessConfidence::iterator iter = clv_rep.get<uniqueness_confidence>().begin(); iter != clv_rep.get<uniqueness_confidence>().end(); iter++)
    {
      if(iter->uniqueness_confidence == 0.0) {
        std::cout << "Uniqueness is zero: Passed!\n";
        continue;
      }
      conf_keys.insert(iter->uniqueness_confidence);
    }
      
    for(std::set<double>::reverse_iterator rev_iter = conf_keys.rbegin(); rev_iter != conf_keys.rend(); rev_iter++)
    {
      std::pair<CleavageByUniquenessConfidence::iterator, CleavageByUniquenessConfidence::iterator> uni_pair = clv_rep.get<uniqueness_confidence>().equal_range(*rev_iter);

      // 1.1 confidence from the background graph and distribution.
      if(rev_iter != conf_keys.rbegin())
        this->estimateBackgroundConfidence(uni_pair, last_tree);

      // 1.2 confidence from the peer nodes.
      this->estimatePeerConfidence(uni_pair);

      // 1.3 Update graph and modification distribution.
      // Assignments are sorted by their confidence values.
      // and sequentially inserted into the graph.
      typedef std::multimap<double, CleavageByUniquenessConfidence::iterator> ConfidenceMap;
      // The map of modification sites, mod number and confidence.
      //typedef std::map<std::pair<ModificationSites, int>, double> SiteMap;
      ConfidenceMap conf_map;
      //SiteMap site_map;
      for(CleavageByUniquenessConfidence::iterator iter = uni_pair.first; iter != uni_pair.second; iter++)
      {
        conf_map.insert(std::make_pair(iter->getConfidenceValue(), iter));
        //site_map.insert(std::make_pair(std::make_pair(iter->mod_sites, iter->mod_num), iter->getConfidenceValue()));
      }

      for(ConfidenceMap::reverse_iterator rev_iter = conf_map.rbegin(); rev_iter != conf_map.rend(); rev_iter++)
      {
        CleavageItem& clv = const_cast<CleavageItem&>(*(rev_iter->second));
        // Check if the cleavage has been included in the background.
        CleavageNodePtr check_ptr = last_tree.locateModificationSites(clv.mod_sites);
        if(check_ptr != nullptr) {
          std::cout << "Recorded! Update confidence value and distribution\n";
          continue;
        }

        // Check if there is contradictory information.
        if(clv.getConfidenceValue()==-1.0)
          continue;

        this->updateLocalDistribution(clv, last_tree);
        
        last_tree.insertNode(clv);

        // Insert complementary node.
        CleavageItem comp_item(clv);
        comp_item.mod_num = gs->getModificationConstraint(mod_symbol) - clv.mod_num;
        comp_item.mod_sites = this->getComplementaryModificationSites(clv.mod_sites);
        last_tree.insertNode(comp_item);
      }

      this->printModificationDistribution(mod_symbol);
      std::cout << "\n";
    }

  }

  void SequencePrediction::estimateBackgroundConfidence( std::pair<CleavageByUniquenessConfidence::iterator, CleavageByUniquenessConfidence::iterator> uni_pair, CleavageTree& clv_tree)
  {
    std::cout << "\nEstimate from background:\n";

    // Iterate over all assignments.
    for(CleavageByUniquenessConfidence::iterator iter = uni_pair.first; iter != uni_pair.second; iter++)
    {
      CleavageItem& clv = const_cast<CleavageItem&>(*iter);
      
      printModificationSites(clv.mod_sites);
      std::cout << "Mass: " << clv.experimental_mass << "\t" << mod_symbol << "\t" << clv.mod_num << "\n";
      // If current cleavage is contradictory to background records, assign it low score.
      
      CleavageNodePtr record = clv_tree.locateModificationSites(clv.mod_sites);
      int total_num = gs->getModificationConstraint(mod_symbol);
      double confidence_constant = param.getParameter<double>("intactness_constant").first;

      // Something found!
      if(record != nullptr) {
        std::cout << "Record found! Pass graph inference.\n";
        const std::vector<CleavageItem>& record_items = record->getCleavageItems();
        // If the information has been recorded by the background, just inherit the confidence. Otherwise, assign the lowest confidence value.
        bool match = false;
        double confidence_score(0.0);
        for(size_t i = 0; i < record_items.size(); i++) {
          if(record_items[i].mod_num == clv.mod_num) {
            confidence_score = record_items[i].getConfidenceValue();         
            match = true;
            break;
          }
        }

        if(match) {
          clv.background_confidence = confidence_score;
        } else {
          clv.background_confidence = pow(confidence_constant, total_num * 10);
        }

      } else {
        double value = this->inferFromGraph(clv, clv_tree);
        
        //int total_num = gs->getModificationConstraint(mod_symbol);
        if(value == -1.0) {
          clv.background_confidence = pow(confidence_constant, total_num * 10);
        } else {
          clv.background_confidence = value;
        }
      }

      
      std::cout << "Update conference: " << clv.background_confidence << "\n\n";
    }

    std::cout << "End of background estimation!\n\n";
  }

  void SequencePrediction::estimatePeerConfidence( std::pair<CleavageByUniquenessConfidence::iterator, CleavageByUniquenessConfidence::iterator> uni_pair )
  {
    std::cout << "\nEstimate from peer graph:\n";

    double confidence_constant = param.getParameter<double>("intactness_constant").first;

    // Check each node in the current tree and complemenatry tree.
    CleavageTree current_tree(gs, mod_symbol);
    CleavageTree comp_tree(gs, mod_symbol);

    for(CleavageByUniquenessConfidence::iterator iter = uni_pair.first; iter != uni_pair.second; iter++)
    {
      // Insert all assignments into the tree.
      CleavageItem& clv = const_cast<CleavageItem&>(*iter);
        
      //this->updateCleavageTree(clv, current_tree);
      current_tree.insertNode(clv);

      // Insert complementary node into complementary tree.
      CleavageItem comp_item(clv);
      comp_item.mod_num = gs->getModificationConstraint(mod_symbol) - clv.mod_num;
      comp_item.mod_sites = this->getComplementaryModificationSites(clv.mod_sites);

      comp_tree.insertNode(comp_item);
        
    }

    std::vector<double> conf_vec;
    for(CleavageByUniquenessConfidence::iterator iter = uni_pair.first; iter != uni_pair.second; iter++)
    {
      // Estimate the confidence from peers.
      CleavageItem& clv = const_cast<CleavageItem&>(*iter);
      
      printModificationSites(clv.mod_sites);
      std::cout << "Mass: " << clv.experimental_mass << "\t" << mod_symbol << "\t" << clv.mod_num << "\tUniqueness: " << clv.uniqueness_confidence << "\n";

      std::cout << "\nInfer from the tree:\n";
      double val1 = this->inferFromGraph(clv, current_tree);
      
      double total_num = gs->getModificationConstraint(mod_symbol);

      // Zero tolerance with -1.0
      if(val1 == -1.0) {
        conf_vec.push_back(-1.0);
      } else {
        std::cout << "\nInfer from the complementary tree.\n";
        double val2 = this->inferFromGraph(clv, comp_tree, false);

        conf_vec.push_back(std::max(val1, val2));
      }
      
    }

    int index = 0;
    for(CleavageByUniquenessConfidence::iterator iter = uni_pair.first; iter != uni_pair.second; iter++)
    {
      CleavageItem& clv = const_cast<CleavageItem&>(*iter);
      // For peer confidence, it is important to append the uniqueness values.
      if(conf_vec[index] == -1.0)
        clv.background_confidence = -1.0;
      else
        // The reason of using sqrt is to make the likelihood distinctive between different candidate sites.
        clv.background_confidence = sqrt(clv.background_confidence * conf_vec[index]);
        //clv.background_confidence = clv.background_confidence * conf_vec[index];

      index++;

      printModificationSites(clv.mod_sites);
      std::cout << mod_symbol << "\t" << clv.mod_num << "\n";

      std::cout << "Updated confidence: " << clv.background_confidence << "\n";
    }

    std::cout << "End of peer estimation!\n\n";
  }

  void SequencePrediction::initializeModificationDistribution( ModificationDistribution& mod_dist )
  {
    ModificationSites complete_sites = gs->getModificationSitesBySymbol(mod_symbol);
    ModificationSites available_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);
    ModificationSites diff_sites = this->getDiffSites(available_sites, complete_sites);

    std::cout << "Complete: " << complete_sites.size() << "\tAvailable: " << available_sites.size() << "\tDifference: " << diff_sites.size() << "\n";

    int total_num = gs->getModificationConstraint(mod_symbol);

    // The distribution is not uniform on sites, but based on residues. However, for highly sulfated HS, a single site might have number larger than 1. Therefore, calculate the distance to max sites, and average them along the residue.
    // The residue ID and the positions on the residue.
    std::map<size_t, int> residue_map;

    for(ModificationSites::iterator site_iter = available_sites.begin(); site_iter != available_sites.end(); site_iter++) {
      size_t mono_id = site_iter->getMonosaccharideID();
      if(residue_map.find(mono_id)== residue_map.end()) {
        residue_map[mono_id] = 1;
      } else {
        residue_map[mono_id] += 1;
      }
    }

    // If the sulfate number is larger than 0.
    double avg0 = (double)total_num / (double)residue_map.size();

    if(avg0 > 1.0) {
      // Highly sulfated!
      double avg = (double)total_num / (double)available_sites.size();
      for(ModificationSites::iterator site_iter = available_sites.begin(); site_iter != available_sites.end(); site_iter++) {
        //size_t mono_id = site_iter->getMonosaccharideID();
        mod_dist[*site_iter] = avg;
      }
      std::cout << "Highly modified! AVG NUMBER: " << avg << std::endl;
    } else {

      for(ModificationSites::iterator site_iter = available_sites.begin(); site_iter != available_sites.end(); site_iter++) {
        size_t mono_id = site_iter->getMonosaccharideID();
        mod_dist[*site_iter] = avg0 / (double)residue_map[mono_id];
        std::cout << "Lowly modified! AVG NUMBER: " << mod_dist[*site_iter] << std::endl; 
      }

      
    }
    



    for(ModificationSites::iterator site_iter = diff_sites.begin(); site_iter != diff_sites.end(); site_iter++)
    {
      mod_dist[*site_iter] = 1e-6;
    }

    //mod_pattern[mod_symbol] = mod_dist;

    std::cout << mod_symbol << "\t" << total_num << "\t" <<  complete_sites.size() << "\n\n";
  }

  double SequencePrediction::inferFromGraph( CleavageItem& clv, CleavageTree& clv_tree, bool complemenatry /*= true*/ )
  {
    double confidence_constant = param.getParameter<double>("intactness_constant").first;
    int total_mod_num = gs->getModificationConstraint(mod_symbol);

    CleavageNodePtr parent = *(clv_tree.getSupersets(clv).begin());
    CleavageNodePtr child = *(clv_tree.getSubsets(clv).begin());

    // Arbitrarily choose the one with the largest confidence.
    const std::vector<CleavageItem>& child_items = child->getCleavageItems();
      
    double child_conf(0.0);
    int child_index(0);
    double max_child_conf(0.0);
    std::cout << child_items.size() << " Children!\n";

    int total_size(0);
    for(size_t i = 0; i < child_items.size(); i++)
      total_size += child_items[i].getRepeatNumber();

    std::cout << "Total size: " << total_size << "\n";

    for(size_t i = 0; i < child_items.size(); i++)
    {
      // Get the confidence value.
      double conf = child_items[i].getConfidenceValue();
      int mem_size = child_items[i].getRepeatNumber();

      std::cout << "Child " << i << ": " << conf << "\n";
      std::cout << "Member size: " << mem_size << "\n";
      printModificationSites(child_items[i].mod_sites);
      std::cout << mod_symbol << "\t" << child_items[i].mod_num << "\n";

      if(conf > max_child_conf) {
        child_index = (int)i;
        max_child_conf = conf * (mem_size/total_size);
      }      
    }
    
    const CleavageItem& child_item = child_items[child_index];
       
    int up_bound = child_item.mod_num + (int)(clv.mod_sites.size() - child_item.mod_sites.size());

    int low_bound = child_item.mod_num;
              
    int child_dist(0);
    if(clv.mod_num < low_bound || clv.mod_num > up_bound) {
      // Unaccepted situation. Punish the confidence.
      child_dist = -1;
    } else if(up_bound > total_mod_num) {
      child_dist = total_mod_num - clv.mod_num;
    } else {
      child_dist = up_bound - clv.mod_num;
    }

    //double conf = const_iter->getConfidenceValue() / conf_sum;
    //std::cout << "Max child confidence: " << max_child_conf << "\tChild dist: " << child_dist << " ";
      
    if(child_dist < 0) {
      //child_dist = total_mod_num * 10;
      //child_conf += max_child_conf * pow(confidence_constant, child_dist);
      child_conf = -1.0;
    } else if(child_dist == 0) {
      child_conf += max_child_conf * pow(confidence_constant, (double)child_dist + 1e-3);
    } else {
      child_conf += max_child_conf * pow(confidence_constant, child_dist);
    }

    // Reset conf_sum;
    double parent_conf(0.0);
    int parent_index(0);
    //max_conf = 0.0;
    double max_parent_conf(0.0);
    
    // Reset total_size
    total_size = 0;

    const std::vector<CleavageItem>& parent_items = parent->getCleavageItems();
    for(size_t i = 0; i < parent_items.size(); i++)
      total_size += parent_items[i].getRepeatNumber();

    for(size_t i = 0; i < parent_items.size(); i++)
    {
      // Get the confidence value.
      double conf = parent_items[i].getConfidenceValue();
      int mem_size = parent_items[i].getRepeatNumber();

      if(conf > max_parent_conf) {
        parent_index = (int)i;
        max_parent_conf = conf * (mem_size/total_size);
      }      
    }

    up_bound = parent_items[parent_index].mod_num;

    low_bound = parent_items[parent_index].mod_num - (int)(parent_items[parent_index].mod_sites.size() - clv.mod_sites.size());

    int parent_dist(0);
    if(clv.mod_num < low_bound || clv.mod_num > up_bound) {
      // Unaccepted situation. Punish the confidence.
      parent_dist = -1;
    } else {
      parent_dist = up_bound - clv.mod_num; 
    }

    if(parent_dist < 0) {
      //parent_dist = total_mod_num * 10;
      //parent_conf += max_parent_conf * pow(confidence_constant, parent_dist);
      parent_conf = -1.0;
    } else if(parent_dist == 0) {
      parent_conf += max_parent_conf * pow(confidence_constant, (double)parent_dist + 1e-3);
    } else {
      parent_conf += max_parent_conf * pow(confidence_constant, parent_dist);
    }

    //std::cout << "Parent dist: " << parent_dist << " ";

    // Check the golden pair.
    double comp_conf(0.0);
    
    // If complementary is false, the comp node created by clv will be ignored.
    if(!complemenatry) {
      // Do not check the complementary node, simply assign it a high value.
      int comp_dist = 0;
      comp_conf = pow(confidence_constant, (double)comp_dist);
      //std::cout << "Comp dist: " << comp_dist << "\n";  
    
    } else {
      ModificationSites comp_sites = this->getComplementaryModificationSites(clv.mod_sites);

      CleavageNodePtr comp_node = clv_tree.locateModificationSites(comp_sites);


      if(comp_node == nullptr) {
        std::cout << "No complementary node!\n";
        double comp_dist = (double)total_mod_num - (double)clv.mod_num;

        if(comp_dist == 0.0) comp_dist = 1e-3;

        comp_conf = pow(confidence_constant, comp_dist);
        //std::cout << "Comp dist: " << comp_dist << " ";
      } else {
        const std::vector<CleavageItem>& comp_items = comp_node->getCleavageItems();

        double max_comp_conf(0.0);
        int comp_index(0);
        for(size_t i = 0; i < comp_items.size(); i++)
        {
          // Get the confidence value.
          double conf = comp_items[i].getConfidenceValue();
          if(conf > max_comp_conf) {
            comp_index = (int)i;
            max_comp_conf = conf;
          }      
        }

        int sum = comp_items[comp_index].mod_num + clv.mod_num;

        double comp_dist(1e-3);
        if(sum > total_mod_num) {
          // Not golden pair. Not interested.
          std::cout << "Unaccepted complementary node!\n";
          comp_dist = -1.0;
          //std::cout << "Comp dist: " << comp_dist << " ";
          comp_conf = -1.0;
        } else if(sum < total_mod_num) {
          comp_dist = (double) total_mod_num * 10;
          //std::cout << "Comp dist: " << comp_dist << " ";
          comp_conf += max_comp_conf * pow(confidence_constant, comp_dist);
        } else {
          //comp_dist = 1e-3;
          //std::cout << "Comp dist: " << comp_dist << " ";
          comp_conf += max_comp_conf * pow(confidence_constant, comp_dist);
        }

        
      }
      //std::cout << "\n";
    }

    double dist_array[] = {child_conf, parent_conf, comp_conf};
    
    double min_conf = *std::min_element(dist_array, dist_array+3);
    //printModificationSites(clv.mod_sites);
    std::cout << "\nChild confidence: " << child_conf << "\tParent confidence: " << parent_conf << "\tComp confidence: " << comp_conf << "\n" << std::endl;
    std::cout << "Min confidence: " << min_conf << "\n";
    
    // The constant should be able to be controlled by the parameter table.
    return min_conf;
  }

  AssignmentSpace SequencePrediction::groupAssignmentsByFragmentType( const AssignmentSpace& matched_results, bool missing /* = true */ )
  {
    std::cout << "Group assignments by fragment type...\n";
    
    // Notice that B and C ion should not be differentiated.
    typedef std::pair<std::pair<std::string, FragmentPosition>, int> KeyType;
    typedef std::map<KeyType, std::map<MonoPeakPtr, NodeItem>> FragmentCluster;
    FragmentCluster iid_evidence;

    for(AssignmentSpace::const_iterator assign_iter = matched_results.begin(); assign_iter != matched_results.end(); assign_iter++) {
      const FragmentPtr temp_frag = assign_iter->second.getFragment();
      
      // Remove unqualified assignments.
      if(!this->isQualifiedAssignment(*temp_frag)) continue;

      int mod_num = (missing ? -1 : assign_iter->second.getModificationNum(mod_symbol));

      FragmentPosition frag_pos = temp_frag->getCleavageCollection().begin()->second;
      std::string frag_type = temp_frag->getCleavageCollection().begin()->first;
      if(frag_type == "B") {
        frag_type = "C";
      } else if(frag_type == "Z") {
        frag_type = "Y";
      }
      KeyType cluster_key = std::make_pair(std::make_pair(frag_type, frag_pos), mod_num);
      iid_evidence[cluster_key].insert(*assign_iter);
    }

    // For each group, take the one with the smallest value of shift.
    AssignmentSpace frag_center;
    for(FragmentCluster::iterator iter = iid_evidence.begin(); iter!= iid_evidence.end(); iter++) {
      std::map<MonoPeakPtr, NodeItem>& satellite = iter->second;

      std::map<MonoPeakPtr, NodeItem>::iterator best_choice = satellite.begin();
      
      for(std::map<MonoPeakPtr, NodeItem>::iterator sate_iter = satellite.begin(); sate_iter != satellite.end(); sate_iter++) {
        std::cout << "MZ: " << sate_iter->first->mz << "\tCharge: " << sate_iter->first->z << "\tType: " << sate_iter->second.getCleavageType() << "\n";

        if(!missing) {  
          // C type over B type, Y over Z type.
          if((sate_iter->second.getCleavageClass() == "C" && best_choice->second.getCleavageClass() == "B") || (sate_iter->second.getCleavageClass() == "Y" && best_choice->second.getCleavageClass() == "Z")) {
            std::cout << "Type privilege\n";
            best_choice = sate_iter;
          } else if(sate_iter->second.getCleavageClass() == best_choice->second.getCleavageClass()) {
            if(sate_iter->second.getCompositionShiftList().size() < best_choice->second.getCompositionShiftList().size()) {
              std::cout << "Better shift!\n";
              best_choice = sate_iter;
            }
          }
        } else {
          // The one with the highest number of modification.
          if(sate_iter->second.getModificationNum(mod_symbol) > best_choice->second.getModificationNum(mod_symbol)) {
            std::cout << "Higher mod number: " << sate_iter->second.getModificationNum(mod_symbol) << " -> " << best_choice->second.getModificationNum(mod_symbol);
            best_choice = sate_iter;
          } else if(sate_iter->second.getCompositionShift() < best_choice->second.getCompositionShift()) {
            std::cout << "Better shift!\n";
            best_choice = sate_iter;
          }
        }
        
      }

      frag_center.insert(*best_choice);
    }

    return frag_center;

  }

  AssignmentSpace SequencePrediction::groupAssignmentsByModificationSites( AssignmentSpace& evidence, bool missing /*= true*/ )
  {
    std::cout << "Group assignments by modification sites...\n";
    
    // Group assignments by modification sites.  The assignments for grouping should have identical uniqueness value.
    // For Ac, the modification number should be recorded; For SO3, the top one is taken.
    typedef std::map<ModificationSites, std::multimap<MonoPeakPtr, NodeItem>> ModificationCluster;
    ModificationCluster fold_evidence;

    for(AssignmentSpace::iterator clv_iter = evidence.begin(); clv_iter != evidence.end(); clv_iter++) {
      
      ModificationSites mod_sites = clv_iter->second.getFragment()->getModificationSitesBySymbol(mod_symbol,1);
      
      fold_evidence[mod_sites].insert(*clv_iter);
    }

    AssignmentSpace site_center;
    for(ModificationCluster::iterator mod_iter = fold_evidence.begin(); mod_iter != fold_evidence.end(); mod_iter++) {
      AssignmentSpace& mod_cluster = mod_iter->second;
      if(!missing) {
        // Simply append the information to site_center.
        for(AssignmentSpace::iterator assign_iter = mod_cluster.begin(); assign_iter!= mod_cluster.end();assign_iter++)
          site_center.insert(*assign_iter);
      } else {
        AssignmentSpace::iterator best_choice = mod_cluster.begin();
        for(AssignmentSpace::iterator assign_iter = mod_cluster.begin(); assign_iter != mod_cluster.end(); assign_iter++) {
          if(assign_iter->second.getModificationNum(mod_symbol) > best_choice->second.getModificationNum(mod_symbol)) {
            best_choice = assign_iter;
          }
        }
        site_center.insert(*best_choice);
      }
    }
    //if(missing) {
    //  // Simply get the one with the highest number of SO3.
    //  for(ModificationCluster::iterator mod_iter = fold_evidence.begin(); mod_iter != fold_evidence.end(); mod_iter++) {
    //    AssignmentSpace& mod_cluster = mod_iter->second;
    //    AssignmentSpace::iterator best_choice = mod_cluster.begin();
    //    for(AssignmentSpace::iterator assign_iter = mod_cluster.begin(); assign_iter != mod_cluster.end(); assign_iter++) {
    //      if(assign_iter->second.getModificationNum(mod_symbol) > best_choice->second.getModificationNum(mod_symbol)) {
    //        best_choice = assign_iter;
    //      }
    //    }
    //    site_center.insert(*best_choice);
    //  }
    //  
    //} else {
    //  // Separating assignments based on the number of mod.
    //  for(ModificationCluster::iterator mod_iter = fold_evidence.begin(); mod_iter != fold_evidence.end(); mod_iter++) {
    //    //std::map<int, AssignmentSpace> multi_evi;
    //    std::map<int, std::set<FragmentPosition>> score_map;
    //    AssignmentSpace& mod_cluster = mod_iter->second;
    //    for(AssignmentSpace::iterator assign_iter = mod_cluster.begin(); assign_iter != mod_cluster.end(); assign_iter++) {
    //      int mod_num = assign_iter->second.getModificationNum(mod_symbol);

    //      FragmentPosition& frag_pos = assign_iter->second.getFragmentPosition();

    //      // If the fragment position has been recorded, skip the recording.
    //      std::set<FragmentPosition>& frag_record = score_map[mod_num];
    //      if(frag_record.find(frag_pos) == frag_record.end()) {
    //        score_map[mod_num].insert(frag_pos);
    //        site_center.insert(*assign_iter);
    //      } 
    //    }
    //  }
    //}
     
    return site_center;
  }

  double SequencePrediction::estimateUniquenessValue(AssignmentSpace::iterator ele, const AssignmentSpace& pool) 
  {
    std::cout << "\nEstimating uniqueness value:\n\n"; 

    std::cout << "MZ: " << ele->first->mz << "\tCharge: " << ele->first->z << "\n";

    //std::set<std::pair<ModificationSites, int>> gb_site_vec;
    //std::set<std::pair<ModificationSites, int>> cr_site_vec;
    //// Currently the GG, GC, CG are treated as the same.
    //std::set<std::pair<ModificationSites, int>> internal_site_vec;
    //// CC are treated separately.
    //std::set<std::pair<ModificationSites, int>> cc_site_vec;

    AssignmentClassification site_map;

    std::pair<AssignmentSpace::const_iterator, AssignmentSpace::const_iterator> p = pool.equal_range(ele->first);

    while(p.first != p.second) {
      const NodeItem& node = p.first->second;
      const ModificationSites& mod_sites = node.getFragment()->getModificationSitesBySymbol(mod_symbol,1);
      int mod_num = node.getModificationNum(mod_symbol);
      std::cout << "General type: " << node.getGeneralType() << "\n";
      /*if(node.getGeneralType() == "G") {
      gb_site_vec.insert(std::make_pair(mod_sites, mod_num));
      } else if(node.getGeneralType() == "C") { 
      cr_site_vec.insert(std::make_pair(mod_sites, mod_num));
      } else if(node.getGeneralType() == "CC") {
      cc_site_vec.insert(std::make_pair(mod_sites, mod_num));
      } else if(node.getGeneralType() != "SEQ") {
      internal_site_vec.insert(std::make_pair(mod_sites, mod_num));
      }*/

      if(node.getGeneralType() != "SEQ") {
        site_map[node.getGeneralType()].insert(std::make_pair(mod_sites, mod_num));
      }
      p.first++;
    }

    if(site_map["G"].size() == 0 && site_map["C"].size() == 0) return 0;
    
    //if(gb_site_vec.size() == 0 && cr_site_vec.size() == 0) return 0.0;

    std::string type = ele->second.getGeneralType();
    double uni_score = calculateSiteUniqueness(type, site_map);

    return uni_score;

  }

  double SequencePrediction::calculateSiteUniqueness(std::string type,  const AssignmentClassification& assignments )
  {
    // Check glycosidic-bond cleavage.
    double total_score = 0.0;
    //std::set<std::pair<ModificationSites, int>> total_mod_set;
    for(std::vector<std::string>::iterator iter = type_vec.begin(); iter != type_vec.end(); iter++)
    {
      std::string temp_type = *iter;
      AssignmentClassification::const_iterator class_iter = assignments.find(temp_type);
      if(class_iter == assignments.end()) continue;

      for(std::set<std::pair<ModificationSites, int>>::const_iterator assign_iter = class_iter->second.begin(); assign_iter != class_iter->second.end(); assign_iter++) {
        //const std::pair<ModificationSites,int>& key = *assign_iter;
        //std::set<std::pair<ModificationSites, int>>::iterator mod_iter = total_mod_set.find(key);
        
        total_score += type_score[temp_type];

        /*if(mod_iter == total_mod_set.end()) {
          total_mod_set.insert(key);
          total_score += type_score[temp_type];
        } else {
          std::cout << "Redundant and reduced!\n";
        }*/
      }

    }

    return type_score[type] / (1e-3 + total_score);
    //double weight0 = param.getParameter<double>("gb_cleavage").first;
    //for(std::set<std::pair<ModificationSites, int>>::const_iterator iter = assignments['G'].begin(); iter != assignments['G'].end();  iter++) {
    //  total_mod_set.insert(*iter);
    //  total_score += weight0;
    //}

    //// For the rest of the assignments, check and see if they can be reduced to the total_mod_set.
    //double weight1 = param.getParameter<double>("cr_cleavage").first; // Weight for cross-ring cleavage.
    //for(std::set<std::pair<ModificationSites, int>>::const_iterator iter = assignments['C'].begin(); iter != assignments['C'].end();  iter++) {
    //  std::set<std::pair<ModificationSites, int>>::iterator mod_iter = total_mod_set.find(*iter);
    //  // If found, do nothing. Otherwise, add to the total score.
    //  if(mod_iter == total_mod_set.end()) {
    //    total_mod_set.insert(*iter);
    //    total_score += weight1;
    //  } else {
    //    std::cout << "Redundant and reduced!\n";
    //  }
    //}

    //double weight2 = param.getParameter<double>("int_cleavage").first; // Weight for internal cleavage.
    //for(std::set<std::pair<ModificationSites, int>>::iterator iter = assignments[].begin(); iter != internal_site.end();  iter++) {
    //  std::set<std::pair<ModificationSites, int>>::iterator mod_iter = total_mod_set.find(*iter);
    //  if(mod_iter == total_mod_set.end()) {
    //    total_mod_set.insert(*iter);
    //    total_score += weight2;
    //  } else {
    //    std::cout << "Redundant and reduced!\n";
    //  }
    //}


    //if(type == "G")
    //  return weight0/total_score;
    //else if(type == "C")
    //  return weight1/total_score;
    //else
    //  throw std::runtime_error("Strange cleavage type.");


  }

  bool SequencePrediction::generateModificationDistributionV2( std::string symbol, const AssignmentSpace& pool, bool missing /*= false*/ )
  {
    // Update the internal mod status.
    mod_symbol = symbol;

    std::cout << "Generate modification distribution V2 -- " << mod_symbol << "\n";

    int total_mod_num = gs->getModificationConstraint(mod_symbol);
    if(total_mod_num == 0) {
      std::cout << "No " << mod_symbol << " groups\n";
      return false;
    }

    ModificationSites complete_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);
    //ModificationSites empty_sites;

    // 1. Group assignments by fragment types. For Ac, assignments with the same site but different mod number are grouped separately. B and C, Z and Y are treated as the same. The multiple charge problem has also been processed at this point.
    AssignmentSpace assignment_rep = this->groupAssignmentsByFragmentType(pool, missing);

    std::set<MonoPeakPtr> mono_pk_set;
    for(AssignmentSpace::iterator iter = assignment_rep.begin(); iter != assignment_rep.end(); iter++) {
      mono_pk_set.insert(iter->first);
    }

    // 2. Grouping assignments by the mono pk, and Evaluate the uniqueness of assignments
    // The first key is the uniqueness value;
    // The second is the mass value.
    //typedef std::map<double, std::map<MonoPeakPtr, std::vector<NodeItem>>> AssignmentCluster;

    //AssignmentCluster uni_assignments;

    // Uniqueness and 
    std::map<double, std::set<MonoPeakPtr>> uni_map;
    std::map<MonoPeakPtr, AssignmentClassification> uni_abstract_map;

    // For each peak, get the node items from pool.
    for(std::set<MonoPeakPtr>::iterator mono_iter = mono_pk_set.begin(); mono_iter != mono_pk_set.end(); mono_iter++) {

      std::cout << "MZ: " << (*mono_iter)->mz << "\tCharge: " << (*mono_iter)->z << "\n";
      std::pair<AssignmentSpace::const_iterator, AssignmentSpace::const_iterator> p = pool.equal_range(*mono_iter);

      AssignmentClassification abstract_assigns = this->generateUniqueAssignments(p);

      uni_abstract_map.insert(std::make_pair(*mono_iter, abstract_assigns));

      for(AssignmentSpace::const_iterator assign_iter = p.first; assign_iter != p.second; assign_iter++)
      {
        std::string type = assign_iter->second.getGeneralType();

        if(type == "CC") continue;

        // Calculate the uniqueness value.
        double uniqueness_value = this->calculateSiteUniqueness(type, abstract_assigns);

        uni_map[uniqueness_value].insert(*mono_iter);

        //p.first++;
      }
      

    }

    // Initialize modification distribution.
    this->initialize();

    // Convert the distribution into log version for multiplication.
    ModificationDistribution& prior = mod_pattern[mod_symbol];
    //ModificationDistribution log_prior = mod_pattern[mod_symbol];
    
    
    /*for(ModificationDistribution::iterator dist_iter = log_prior.begin(); dist_iter != log_prior.end(); dist_iter++) {
      dist_iter->second = log(dist_iter->second);
    }*/
    int current_size(0);
    // 3. Sort mass values by their uniqueness values.
    for(std::map<double, std::set<MonoPeakPtr>>::reverse_iterator rev_iter = uni_map.rbegin(); rev_iter != uni_map.rend(); rev_iter++) {
      std::cout << "Uniqueness value: " << rev_iter->first << "\n";
      std::set<MonoPeakPtr>& mono_pk_list = rev_iter->second;
      
      current_size += (int) mono_pk_list.size();
      ModificationDistribution temp_dist;
      
      for(std::set<MonoPeakPtr>::iterator mass_iter = mono_pk_list.begin(); mass_iter != mono_pk_list.end(); mass_iter++) 
      {
        std::cout << "MZ: " << (*mass_iter)->mz << "\tCharge: " << (*mass_iter)->z << "\n";
        AssignmentClassification& abstract_assigns = uni_abstract_map[*mass_iter];
        //std::cout << "The number of abstract assignments: " << abstract_assigns.size() << "\n";

        if(abstract_assigns.size() == 0) continue;
        // 
        if(abstract_assigns.size() == 1 && abstract_assigns.begin()->second.size() == 1){
          if(*(abstract_assigns.begin()->second.begin()) == std::make_pair(complete_sites, total_mod_num) || abstract_assigns.begin()->second.begin()->first.size()==0) {
            std::cout << "No contribution to the understanding!\n";
            continue;
          }
        }

        // 3.1 Convert the mass into a distribution.
        ModificationDistribution likelihood = this->updateDistributionV3(abstract_assigns, prior);

        printModificationDistribution(likelihood);

        // 3.2 Integrate the distribution into the background distribution.
        if(temp_dist.empty() == true) {
          temp_dist = likelihood;
        } else {
          for(ModificationDistribution::iterator dist_iter = temp_dist.begin(); dist_iter != temp_dist.end(); dist_iter++)
          {
            dist_iter->second += likelihood[dist_iter->first];
          }
        }
        //for(ModificationDistribution::iterator dist_iter = log_prior.begin(); dist_iter != log_prior.end(); dist_iter++) 
        //{
        //  dist_iter->second = likelihood[dist_iter->first] + dist_iter->second;
        //}
      
      }

      // 4. Normalization for next circle of updating.
      double sum(0.0);
      for(ModificationDistribution::iterator dist_iter = temp_dist.begin(); dist_iter != temp_dist.end(); dist_iter++) {
        //dist_iter->second = exp(dist_iter->second);
        sum += dist_iter->second;
      }

      double origial_sum(0.0);
      for(ModificationDistribution::iterator dist_iter = temp_dist.begin(); dist_iter != temp_dist.end(); dist_iter++) {
        //ModificationPosition pos = dist_iter->first;
        prior[dist_iter->first] += dist_iter->second/(sum * (double)current_size);
        origial_sum += prior[dist_iter->first];
      }

      for(ModificationDistribution::iterator dist_iter = temp_dist.begin(); dist_iter != temp_dist.end(); dist_iter++) {
        prior[dist_iter->first] /= origial_sum;
      }


      std::cout << "\nUpdated modification distribution!\n";
      // Display the modification distribution.
      printModificationDistribution(mod_symbol);
    }
    return true;
  }

  ModificationDistribution SequencePrediction::updateDistributionV2( const AssignmentClassification& uni_set, const ModificationDistribution& dist )
  {
    int total_num = gs->getModificationConstraint(mod_symbol);
    ModificationSites& complete_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);

    ModificationDistribution likelihood_distribution;

    // Collect the covered sites by the abstract assignments.
    ModificationSites covered_sites;
    
    double total_size(0.0);
    for(AssignmentClassification::const_iterator const_iter = uni_set.begin(); const_iter != uni_set.end(); const_iter++) {
      total_size += (double)const_iter->second.size() * type_score[const_iter->first];
    }
    std::cout << "Total score: " << total_size << "\n";

    // 1. Evaluate the relative importance of elements.
    for(std::vector<std::string>::const_iterator const_iter = type_vec.begin(); const_iter != type_vec.end(); const_iter++) 
    {
      // Consider only G and C type of cleavage.
      if(*const_iter == "CC") continue;
           
      // Get the assignments corresponding to specific type.
      
      AssignmentClassification::const_iterator uni_iter = uni_set.find(*const_iter);
      if(uni_iter == uni_set.end()) continue;

      std::cout << "Type: " << uni_iter->first << "\n";
      
      BOOST_FOREACH(auto p, uni_iter->second) {

        ModificationSites assign_sites = p.first;
        int current_num = p.second;
        std::cout << "Current sites:\n";
        printModificationSites(assign_sites);
        std::cout << "Current num: " << current_num << "\n";

        // Conversion to the complementary part.
        if(p.second == 0) {
          assign_sites = getDiffSites(assign_sites, complete_sites);
          current_num = total_num - current_num;
          std::cout << "Converted!\n";
        }
        
        if(*const_iter != "G" && *const_iter != "C") {
          std::cout << "Pass internal cleavages!\n";
          continue;
        }

        double avg = (double) (current_num) * type_score[uni_iter->first]/ (total_size * (double)(assign_sites.size()));
        
        // Append assign_sites to covered_sites.
        covered_sites.insert(assign_sites.begin(), assign_sites.end());
        // For each position covered by the assignments, sum the avg number.
        BOOST_FOREACH(auto s, assign_sites)
        {
          if(likelihood_distribution.find(s) == likelihood_distribution.end()) {
            likelihood_distribution.insert(std::make_pair(s, avg));
          } else {
            likelihood_distribution[s] += avg;
          } 
        }
      }  
      
    }

    // For the rest of the assignments, keep the local distribution from the background.
    ModificationSites rest_sites = this->getDiffSites(covered_sites, complete_sites);
    double covered_sum(0.0);
    BOOST_FOREACH(auto s, covered_sites)
    {
      covered_sum += likelihood_distribution[s];
    }

    double rest_sum = (double)total_num - covered_sum;

    double sum(0.0);
    BOOST_FOREACH(auto s, rest_sites)
    {
      ModificationDistribution::const_iterator const_iter = dist.find(s);

      sum += const_iter->second + 1e-3;
    }

    BOOST_FOREACH(auto s, rest_sites)
    {
      ModificationDistribution::const_iterator const_iter = dist.find(s);

      double relative_dist = (const_iter->second + 1e-3)/sum;
      likelihood_distribution[s] = rest_sum * relative_dist;
    }

    return likelihood_distribution;
  }

  AssignmentClassification SequencePrediction::generateUniqueAssignments( const std::pair<AssignmentSpace::const_iterator, AssignmentSpace::const_iterator>& assign_range )
  {
    
    std::map<std::string, AssignmentSpace> type_cluster;
    
    // Organize assignments by cleavage type.
    for(AssignmentSpace::const_iterator assign_iter = assign_range.first; assign_iter != assign_range.second; assign_iter++) {
      const NodeItem& node = assign_iter->second;
      //const ModificationSites& mod_sites = node.getFragment()->getModificationSitesBySymbol(mod_symbol,1);
      //int mod_num = node.getModificationNum(mod_symbol);
      std::string current_type = node.getGeneralType();
      std::cout << "General type: " << current_type << "\n";

      if(current_type == "SEQ") continue;

      type_cluster[current_type].insert(*assign_iter);

      /*std::pair<ModificationSites, int> key = std::make_pair(mod_sites, mod_num);
      if(total_mod_set.find(key) == total_mod_set.end())
        site_map[node.getGeneralType()].insert(std::make_pair(mod_sites, mod_num));*/
      
      //assign_range.first++;
    }

    ModificationSites& complete_sites = gs->getModificationSitesBySymbol
(mod_symbol, 1);
    int total_num = gs->getModificationConstraint(mod_symbol);

    AssignmentClassification site_map;
    std::set<std::pair<ModificationSites, int>> total_mod_set;
    BOOST_FOREACH(auto t, type_vec)
    {
      auto iter = type_cluster.find(t);
      if(iter == type_cluster.end()) continue;

      AssignmentSpace& storage = iter->second;
      BOOST_FOREACH(auto i, storage)
      {
        const NodeItem& node = i.second;
        const ModificationSites& mod_sites = node.getFragment()->getModificationSitesBySymbol(mod_symbol,1);
        int mod_num = node.getModificationNum(mod_symbol);

        // Unlikely explanation.
        if(mod_sites == complete_sites && mod_num != total_num)
          continue;

        std::pair<ModificationSites, int> key = std::make_pair(mod_sites, mod_num);
        if(total_mod_set.find(key) == total_mod_set.end()) {
          total_mod_set.insert(key);
          site_map[node.getGeneralType()].insert(key);
        } else {
          std::cout << "Redundant type information!\n";
        }
      }

    }

    return site_map;
  }

  ModificationDistribution SequencePrediction::updateDistributionV3( const AssignmentClassification& uni_set, const ModificationDistribution& dist )
  {
    int total_num = gs->getModificationConstraint(mod_symbol);
    ModificationSites& complete_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);

    ModificationDistribution likelihood_distribution;

    typedef std::map<std::pair<ModificationSites, int>, double> DistanceMap;
    DistanceMap distance_map;

    double total_distance(0.0);
    // Iterate over all cleavage types.
    for(std::vector<std::string>::iterator type_iter = type_vec.begin(); type_iter != type_vec.end(); type_iter++)
    {
      std::cout << "Type: " << *type_iter << "\n";
      if(*type_iter == "CC") continue;

      AssignmentClassification::const_iterator const_iter = uni_set.find(*type_iter);
      if(const_iter == uni_set.end()) continue;

      std::set<std::pair<ModificationSites, int>> clv_set = const_iter->second;

      // Calculate the distance to the prior.
      BOOST_FOREACH(auto clv, clv_set)
      {
        printModificationSites(clv.first);
        std::cout << "Mod num: " << clv.second << "\n";
        // Iterate over all position.
        double bg_sum(0.0);
        BOOST_FOREACH(auto s, clv.first)
        {
          ModificationDistribution::const_iterator const_iter = dist.find(s);
          bg_sum += const_iter->second;
        }
        double bg_distance = type_score[*type_iter]/(abs(clv.second - bg_sum) + 0.01);
        std::cout << "Converted distance: " << bg_distance << "\n";
        distance_map.insert(std::make_pair(clv, bg_distance));
        total_distance += bg_distance;
      }
    }

    std::cout << "Total distance sum: " << total_distance << "\n";

    //ModificationDistribution covered_map;
    ModificationSites covered_sites;
    for(DistanceMap::iterator dist_iter = distance_map.begin(); dist_iter != distance_map.end(); dist_iter++)
    {
      std::pair<ModificationSites, int> assign = dist_iter->first;
      double weight = (dist_iter->second + 1e-3)/total_distance;
      ModificationSites assign_sites = assign.first;
      int current_num = assign.second;
      
      // The avg information on each modification position.
      if(current_num == 0) {
        assign_sites = getDiffSites(assign_sites, complete_sites);
        current_num = total_num - current_num;
        std::cout << "Converted!\n";
      }
      
      double avg = (double)current_num * weight /(double) assign_sites.size();

      printModificationSites(assign_sites);
      std::cout << "Num: " << current_num << "\tWeight: " << weight << "\tAvg: " << avg << "\n";

      // Iterate over all position.
      BOOST_FOREACH(auto pos, assign_sites) {
        covered_sites.insert(pos);
        if(likelihood_distribution.find(pos) == likelihood_distribution.end()) {
          likelihood_distribution.insert(std::make_pair(pos, avg));
        } else {
          likelihood_distribution[pos] += avg;
        }
      }
    }

    ModificationSites rest_sites = this->getDiffSites(covered_sites, complete_sites);
    double covered_sum(0.0);
    BOOST_FOREACH(auto s, covered_sites)
    {
      covered_sum += likelihood_distribution[s];
    }

    double rest_sum = (double)total_num - covered_sum;

    double sum(0.0);
    BOOST_FOREACH(auto s, rest_sites)
    {
      ModificationDistribution::const_iterator const_iter = dist.find(s);

      sum += const_iter->second + 1e-3;
    }

    BOOST_FOREACH(auto s, rest_sites)
    {
      ModificationDistribution::const_iterator const_iter = dist.find(s);

      double relative_dist = (const_iter->second + 1e-3)/sum;
      likelihood_distribution[s] = rest_sum * relative_dist;
    }

    return likelihood_distribution;


    //return likelihood_distribution;

    
  }

  bool SequencePrediction::generateModificationDistributionV3( std::string symbol, const AssignmentSpace& pool, bool missing /*= false*/ )
  {
    // Update the internal mod status.
    mod_symbol = symbol;

    std::cout << "Generate modification distribution V3 -- " << mod_symbol << "\n";

    int total_mod_num = gs->getModificationConstraint(mod_symbol);
    if(total_mod_num == 0) {
      std::cout << "No " << mod_symbol << " groups! Skip the step\n";
      return false;
    }

    ModificationSites complete_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);
    //ModificationSites empty_sites;

    // 1. Group assignments by fragment types. For Ac, assignments with the same site but different mod number are grouped separately. B and C, Z and Y are treated as the same. The multiple charge problem has also been processed at this point.
    AssignmentSpace assignment_rep = this->groupAssignmentsByFragmentType(pool, missing);

    std::set<MonoPeakPtr> mono_pk_set;
    for(AssignmentSpace::iterator iter = assignment_rep.begin(); iter != assignment_rep.end(); iter++) {
      mono_pk_set.insert(iter->first);
    }

    // 2. Grouping assignments by the mono pk, and Evaluate the uniqueness of assignments
    // The first key is the uniqueness value;
    // The second is the mass value.
    //typedef std::map<double, std::map<MonoPeakPtr, std::vector<NodeItem>>> AssignmentCluster;

    // Uniqueness value and corresponding peaks.
    std::map<double, std::set<MonoPeakPtr>> uni_map;
    std::map<MonoPeakPtr, AssignmentClassification> uni_abstract_map;

    // For each peak, get the node items from pool.
    for(std::set<MonoPeakPtr>::iterator mono_iter = mono_pk_set.begin(); mono_iter != mono_pk_set.end(); mono_iter++) {

      std::cout << "MZ: " << (*mono_iter)->mz << "\tCharge: " << (*mono_iter)->z << "\n";
      std::pair<AssignmentSpace::const_iterator, AssignmentSpace::const_iterator> p = pool.equal_range(*mono_iter);

      AssignmentClassification abstract_assigns = this->generateUniqueAssignments(p);

      uni_abstract_map.insert(std::make_pair(*mono_iter, abstract_assigns));

      for(AssignmentSpace::const_iterator assign_iter = p.first; assign_iter != p.second; assign_iter++)
      {
        std::string type = assign_iter->second.getGeneralType();

        if(type == "CC") continue;

        // Calculate the uniqueness value.
        double uniqueness_value = this->calculateSiteUniqueness(type, abstract_assigns);

        uni_map[uniqueness_value].insert(*mono_iter);
      }
      

    }

    // Initialize modification distribution.
    this->initialize();

    // Convert the distribution into log version for multiplication.
    ModificationDistribution& prior = mod_pattern[mod_symbol];
    
    /*for(ModificationDistribution::iterator dist_iter = log_prior.begin(); dist_iter != log_prior.end(); dist_iter++) {
      dist_iter->second = log(dist_iter->second);
    }*/
    int current_size(0);
    // 3. Sort mass values by their uniqueness values.
    for(std::map<double, std::set<MonoPeakPtr>>::reverse_iterator rev_iter = uni_map.rbegin(); rev_iter != uni_map.rend(); rev_iter++) {
      std::cout << "Uniqueness value: " << rev_iter->first << "\n";
      std::set<MonoPeakPtr>& mono_pk_list = rev_iter->second;
      
      current_size += (int) mono_pk_list.size();
      ModificationDistribution temp_dist;
      
      for(std::set<MonoPeakPtr>::iterator mass_iter = mono_pk_list.begin(); mass_iter != mono_pk_list.end(); mass_iter++) 
      {
        std::cout << "MZ: " << (*mass_iter)->mz << "\tCharge: " << (*mass_iter)->z << "\n";
        AssignmentClassification& abstract_assigns = uni_abstract_map[*mass_iter];

        if(abstract_assigns.size() == 0) continue;
        // 
        if(abstract_assigns.size() == 1 && abstract_assigns.begin()->second.size() == 1){
          if(*(abstract_assigns.begin()->second.begin()) == std::make_pair(complete_sites, total_mod_num) || abstract_assigns.begin()->second.begin()->first.size()==0) {
            std::cout << "No contribution to the understanding!\n";
            continue;
          }
        }

        // 3.1 Convert the mass into a distribution.
        ModificationDistribution likelihood = this->updateDistributionV3(abstract_assigns, prior);

        printModificationDistribution(likelihood);

        // 3.2 Integrate the distribution into the background distribution.
        if(temp_dist.empty() == true) {
          temp_dist = likelihood;
        } else {
          for(ModificationDistribution::iterator dist_iter = temp_dist.begin(); dist_iter != temp_dist.end(); dist_iter++)
          {
            dist_iter->second += likelihood[dist_iter->first];
          }
        }
       
      }

      // 4. Normalization for next circle of updating.
      double sum(0.0);
      for(ModificationDistribution::iterator dist_iter = temp_dist.begin(); dist_iter != temp_dist.end(); dist_iter++) {
        //dist_iter->second = exp(dist_iter->second);
        sum += dist_iter->second;
      }

      double origial_sum(0.0);
      for(ModificationDistribution::iterator dist_iter = temp_dist.begin(); dist_iter != temp_dist.end(); dist_iter++) {
        //ModificationPosition pos = dist_iter->first;
        prior[dist_iter->first] += dist_iter->second/(sum * (double)current_size);
        origial_sum += prior[dist_iter->first];
      }

      for(ModificationDistribution::iterator dist_iter = temp_dist.begin(); dist_iter != temp_dist.end(); dist_iter++) {
        prior[dist_iter->first] /= origial_sum;
      }


      std::cout << "\nUpdated modification distribution!\n";
      // Display the modification distribution.
      printModificationDistribution(mod_symbol);
    }
    return true;
  }

  void SequencePrediction::equalInitialization(ModificationDistribution& mod_dist)
  {
    ModificationSites complete_sites = gs->getModificationSitesBySymbol(mod_symbol);
    ModificationSites available_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);
    ModificationSites diff_sites = this->getDiffSites(available_sites, complete_sites);

    std::cout << "Complete: " << complete_sites.size() << "\tAvailable: " << available_sites.size() << "\tDifference: " << diff_sites.size() << "\n";

    int total_num = gs->getModificationConstraint(mod_symbol);

    for(ModificationSites::iterator site_iter = diff_sites.begin(); site_iter != diff_sites.end(); site_iter++)
    {
      mod_dist[*site_iter] = 1e-6;
    }

    double avg = (double)total_num / (double)available_sites.size();
    for(ModificationSites::iterator site_iter = available_sites.begin(); site_iter != available_sites.end(); site_iter++) {
      mod_dist[*site_iter] = avg;
    }

    std::cout << mod_symbol << "\t" << total_num << "\t" <<  complete_sites.size() << "\n\n";
  }

  void SequencePrediction::averageBySite( ModificationDistribution& mod_dist, const ModificationSites& mod_sites, double total_num )
  {
    double avg = total_num / (double)mod_sites.size();
    for(ModificationSites::const_iterator iter = mod_sites.begin(); iter != mod_sites.end(); iter++) {
     mod_dist[*iter] = avg;
    }
  }

  void SequencePrediction::averageByResidue( ModificationDistribution& mod_dist, const ModificationSites& mod_sites, double total_num )
  {
        // Count the number of sites on each residue.
    //std::map<size_t, int> residue_map;

    //for(ModificationSites::iterator site_iter = bg_sites.begin(); site_iter != bg_sites.end(); site_iter++) {
    //  size_t mono_id = site_iter->getMonosaccharideID();
    //  if(residue_map.find(mono_id)== residue_map.end()) {
    //    residue_map[mono_id] = 1;
    //  } else {
    //    residue_map[mono_id] += 1;
    //  }
    //}

    //double mono_avg1 = updated_mod_num / (double)residue_map.size();
    //
    //// Reallocate the mod number.
    //if(mod_status == "G" && mono_avg1 <= 1.0) { // The averaging is on the residue level.
    //  /*std::cout << "Scale the bg number!\n";
    //  double scale_factor1 = updated_mod_num / bg_num;
    //  for(ModificationSites::iterator mod_iter = bg_sites.begin(); mod_iter != bg_sites.end(); mod_iter++)
    //    mod_dist[*mod_iter] *= scale_factor1;

    //  ModificationSites comp_sites = this->getDiffSites(current_sites, parent_sites);
    //  double scale_factor2 = (total_num - updated_mod_num) / (total_num - bg_num);

    //  for(ModificationSites::iterator mod_iter = comp_sites.begin(); mod_iter != comp_sites.end(); mod_iter++)
    //    mod_dist[*mod_iter] *= scale_factor2;*/

    //  std::cout << "Averaging over the monosaccharide residue." << std::endl;
    //  
    //  for(ModificationSites::iterator site_iter = bg_sites.begin(); site_iter != bg_sites.end(); site_iter++) {
    //    size_t mono_id = site_iter->getMonosaccharideID();
    //    mod_dist[*site_iter] = mono_avg1 / (double)residue_map[mono_id];

    //    std::cout << "ID: " << mono_id << "\t" << mod_dist[*site_iter] << std::endl;
    //  }

    //  std::map<size_t, int> residue_map2;

    //  for(ModificationSites::iterator site_iter = comp_sites.begin(); site_iter != comp_sites.end(); site_iter++) {
    //    size_t mono_id = site_iter->getMonosaccharideID();
    //    if(residue_map2.find(mono_id)== residue_map2.end()) {
    //      residue_map2[mono_id] = 1;
    //    } else {
    //      residue_map2[mono_id] += 1;
    //    }
    //  }

    //  double mono_avg2 = (total_num - updated_mod_num) / (double) residue_map2.size();
    //  for(ModificationSites::iterator site_iter = comp_sites.begin(); site_iter != comp_sites.end(); site_iter++) {
    //    size_t mono_id = site_iter->getMonosaccharideID();
    //    mod_dist[*site_iter] = mono_avg2 / (double)residue_map2[mono_id];

    //    std::cout << "ID: " << mono_id << "\t" << mod_dist[*site_iter] << std::endl;
    //  }

    //} else  { // The averaging is on the site level.
    //  std::cout << "Average mod number over the sites!\n";

    //  for(ModificationSites::iterator mod_iter = bg_sites.begin(); mod_iter != bg_sites.end(); mod_iter++) {
    //    double avg1 = updated_mod_num / (double)bg_sites.size();
    //    mod_dist[*mod_iter] = avg1;
    //    std::cout << "Average 1: " << avg1 << "\n";
    //  }

    //  for(ModificationSites::iterator mod_iter = comp_sites.begin(); mod_iter != comp_sites.end(); mod_iter++) {
    //    double avg2 = (total_num - updated_mod_num) / (double)comp_sites.size();
    //    mod_dist[*mod_iter] = avg2;
    //    std::cout << "Average 2: " << avg2 << "\n";
    //  }

    //} 
  }



 

  //std::multimap<MonoPeakPtr, NodeItem> SequencePrediction::filterData(std::multimap<MonoPeakPtr, NodeItem>& matched_results, const ModificationSequence& seq, bool status )
  //{
  //  std::multimap<MonoPeakPtr, NodeItem>::iterator iter = matched_results.begin();
  //  std::multimap<MonoPeakPtr, NodeItem> new_data;

  //  // Get the modification types to check.
  //  std::set<std::string> mod_set = gs->getModificationTypes();
  //  //for(ModificationSequence::const_iterator const_iter = seq.begin(); const_iter != seq.end(); const_iter++)
  //    //mod_set.insert(const_iter->first);
  //  
  //  for(; iter != matched_results.end(); iter++)
  //  {
  //    // Do not count internal cleavages.
  //    if(iter->second.getCleavageNum() == 2) continue;

  //    // Decide if the number of modification meets with the sequence.
  //    bool pass = true;
  //   
  //    // Iterate over all modification types.
  //    for(std::set<std::string>::iterator mod_iter = mod_set.begin(); mod_iter!=mod_set.end(); mod_iter++) {
  //      // The modification number contained by the fragment.
  //      int mod_num = iter->second.getModificationNum(*mod_iter);

  //      // Check intersection of modification sites from seq and frag. 
  //      ModificationSequence::const_iterator seq_iter = seq.find(*mod_iter);
  //      
  //      // Available sites of the candidate sequence.
  //      ModificationSites mod_sites = seq_iter->second;

  //      // Available sites at this region.
  //      //ModificationSites region_sites = iter->second.getFragment()->getModificationSitesBySymbol(*mod_iter, 1);
		//    ModificationSites region_sites = iter->second.getFragment()->getModificationSitesBySymbol(*mod_iter, 1);

		//    ModificationSites region_mod_sites = getSiteIntersection(mod_sites, region_sites);

  //      int max_mod_num = (int)region_mod_sites.size();

		//    // If status == true, need to specify a range of the modification number.
  //      if(*mod_iter == "Ac") {// No missing 
  //        if(mod_num != max_mod_num) {
  //          pass = false;
  //          break;
  //        }
  //      } else if(*mod_iter == "SO3") { // Missing.
  //        
		//	    if(status) { // Calculate based on residue level.
  //          // Adjust the upper_limit of sulfate number.
  //          std::string clv_type = iter->second.getFragment()->getFragmentType();
  //          if(clv_type == "A" || clv_type == "X") {
  //            ModificationSites expanded_sites = iter->second.getFragment()->getExpandedModificationSites(*mod_iter, 1);
  //            ModificationSites expanded_mod_sites = getSiteIntersection(mod_sites, expanded_sites);
  //            max_mod_num = (int)expanded_mod_sites.size();
  //          }
		//	    } 
  //        if(mod_num > max_mod_num) {
  //          pass = false;
  //          break;
  //        }
  //      }

  //    }

  //    if(pass) {
  //      new_data.insert(*iter);	
  //    }
  //  }
  //  return new_data;
  //}

  double SequencePrediction::calculateCost( const ModificationSequence& seq )
  {
    std::set<std::string> mod_set = gs->getModificationTypes();
    double total_cost = 0.0;
    
    for(std::set<std::string>::iterator iter = mod_set.begin(); iter != mod_set.end(); iter++)
    {
      ModificationSites theo_sites = gs->getModificationSitesBySymbol(*iter,1);
      ModificationSequence::const_iterator seq_iter = seq.find(*iter);
      ModificationPattern::iterator pat_iter = mod_pattern.find(*iter);

      if(seq_iter != seq.end() && pat_iter != mod_pattern.end()) { // Found!
        ModificationSites mod_sites = seq_iter->second;
        // Calculate the distance between the specified modification sites and the distribution.

        for(ModificationSites::iterator site_iter = theo_sites.begin(); site_iter != theo_sites.end(); site_iter++)
        {
          // Check if the candidate sequence records the modified site.
          ModificationSites::iterator it1 = mod_sites.find(*site_iter);
          // Check if the modification distribution records the modified site.
          ModificationDistribution::iterator it2 = pat_iter->second.find(*site_iter);

          double score1 = (it1 != mod_sites.end() ? 1.0 : 0.0);
          double score2 = (it2 != pat_iter->second.end() ? it2->second : 0.0);

          total_cost += pow(abs(score1 - score2), 2.0);
        }

      } else {
        throw std::runtime_error("Modification type not found!");
      }
    }

    return sqrt(total_cost);
  }

  void SequencePrediction::correctSulfateDistribution()
  {
    //bool applied = param.getParameter<bool>("biosynthetic_rule").first;
    //if(!applied) {
    //  std::cout << "Not applied!\n";
    //  return;
    //}

    //// 1. Get the sulfate distribution.
    //ModificationPattern::iterator sulfate_iter = mod_pattern.find("SO3");
    //if(sulfate_iter == mod_pattern.end())
    //  throw std::runtime_error("SequencePrediction.cpp: No sulfation recorded!");

    //// Format:
    //// macro_position(branch_id, mono_id) =>
    ////    site_id => likelihood
    //std::map<MacroPosition, std::map<size_t, ModificationDistribution::iterator>> residue_sulfate; 

    //for(ModificationDistribution::iterator site_iter = sulfate_iter->second.begin(); site_iter != sulfate_iter->second.end(); site_iter++)
    //{
    //  //auto residue_iter = residue_sulfate.find(site_iter->first.macro_pos);
    //  /*if(residue_iter == residue_sulfate.end())
    //    residue_iter.insert(std::make_pair(site_iter->first.macro_pos, std::map<size_t, ModificationDistribution::iterator>()));*/
    //    
    //  residue_sulfate[site_iter->first.macro_pos].insert(std::make_pair(site_iter->first.site_id, site_iter));

    //}

    //// Sequentially check the data structure.
    //for(auto seq_iter = residue_sulfate.begin(); seq_iter != residue_sulfate.end(); seq_iter++)
    //{
    //  // Skip GlcA.
    //  if(seq_iter->second.size() == 1) continue;

    //  // 3-O cannot be higher than 6-O and 2-N.
    //  // 6-O cannot be higher than 2-N
    //  auto second_iter = seq_iter->second.find(2);
    //  auto six_iter = seq_iter->second.find(6);
    //  auto third_iter = seq_iter->second.find(3);

    //  if(second_iter == seq_iter->second.end() || six_iter == seq_iter->second.end() || third_iter == seq_iter->second.end())
    //    throw std::runtime_error("CorrectSulfateDistribution: Missing modification!");

    //  if(third_iter->second->second > six_iter->second->second) {
    //    // Swap 3-O and 6-O
    //    swap(third_iter->second->second, six_iter->second->second);
    //  }
    //  if(six_iter->second->second > second_iter->second->second) {
    //    // Swap 2-N and 6-O
    //    swap(six_iter->second->second, second_iter->second->second);
    //  }
    //}
  }

  double SequencePrediction::calculateMergedCost( const ModificationSequence& seq )
  {
    double cost = 0.0;

    // Get the mod symbol set.
    std::set<std::string> mod_symbols = gs->getModificationTypes();
    for(std::set<std::string>::iterator iter = mod_symbols.begin(); iter != mod_symbols.end(); iter++)
    {
      // Theoretical modification sites. Organize the sites in terms of their mono id.
      ModificationSites theo_sites = gs->getModificationSitesBySymbol(*iter);

      std::multimap<MacroPosition, ModificationPosition> merged_sites;
      std::set<MacroPosition> key_set;
      for(ModificationSites::iterator site_iter = theo_sites.begin(); site_iter != theo_sites.end(); site_iter++) {
        merged_sites.insert(std::make_pair(site_iter->macro_pos, *site_iter));
        key_set.insert(site_iter->macro_pos);
      }

      ModificationSequence::const_iterator seq_iter = seq.find(*iter);

      // Modification sites from the given sequence.
      ModificationSites mod_sites = seq_iter->second;

      // Recorded modification distribution for specified mod symbol.
      ModificationPattern::iterator pat_iter = mod_pattern.find(*iter);
      ModificationDistribution recorded_dist;
      if(pat_iter != mod_pattern.end())
        recorded_dist = mod_pattern[*iter];

      // Iterate over all residues.
      for(std::set<MacroPosition>::iterator mp_iter = key_set.begin(); mp_iter != key_set.end(); mp_iter++)
      {
        ModificationSites temp_sites;
        std::pair<std::multimap<MacroPosition, ModificationPosition>::iterator, std::multimap<MacroPosition, ModificationPosition>::iterator> p = merged_sites.equal_range(*mp_iter);

        double score1 = 0.0; double score2 = 0.0;
        while(p.first != p.second)
        {
          // Check if mod_sites and recorded_dist has the corresponding information, and calculate the distance.
          ModificationSites::iterator it1 = mod_sites.find(p.first->second);
          ModificationDistribution::iterator it2 = recorded_dist.find(p.first->second);

          if(it1 != mod_sites.end()) // Modified.
            score1 += 1.0;

          if(it2 != recorded_dist.end()) 
            score2 += (it2->second > 1.0 ? 1.0 : it2->second);

          ++p.first;
        }

#ifdef _DEBUG
        std::cout << "Mono ID: " << mp_iter->mono_id << " ";
        std::cout << "Score 1: " << score1 << "\tScore 2: " << score2 << "\n"; 
#endif
        cost += pow(abs(score1 - score2), 2.0);
      }

    }
    return sqrt(cost);
  }

  void printModificationDistribution( const ModificationDistribution& dist, const std::string& mod_name )
  {
      for(ModificationDistribution::const_iterator iter = dist.begin(); 
        iter != dist.end(); iter++)
      {
        std::cout << iter->first.printString() << "\t" << mod_name << "\t" << iter->second << "\n";
      }

  }

}

