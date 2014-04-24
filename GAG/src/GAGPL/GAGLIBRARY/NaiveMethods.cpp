/********************************************************************
	created:	2014/01/15
	created:	15:1:2014   23:02
	filename: 	NaiveMethods.cpp
	file path:	GAGPL\GAGLIBRARY
	file base:	NaiveMethods
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include <GAGPL/GAGLIBRARY/NaiveMethods.h>
#include <GAGPL/MATH/MassConversion.h>

namespace gag
{
  double NaiveMethods::calculateGoldenPairNum(const ModificationSequence& mod_seq)
  {
    std::multimap<MonoPeakPtr, NodeItem> filtered_data = this->filterData(mod_seq);

    int chain_length = gs->getBranchByID(0).getUnitNum();

    // The key consists of a vector of <mod_sites, mod_num> pairs.
    ModificationKeySet key_set;
    std::set<std::string> symbol_set = gs->getModificationTypes();

    int count = 0;

#ifdef _DEBUG
    std::vector<std::pair<ModificationKey, ModificationKey>> key_pool;
#endif // _DEBUG
    for(std::multimap<MonoPeakPtr, NodeItem>::iterator iter = filtered_data.begin(); iter != filtered_data.end(); iter++)
    {
      // Consider only terminal cleavage.
      if (iter->second.getCleavageNum() > 1) continue;

      const CleavageCollection& cc = iter->second.getFragment()->getCleavages();
      if(cc.size() == 0) continue;

      //std::string type = cc.begin()->first;

      // Consider both glycosidic-bond and cross-ring cleavages.
      //if(type == "A" || type == "X") continue;
#ifdef _DEBUG
      iter->second.getFragment()->printFragment();
#endif // _DEBUG

      ModificationKey mod_key, complement_key;

      double pass = true;
      // Store the results into key_vec.
      BOOST_FOREACH(std::string symbol, symbol_set)
      {
        int total_num = gs->getModificationConstraint(symbol);
        int mod_num = iter->second.getModificationNum(symbol);

        // Modification sites covered by the fragment.
        ModificationSites mod_sites = iter->second.getFragment()->getModificationSitesBySymbol(symbol, 1);
        
        // Should be all possible sites including the Ac sites.
        ModificationSites total_sites = gs->getModificationSitesBySymbol(symbol, 1);

#ifdef _DEBUG
        std::cout << "Symbol: " << symbol << "\n";
        std::cout << "Mod number: " << mod_num << "\n";
        printModificationSites(mod_sites);
#endif // _DEBUG
        // Skip meaningless fragments.
        if(mod_sites.size() == total_sites.size() || mod_sites.size() == 0) {
#ifdef _DEBUG
          std::cout << "Unqualified modification sites.  Pass the fragment!\n";
#endif // _DEBUG
          pass = false;
          break;
        }

        mod_key.addKey(mod_sites, mod_num);
        complement_key.addKey(getSiteDifference(total_sites, mod_sites), total_num - mod_num);      

      }
#ifdef _DEBUG
      std::cout << "Key!\n";
      mod_key.print();
      std::cout << "\n";
      std::cout << "Complement key!\n";
      complement_key.print();
      std::cout << "\n";
#endif // _DEBUG

      if(!pass) continue;

      // 1. Determine if the key has been recorded.
      if(key_set.find(mod_key) == key_set.end()) {
#ifdef _DEBUG
        std::cout << "Key not recorded!\n";
#endif // _DEBUG
        // 2. If not, determine if there is a complementary one recorded.
        if(key_set.find(complement_key) == key_set.end()) {
#ifdef _DEBUG
          std::cout << "Comp key not recorded!\n";
#endif // _DEBUG
          // 3. If not, store the key.
          
        } else {
          // This is a match!
#ifdef _DEBUG
          std::cout << "There is a match!\n";
          key_pool.push_back(std::make_pair(complement_key, mod_key));
#endif // _DEBUG
          count++;
        }
        key_set.insert(mod_key);

      } else {
#ifdef _DEBUG
        std::cout << "Key recorded, skip it!\n";
#endif // _DEBUG
        continue;
      } 

    }
// Check the key pair and see if there is overlapping.
#ifdef _DEBUG
    for(auto iter = key_pool.begin(); iter != key_pool.end(); iter++)
    {
      std::cout << iter->first.mod_key << "\t" << iter->second.mod_key << "\n";
    }
#endif // _DEBUG
    return (double)count;
  }

  double NaiveMethods::calculateCoverage(const ModificationSequence& mod_seq, int cleavage_num /*= 2*/ )
  {
    std::multimap<MonoPeakPtr, NodeItem> filtered_data = this->filterData(mod_seq);

    std::multimap<MonoPeakPtr, NodeItem>::iterator iter = filtered_data.begin();

    std::set<double> mass_set;
    for(; iter != filtered_data.end(); iter++)
    {
      if(iter->second.getCleavageNum() > (size_t)cleavage_num)
        continue;

      double mass = msmath::calculateMass(iter->first->mz, -1 * iter->first->z);
      mass_set.insert(mass);
    }

    return (double)mass_set.size();
  }

  std::multimap<MonoPeakPtr, NodeItem> NaiveMethods::filterData( const ModificationSequence& mod_seq )
  {
    std::multimap<MonoPeakPtr, NodeItem>::iterator iter = data.begin();
    std::multimap<MonoPeakPtr, NodeItem> new_data;

    // Get the modification types to check.
    std::set<std::string> mod_set = gs->getModificationTypes();
    //for(ModificationSequence::const_iterator const_iter = seq.begin(); const_iter != seq.end(); const_iter++)
    //mod_set.insert(const_iter->first);

    for(; iter != data.end(); iter++)
    {
#ifdef _DEBUG
      std::cout << "Fragment information:\n";
      iter->second.getFragment()->printFragment();
#endif // _DEBUG

      // Do not count internal cleavages.
      if(iter->second.getCleavageNum() == 2) {
#ifdef _DEBUG
        std::cout << "Unqualified cleavage number." << "\n\n";
#endif // _DEBUG
        continue;

      }

      // Decide if the number of modification meets with the sequence.
      bool pass = true;

      // Iterate over all modification types.
      for(std::set<std::string>::iterator mod_iter = mod_set.begin(); mod_iter!=mod_set.end(); mod_iter++) {
        // The modification number contained by the fragment.
        int mod_num = iter->second.getModificationNum(*mod_iter);

        // Check intersection of modification sites from seq and frag. 
        ModificationSequence::const_iterator seq_iter = mod_seq.find(*mod_iter);

        // Available sites of the candidate sequence.
        ModificationSites mod_sites = seq_iter->second;

        // Available sites at this region.
        //ModificationSites region_sites = iter->second.getFragment()->getModificationSitesBySymbol(*mod_iter, 1);
        ModificationSites region_sites = iter->second.getFragment()->getModificationSitesBySymbol(*mod_iter, 1);

        ModificationSites region_mod_sites = getSiteIntersection(mod_sites, region_sites);

        int max_mod_num = (int)region_mod_sites.size();

#ifdef _DEBUG
        std::cout << "Symbol: " << *mod_iter << "\tNumber: " << mod_num << "\tMax: " << max_mod_num << "\n";
        std::cout << "Modification sites from specified sequence:\n";
        printModificationSites(mod_sites);
        std::cout << "Fragment sites:\n";
        printModificationSites(region_sites);
        std::cout << "Intersection:\n";
        printModificationSites(region_mod_sites);
#endif // _DEBUG

        // If status == true, need to specify a range of the modification number.
        if(*mod_iter == "Ac") {// No missing 
          if(mod_num != max_mod_num) {
#ifdef _DEBUG
            std::cout << "Ac number does not match!\n";
#endif // _DEBUG
            pass = false;
            break;
          }
        } else if(*mod_iter == "SO3") { // Missing.

          if(status) { // Calculate based on residue level.
            // Adjust the upper_limit of sulfate number.
            std::string clv_type = iter->second.getFragment()->getFragmentType();
            if(clv_type == "A" || clv_type == "X") {
              ModificationSites expanded_sites = iter->second.getFragment()->getExpandedModificationSites(*mod_iter, 1);
              ModificationSites reduced_sites = iter->second.getFragment()->getReducedModificationSites(*mod_iter, 1);
              // The cross-ring part of the residue.
              ModificationSites cr_residue_sites = getSiteDifference(region_sites, reduced_sites);
              // The whole sites of the residue.
              ModificationSites residue_diff_sites = getSiteDifference(expanded_sites, reduced_sites);

              // Lower_bound.
              ModificationSites expanded_mod_sites1 = getSiteIntersection(mod_sites, reduced_sites);
              // Modification sites on the boundary residue.
              ModificationSites expanded_mod_sites2 = getSiteIntersection(mod_sites, residue_diff_sites);
              max_mod_num = (int)expanded_mod_sites1.size() + (int)std::min(cr_residue_sites.size(), expanded_mod_sites2.size());

#ifdef _DEBUG
              std::cout << "Fragment type: " << clv_type << "\n";
              std::cout << "Expanded fragment sites:\n";
              printModificationSites(expanded_sites);
              std::cout << "Lower_bound modified sites:\n";
              printModificationSites(expanded_mod_sites1);
              std::cout << "Cross-ring part modification sites:\n";
              printModificationSites(cr_residue_sites);
              std::cout << "Boundary Modified sites:\n";
              printModificationSites(expanded_mod_sites2);
              std::cout << "Updated max number: " << max_mod_num << "\n";
#endif // _DEBUG
            }
          } 
          if(mod_num > max_mod_num) {
#ifdef _DEBUG
            std::cout << "SO3 number beyond boundary " << max_mod_num << "\n";
#endif // _DEBUG
            pass = false;
            break;
          }
        }

      }

      if(pass) {
#ifdef _DEBUG
        std::cout << "Accept new data!\n\n";
#endif // _DEBUG

        new_data.insert(*iter);	
      } else {
#ifdef _DEBUG
        std::cout << "Reject new data!\n\n";
#endif // _DEBUG
      }
    }
#ifdef _DEBUG
    std::cout << "Old data size: " << data.size() << "\n";
    std::cout << "New data size: " << new_data.size() << "\n"; 
#endif // _DEBUG
    return new_data;
  }
}