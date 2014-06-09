#include <GAGPL/SEQUENCING/DivisionCluster.h>

namespace gag
{
  void DivisionCluster::addAssignment(AssignmentPtr assign)
  {
    ModificationSites sites = assign->getBackboneModificationSites(_symbol);
    int num = assign->getModificationNumber(_symbol);

#ifdef _DEBUG
      cout << "Assignment information:\n";
      cout << "Mod sites: " << sites << "\n";
      cout << "Mod num: " << num << "\n";
#endif // _DEBUG

    if(assign->isRECleavage()) {
      // Convert to NRE cleavage.
#ifdef _DEBUG
      cout << "Convert to NRE cleavage!\n";
#endif // _DEBUG

      sites = getSiteDifference(full_sites, sites);
      num = full_num - num;
    }
    DivisionPtr div_ptr = this->locateDivision(sites, num);

    if(div_ptr == nullptr) {
      DivisionPtr new_ptr = boost::make_shared<Division>(sites, num);
      new_ptr->addAssignment(assign);

      _backbone[sites].insert(std::make_pair(num, new_ptr));        

    } else {
      div_ptr->addAssignment(assign);
    }

    _size++;
  }

  DivisionPtr DivisionCluster::locateDivision(const ModificationSites& sites, int num)
  {
    auto site_iter = _backbone.find(sites);
    if(site_iter != _backbone.end()) {
      auto temp_map = site_iter->second;
      auto num_iter = temp_map.find(num);
      if(num_iter != temp_map.end()) {
        return num_iter->second;
      } else {
        return nullptr;
      }
    } else {
      return nullptr;
    }
  }

  size_t DivisionCluster::size() const
  {
    return _size;
  }

  set<DivisionPtr> DivisionCluster::getSelectedDivision(bool highest)
  {
    set<DivisionPtr> div_set;
    for(auto site_iter = _backbone.begin(); site_iter != _backbone.end(); site_iter++)
    {
      auto& temp_map = site_iter->second;

      // Get the Division object with the highest mod number.
      auto num_iter = temp_map.rbegin();

      if(highest) {
     
        if(this->isQualifiedDivision(num_iter->second))
          div_set.insert(num_iter->second);
     
      } else {
        
        // Keep all the divisions.
        while(num_iter != temp_map.rend()) {
         
          if(this->isQualifiedDivision(num_iter->second))
            div_set.insert(num_iter->second);
          
          num_iter++;
        }
      }

    }

    return div_set;
  }

  bool DivisionCluster::isQualifiedDivision(DivisionPtr div_ptr)
  {
    ModificationSites div_sites = div_ptr->getModificationSites();
    int div_num = div_ptr->getModificationNumber();

    // No useful information.
    if(div_sites.size() == 0 || div_sites.size() == full_sites.size()) return false;

    ModificationSites comp_sites = getSiteDifference(full_sites, div_sites);
    int comp_num = full_num - div_num;

    // Overload.
    if(div_num > (int)div_sites.size() || comp_num > (int)comp_sites.size())
      return false;

    return true;
  }

}