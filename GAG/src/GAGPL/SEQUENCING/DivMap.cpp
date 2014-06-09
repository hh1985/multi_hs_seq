#include <GAGPL/SEQUENCING/DivMap.h>

namespace gag
{

  DivMap::DivMap(GlycanSequencePtr gs, string mod_symbol)
    : full_sites(gs->getModificationSitesBySymbol(mod_symbol, 1)), full_num(gs->getModificationConstraint(mod_symbol)), mod_symbol(mod_symbol)
  {
    this->initilize();
  }

  DivMap::DivMap(const ModificationSites& mod_sites, int mod_num, string mod_symbol)
    : full_sites(mod_sites), full_num(mod_num), mod_symbol(mod_symbol)
  {
    this->initilize();
  }

  void DivMap::initilize()
  {
    // Define the dummy nodes.
    _empty_node = boost::make_shared<Division>(ModificationSites(), 0);
    _full_node = boost::make_shared<Division>(full_sites, full_num);

    _empty_node->addParent(_full_node);
    _full_node->addChild(_empty_node);
  }

  //void DivMap::addAssignment(AssignmentPtr assign)
  //{
  //  // 0. Background clearance.
  //  if(!qualityCheck(assign)) return;

  //  // 1. Get the modification sites of the assignment
  //  ModificationSites mod_sites = assign->getBackboneModificationSites(mod_symbol);
  //  int mod_num = assign->getModificationNumber(mod_symbol);
  //  
  //  // 2. Check if this type of assignment has been recorded.
  //  // TBD: Convert the type.
  //  if(1) {// If the assignment is RE.
  //    mod_sites = getSiteDifference(full_sites, mod_sites);
  //    mod_num = full_num - mod_num;
  //  }

  //  DivisionPtr check_div = this->checkRecord(mod_sites, mod_num);
  //  if(check_div == nullptr) {
  //    check_div = boost::make_shared<Division>(mod_sites, mod_num);
  //    // Locate the div in the network.
  //    this->addDivisionNode(check_div);
  //  } else {
  //    // Insert the assignment into the Division object.
  //    check_div->addAssignment(assign);
  //  }
  //}

  void DivMap::addDivisionNode(DivisionPtr div_node)
  {
    this->exploreMap(_empty_node, _empty_node, div_node);
  }

  void DivMap::exploreMap(DivisionPtr child_node, DivisionPtr check_node, DivisionPtr div_node)
  {
    set<DivisionPtr> parents = check_node->getParents();
    for(auto iter = parents.begin(); iter != parents.end(); iter++)
    {
      if(!checkCompatibility(div_node, *iter)) {

        // Not compatible.  Bypass the checked node.
        this->exploreMap(child_node, *iter, div_node);

      } else {

        if(div_node->isLargerThan(*iter)) {

            // Choose a new child node.
            this->exploreMap(*iter, *iter, div_node);

        } else if(div_node->isSmallerThan(*iter)) {
            
          div_node->addParent(*iter);
          div_node->addChild(child_node);

          child_node->replaceParent(*iter, div_node);
          (*iter)->replaceChild(child_node, div_node);

        } else {
          // Sibling.  Bypass the checked node.
          this->exploreMap(child_node, *iter, div_node);

        }
      }
    }
  }

  bool DivMap::checkCompatibility(DivisionPtr div1, DivisionPtr div2)
  {
    int num1 = div1->getModificationNumber();
    int num2 = div2->getModificationNumber();
    int lb, hb;
    if(div1->isSmallerThan(div2)) {
      lb = num1; 
      hb = lb + (div2->getModificationSitesNumber() - div1->getModificationSitesNumber());
    } else if(div1->isLargerThan(div2)) {
      hb = num1;
      lb = max(0, num1-abs(div2->getModificationSitesNumber() - div1->getModificationSitesNumber()));

    } else {
      // Sibling.
      int diff_ab = div1->getDiffSitesNumber(div2);
      int diff_ba = div2->getDiffSitesNumber(div1);
      int int_ab = div1->getInterSitesNumber(div2);

      lb = num1 - diff_ab;
      hb = min(int_ab, num1) + diff_ba;
    }

    return div2->getModificationNumber() >= lb && div2->getModificationNumber() <= hb;
  }

  DivisionPtr DivMap::checkRecord(const ModificationSites& sites, int num) const
  {
     auto site_iter = _map.find(sites);
     if(site_iter != _map.end()) {
       auto num_iter = site_iter->second.find(num);
       if(num_iter != site_iter->second.end()) {
         return num_iter->second;
       } else {
         return nullptr;
       }
     } else {
       return nullptr;
     }
  }

  bool DivMap::qualityCheck(AssignmentPtr assign) const
  {
     // Consider only the assignment is internal assignment.
    string type = assign->getCleavageType();
    return type == "C" || type == "G";
  }

  DivPath DivMap::getMaximumPath()
  {
    DivPath paths;

    return paths;
  }

  DivPath DivMap::getMinimumPath()
  {
    DivPath paths;

    return paths;
  }

  ostream& operator<<(ostream& os, const DivMap& div_map)
  {
    // Iterate over the map.
    for(auto mod_iter = div_map._map.begin(); mod_iter != div_map._map.end(); mod_iter++)
    {
      auto div_int = mod_iter->second;
      for(auto num_iter = div_int.begin(); num_iter != div_int.end(); num_iter++)
      {
        os << num_iter->second << "\n";
      }
    }
    return os;
  }

  ostream& operator<<(ostream& os, const DivPath& path)
  {
    for(auto iter = path.path_set.begin(); iter != path.path_set.end(); iter++)
    {
      for(auto vec_iter = (*iter).begin(); vec_iter != (*iter).end(); vec_iter++)
      {
        os << **vec_iter << "\n";
      }
    }
    return os;
  }



}
