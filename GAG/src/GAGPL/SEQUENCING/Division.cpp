#include <GAGPL/SEQUENCING/Division.h>

namespace gag
{
  
  ostream& operator<<(ostream& os, const Division& div)
  {
    // Print the modification sites and modification number.
    // 1. Structure information.
    os << "\n";
    os << "General type: " << div.getGeneralType() << "\n";
    os << "Modification sites: " << div.getModificationSites() << "\n";
    os << "Modification number: " << div.getModificationNumber() << "\tCorrection number: " << div.getCorrectionNumber() << "\n";

    // 2. Support information.
    os << "Support assignments:\n";
    for(auto iter = div.getSupportAssignments().begin(); iter != div.getSupportAssignments().end(); iter++)
    {
      os << *iter;
    }

    // 3. Environment information (included in DivMap)
    return os;
  }


  void Division::addParent(DivisionPtr p_div)
  {
    parent_set.insert(p_div);
  }

  void Division::addChild(DivisionPtr c_div)
  {
    child_set.insert(c_div);
  }

  void Division::replaceParent(DivisionPtr old_div, DivisionPtr new_div)
  {
    // 1. Check if the parent exist.
    auto iter = parent_set.find(old_div);
    if(iter != parent_set.end()) {
      parent_set.erase(old_div);
      parent_set.insert(new_div);
    } else {
#ifdef _DEBUG
      cout << "Parent not found!\n";
#endif // _DEBUG
    }
  }

  void Division::replaceChild(DivisionPtr old_div, DivisionPtr new_div)
  {
    // 1. Check if the parent exist.
    auto iter = child_set.find(old_div);
    if(iter != child_set.end()) {
      child_set.erase(old_div);
      child_set.insert(new_div);
    } else {
#ifdef _DEBUG
      cout << "Child not found!\n";
#endif // _DEBUG
    }
  }

  set<DivisionPtr>& Division::getParents()
  {
    return parent_set;
  }

  set<DivisionPtr>& Division::getChildren()
  {
    return child_set;
  }

  int Division::getInDegree() const
  {
    return (int)child_set.size();
  }

  int Division::getOutDegree() const
  {
    return (int)parent_set.size();
  }

  bool Division::isLargerThan(DivisionPtr div_ptr)
  {
    return containSubset(mod_sites, div_ptr->getModificationSites()) && ((int)mod_sites.size() > div_ptr->getModificationSitesNumber());
  }

  bool Division::isSmallerThan(DivisionPtr div_ptr)
  {
    return containSubset(div_ptr->getModificationSites(), mod_sites) && ((int)mod_sites.size() < div_ptr->getModificationSitesNumber());
  }

  bool Division::isSibling(DivisionPtr div_ptr)
  {
    return !this->isLargerThan(div_ptr) && !this->isSmallerThan(div_ptr);
  }

  void Division::addAssignment(AssignmentPtr assign)
  {
    assign_support.insert(assign);
  }

  int Division::getInterSitesNumber(DivisionPtr div_node) const
  {
    return (int)getSiteIntersection(mod_sites, div_node->getModificationSites()).size();
  }

  int Division::getDiffSitesNumber(DivisionPtr div_node) const
  {
    return (int)getSiteDifference(mod_sites, div_node->getModificationSites()).size();
  }

}

