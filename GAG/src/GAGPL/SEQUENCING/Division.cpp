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
    os << "Support assignments: " << div.getSupportAssignments().size() << "\n";

#ifdef _DEBUG

    for(auto iter = div.getSupportAssignments().begin(); iter != div.getSupportAssignments().end(); iter++)
      os << **iter << "\n";

#endif // _DEBUG

    // 3. Environment information (included in DivMap)
    return os;
  }


  void Division::addParent(DivisionPtr p_div)
  {
    // Check the parents.
    this->updateParent(p_div);
    p_div->updateChild(this->self());
  }

  void Division::addChild(DivisionPtr c_div)
  {
    this->updateChild(c_div);
    c_div->updateParent(this->self());
  }

  void Division::replaceParent(DivisionPtr old_div, DivisionPtr new_div)
  {
    // 1. Check if the parent exist.
    auto iter = parent_set.find(old_div);
    if(iter != parent_set.end()) {
      parent_set.erase(old_div);
      parent_set.insert(new_div);
    } else {
      // Not found: simply add the new div_node.
      parent_set.insert(new_div);
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
      child_set.insert(new_div);
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

  bool Division::isLargerThan(DivisionPtr div_ptr) const
  {
    return containSubset(mod_sites, div_ptr->getModificationSites()) && ((int)mod_sites.size() > div_ptr->getModificationSitesNumber());
  }

  bool Division::isSmallerThan(DivisionPtr div_ptr) const
  {
    return containSubset(div_ptr->getModificationSites(), mod_sites) && ((int)mod_sites.size() < div_ptr->getModificationSitesNumber());
  }

  bool Division::isSibling(DivisionPtr div_ptr) const
  {
    return !this->isLargerThan(div_ptr) && !this->isSmallerThan(div_ptr);
  }

  void Division::addAssignment(AssignmentPtr assign)
  {
    string assign_type = assign->getCleavageType();
    if(_type.empty()) {
      _type = assign_type;
    } else if(_type != "G" && assign_type == "G"){
      _type = assign_type;
    } else {
      // Do nothing.
    }

    assign_support.insert(assign);
  }

  int Division::getInterSitesNumber(const DivisionPtr div_node) const
  {
    return (int)getSiteIntersection(mod_sites, div_node->getModificationSites()).size();
  }

  int Division::getDiffSitesNumber(const DivisionPtr div_node) const
  {
    return (int)getSiteDifference(mod_sites, div_node->getModificationSites()).size();
  }

  void Division::insertNode(DivisionPtr c_div, DivisionPtr p_div)
  {
    parent_set.insert(p_div);
    child_set.insert(c_div);
    c_div->replaceParent(p_div, this->self());
    p_div->replaceChild(c_div, this->self());
  }

  bool Division::isCompatible(DivisionPtr div_ptr) const
  {
    int num1 = this->getModificationNumber();
    int num2 = div_ptr->getModificationNumber();
    int lb, hb;
    if(this->isSmallerThan(div_ptr)) {
      lb = num1; 
      hb = lb + (div_ptr->getModificationSitesNumber() - this->getModificationSitesNumber());
    } else if(this->isLargerThan(div_ptr)) {
      hb = num1;
      lb = max(0, num1-abs(div_ptr->getModificationSitesNumber() - this->getModificationSitesNumber()));

    } else {
      // Sibling.
      ModificationSites cur_sites = this->getModificationSites();
      ModificationSites div_sites = div_ptr->getModificationSites();

      int diff_ab = (int)getSiteDifference(cur_sites, div_sites).size();

      int diff_ba = (int)getSiteDifference(div_sites, cur_sites).size();
      int int_ab = this->getInterSitesNumber(div_ptr);

      lb = num1 - diff_ab;
      hb = min(int_ab, num1) + diff_ba;
    }

    return div_ptr->getModificationNumber() >= lb && div_ptr->getModificationNumber() <= hb;
  }

  void Division::insertParent(DivisionPtr p_div)
  {
    parent_set.insert(p_div);
  }

  void Division::insertChild(DivisionPtr c_div)
  {
    child_set.insert(c_div);
  }

  void Division::updateParent(DivisionPtr new_ptr)
  {
    set<DivisionPtr> check_set;
    for(auto iter = parent_set.begin(); iter != parent_set.end(); iter++)
    {
      if((new_ptr)->isSmallerThan(*iter) && new_ptr->isCompatible(*iter))
        check_set.insert(*iter);
    }
    for(auto iter = check_set.begin(); iter != check_set.end(); iter++)
    {
      this->replaceParent(*iter, new_ptr);
      (*iter)->replaceChild(this->self(), new_ptr);
    }

    if(check_set.size() == 0) this->insertParent(new_ptr);
  }

  void Division::updateChild(DivisionPtr new_ptr)
  {
    set<DivisionPtr> check_set;
   
    for(auto iter = child_set.begin(); iter != child_set.end(); iter++)
    {
      if((new_ptr)->isLargerThan(*iter) && new_ptr->isCompatible(*iter))
        check_set.insert(*iter);
    }
    
    for(auto iter = check_set.begin(); iter != check_set.end(); iter++)
      this->replaceChild(*iter, new_ptr);

    if(check_set.size() == 0) this->insertChild(new_ptr);
   
  }

  void Division::addNode(DivisionPtr c_div, DivisionPtr p_div)
  {
    this->addChild(c_div);
    this->addParent(p_div);
  }



}

