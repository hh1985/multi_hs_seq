#include <GAGPL/GAGLIBRARY/Backbone.h>

namespace gag
{

  void Backbone::addAssignment( AssignmentPtr assignment, const string& mod_symbol )
  {
    // Check if the assignment is qualified.
      
      if(mod_sites.size() == 0)
          mod_sites = assignment->getBackboneModificationSites(mod_symbol);
      else if(!(mod_sites == assignment->getBackboneModificationSites(mod_symbol))) {

#ifdef _DEBUG  
        std::cout << "Assignment doesn't match!\n";
#endif
        return;
      }

    int mod_num = assignment->getModificationNumber(mod_symbol);
    members[mod_num].insert(assignment);

  }

  set<AssignmentPtr> Backbone::getAssignmentsByModNumber( int mod_num )
  {
    if(members.find(mod_num) != members.end())
      return members[mod_num];
    else
      return set<AssignmentPtr>();
  }

  set<int> Backbone::getModNumbers() const
  {
    set<int> mod_num_set;
    for(auto iter = members.begin(); iter != members.end(); iter++)
    {
      mod_num_set.insert(iter->first);
    }
    return mod_num_set;
  }

  int Backbone::getLargestModNumber() const
  {
    return members.size() == 0 ? 0 : members.rbegin()->first;
  }

  void Backbone::addFamily(BackbonePtr child_node, BackbonePtr parent_node)
  {
    _neighbors.insert(NeighborType(child_node, parent_node));
  }

  set<BackbonePtr> Backbone::getParents()
  {
    set<BackbonePtr> parent_set;
    for(auto iter = _neighbors.by<parent>().begin(); iter != _neighbors.by<parent>().end(); iter++)
      parent_set.insert(iter->get<parent>());
    
    return parent_set;
  }

  const set<BackbonePtr> Backbone::getParents() const
  {
    set<BackbonePtr> parent_set;
    for(auto iter = _neighbors.by<parent>().begin(); iter != _neighbors.by<parent>().end(); iter++)
      parent_set.insert(iter->get<parent>());

    return parent_set;
  }

  set<BackbonePtr> Backbone::getChildren()
  {
    set<BackbonePtr> child_set;
    for(auto iter = _neighbors.by<child>().begin(); iter != _neighbors.by<child>().end(); iter++)
      child_set.insert(iter->get<child>());

    return child_set;
  }

  const set<BackbonePtr> Backbone::getChildren() const
  {
    set<BackbonePtr> child_set;
    for(auto iter = _neighbors.by<child>().begin(); iter != _neighbors.by<child>().end(); iter++)
      child_set.insert(iter->get<child>());

    return child_set;
  }

  void Backbone::addParent(BackbonePtr node)
  {
    set<BackbonePtr> child_set = this->getChildren();
    for(auto iter = child_set.begin(); iter != child_set.end(); iter++)
      this->addFamily(*iter, node);
  }

  void Backbone::addChild(BackbonePtr node)
  {
    set<BackbonePtr> parent_set = this->getParents();
    for(auto iter = parent_set.begin(); iter != parent_set.end(); iter++)
      this->addFamily(node, *iter);
  }

  void Backbone::replaceParent(BackbonePtr last_node, BackbonePtr next_node)
  {
    auto p = _neighbors.by<parent>().equal_range(last_node);
    while(p.first != p.second)
    {
      _neighbors.by<parent>().replace_key(p.first, next_node);
      p.first++;
    }
  }

  void Backbone::replaceChild(BackbonePtr last_node, BackbonePtr next_node)
  {
    auto p = _neighbors.by<child>().equal_range(last_node);
    while(p.first != p.second)
    {
      _neighbors.by<child>().replace_key(p.first, next_node);
      p.first++;
    }
  }

  ostream& operator<<(ostream& os, const Backbone& bone)
  {
      os << bone.mod_sites << "\n";
      const set<BackbonePtr> parent_set = bone.getParents();
      const set<BackbonePtr> child_set = bone.getChildren();

      for(auto iter = parent_set.begin(); iter != parent_set.end(); iter++)
        os << "--Parent:" << (*iter)->mod_sites << "\n";

      for(auto iter = child_set.begin(); iter != child_set.end(); iter++)
        os << "--Child:" << (*iter)->mod_sites << "\n";

      return os;
  }
}