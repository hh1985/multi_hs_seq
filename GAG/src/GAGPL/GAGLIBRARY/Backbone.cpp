#include <GAGPL/GAGLIBRARY/Backbone.h>
#include <iterator>
#include <boost/bimap/support/lambda.hpp>

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

  bool Backbone::addFamily( BackbonePtr child_node, BackbonePtr parent_node )
  {
#ifdef _DEBUG
    cout << "Before:" << _neighbors.size() << "\n";
#endif // _DEBUG
    size_t before_size = _neighbors.size();
    _neighbors.insert(NeighborType(child_node, parent_node));
    size_t after_size = _neighbors.size();
#ifdef _DEBUG
    cout << "After:" << _neighbors.size() << "\n";
#endif // _DEBUG
    return before_size != after_size;
    
    
  }

  set<BackbonePtr> Backbone::getParents()
  {
    set<BackbonePtr> parent_set;
    for(auto iter = _neighbors.by<parent>().begin(); iter != _neighbors.by<parent>().end(); iter++) {
      if(iter->get<parent>() == nullptr) continue;
      parent_set.insert(iter->get<parent>());
    }
    
    return parent_set;
  }

  const set<BackbonePtr> Backbone::getParents() const
  {
    set<BackbonePtr> parent_set;
    for(auto iter = _neighbors.by<parent>().begin(); iter != _neighbors.by<parent>().end(); iter++) {
      if(iter->get<parent>() == nullptr) continue;
      parent_set.insert(iter->get<parent>());
      
    }

    return parent_set;
  }

  set<BackbonePtr> Backbone::getChildren()
  {
    set<BackbonePtr> child_set;
    for(auto iter = _neighbors.by<child>().begin(); iter != _neighbors.by<child>().end(); iter++) {
      if(iter->get<child>() == nullptr) continue;
      child_set.insert(iter->get<child>());
    }

    return child_set;
  }

  const set<BackbonePtr> Backbone::getChildren() const
  {
    set<BackbonePtr> child_set;
    for(auto iter = _neighbors.by<child>().begin(); iter != _neighbors.by<child>().end(); iter++) {
      if(iter->get<child>() == nullptr) continue;
      child_set.insert(iter->get<child>());
    }

    return child_set;
  }

  void Backbone::addParent(BackbonePtr node)
  {
    set<BackbonePtr> child_set = this->getChildren();
    if(child_set.size() == 0) {
      //ModificationSites mod_sites;
      //BackbonePtr b_p = boost::make_shared<Backbone>(mod_sites, -1);
      this->addFamily(BackbonePtr(), node);
    } else {
      for(auto iter = child_set.begin(); iter != child_set.end(); iter++)
        this->addFamily(*iter, node);
    }
    
  }

  void Backbone::addChild(BackbonePtr node)
  {
    set<BackbonePtr> parent_set = this->getParents();
    if(parent_set.size() == 0)
      this->addFamily(node, BackbonePtr());
    else {
      for(auto iter = parent_set.begin(); iter != parent_set.end(); iter++)
        this->addFamily(node, *iter);
    } 
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
#ifdef _DEBUG
    cout << "Distance:" << std::distance(p.first, p.second) << "\n";
    for(auto iter = p.first; iter != p.second; ) {
      cout << &(*iter) << "\t" << &(*(p.second));
      iter++;
      cout << &(*iter) << "\t" << &(*(p.second));
    }
#endif // _DEBUG
    

    while(p.first != p.second)
    {
#ifdef _DEBUG

      cout << &(*(p.first)) << "\t" << &(*(p.second)) << "\n";
#endif // _DEBUG
      _neighbors.by<child>().modify_key(p.first++, _key = next_node);
#ifdef _DEBUG
      cout << &(*(p.first)) << "\t" << &(*(p.second)) << "\n";
#endif // _DEBUG

      //p.first++;

#ifdef _DEBUG
      //cout << &(*(p.first)) << "\t" << &(*(p.second)) << "\n";
      
      set<BackbonePtr> parents = this->getParents();
      cout << "Change of parents:\n";
      for(auto iter = parents.begin(); iter != parents.end(); iter++)
        cout << **iter << "\n";

      set<BackbonePtr> children = this->getChildren();
      cout << "Change of children:\n";
      for(auto iter = children.begin(); iter != children.end(); iter++)
        cout << **iter << "\n";

#endif // _DEBUG
    }
  }

  ostream& operator<<(ostream& os, const Backbone& bone)
  {
      os << bone.mod_sites << "\n";
      const set<BackbonePtr> parent_set = bone.getParents();
      const set<BackbonePtr> child_set = bone.getChildren();

      for(auto p_iter = parent_set.begin(); p_iter != parent_set.end(); p_iter++) {
        if(*p_iter == nullptr) continue;

        os << "--Parent:" << (*p_iter)->mod_sites << "\n";
      }

      for(auto c_iter = child_set.begin(); c_iter != child_set.end(); c_iter++) {
        if(*c_iter == nullptr) continue;

        os << "--Child:" << (*c_iter)->mod_sites << "\n";
      }

      return os;
  }
}