#include <GAGPL/GAGLIBRARY/Backbone.h>
#include <iterator>
#include <boost/bimap/support/lambda.hpp>

namespace gag
{

//  void Backbone::addAssignment( AssignmentPtr assignment, const string& mod_symbol )
//  {
//    // Check if the assignment is qualified.
//      
//      if(mod_sites.size() == 0)
//          mod_sites = assignment->getBackboneModificationSites(mod_symbol);
//      else if(!(mod_sites == assignment->getBackboneModificationSites(mod_symbol))) {
//
//#ifdef _DEBUG  
//        std::cout << "Assignment doesn't match!\n";
//#endif
//        return;
//      }
//
//    int mod_num = assignment->getModificationNumber(mod_symbol);
//    members[mod_num].insert(assignment);
//
//  }

  /*set<AssignmentPtr> Backbone::getAssignmentsByModNumber( int mod_num )
  {
    if(members.find(mod_num) != members.end())
      return members[mod_num];
    else
      return set<AssignmentPtr>();
  }*/

  //set<int> Backbone::getModNumbers() const
  //{
  //  set<int> mod_num_set;
  //  for(auto iter = members.begin(); iter != members.end(); iter++)
  //  {
  //    mod_num_set.insert(iter->first);
  //  }
  //  return mod_num_set;
  //}

  //int Backbone::getLargestModNumber() const
  //{
  //  return members.size() == 0 ? 0 : members.rbegin()->first;
  //}

//  bool Backbone::addFamily( BackbonePtr child_node, BackbonePtr parent_node )
//  {
//#ifdef _DEBUG
//    cout << "Before:" << _neighbors.size() << "\n";
//#endif // _DEBUG
//    size_t before_size = _neighbors.size();
//    _neighbors.insert(NeighborType(child_node, parent_node));
//    size_t after_size = _neighbors.size();
//#ifdef _DEBUG
//    cout << "After:" << _neighbors.size() << "\n";
//#endif // _DEBUG
//    return before_size != after_size;
//    
//    
//  }

  set<BackbonePtr>& Backbone::getParents()
  {
    return _parents;
  }

  set<BackbonePtr>& Backbone::getChildren()
  {
    return _children;
  }

  void Backbone::addParent(BackbonePtr node)
  {
    _parents.insert(node); 
  }

  void Backbone::addChild(BackbonePtr node)
  {
    _children.insert(node);
  }

  void Backbone::replaceParent(BackbonePtr last_node, BackbonePtr next_node)
  {
    auto iter = _parents.find(last_node);
    if(iter != _parents.end()) {
      _parents.erase(iter);
      
      next_node->addChild(this->self());
      this->addParent(next_node);
      //_parents.insert(next_node);
    }

  }

  void Backbone::replaceChild(BackbonePtr last_node, BackbonePtr next_node)
  {
    auto iter = _children.find(last_node);
    if(iter != _children.end()) {
      _children.erase(iter);

      next_node->addParent(this->self());
      this->addChild(next_node);
    }
  }

  void Backbone::addAssignment( AssignmentPtr assignment )
  {
    members.insert(assignment);
  }

  void Backbone::removeAssignment( AssignmentPtr assignment )
  {
    members.erase(assignment);
  }

  void Backbone::removeParent( BackbonePtr node )
  {
    _parents.erase(node); // Exception?
  }

  void Backbone::removeChild( BackbonePtr node )
  {
    _children.erase(node);
  }


  ostream& operator<<(ostream& os, const Backbone& bone)
  {
      // Child on the left and parent on the right.
      for(auto iter = bone._neighbors.begin(); iter != bone._neighbors.end(); iter++)
      {

        if(iter->get<child>() == nullptr)
          os << "NULL\t";
        else
          os << iter->get<child>()->mod_sites << "\t" << bone.mod_sites << "\n";

        os << bone.mod_sites << "\t";
        if(iter->get<parent>() == nullptr)
          os << "NULL\n";
        else
          os << iter->get<parent>()->mod_sites << "\n";

      }

      return os;
  }
}