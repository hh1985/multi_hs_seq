#include <GAGPL/GAGLIBRARY/Backbone.h>

namespace gag
{

  void Backbone::addAssignment( AssignmentPtr assignment, const string& mod_symbol )
  {
    // Check if the assignment is qualified.
    if(mod_sites.size() == 0 || mod_sites == assignment->getBackboneModificationSites(mod_symbol)) {
      int mod_num = assignment->getModificationNumber(mod_symbol);
      members[mod_num].insert(assignment);
    } else {
      cout << "Assignment doesn't match!\n";
    }

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

  ostream& operator<<(ostream& os, const Backbone& bone)
  {
      os << bone.mod_sites;

      for(auto iter = bone.members.begin(); iter != bone.members.end(); iter++)
      {
          os << "Number: " << iter->first << "\n";
          for(auto assign_iter = (*iter).second.begin(); assign_iter != (*iter).second.end(); assign_iter++)
          {
              os << *assign_iter << "\n";
          }
      }

      return os;
  }
}