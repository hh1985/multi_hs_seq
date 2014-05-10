#include <GAGPL/GAGLIBRARY/FullMap.h>

namespace gag
{

  void FullMap::connectBackboneSet()
  {
#ifdef _DEBUG
    int count = 0;
    cout << "Backbone size:" << _bone_set.size() << "\n";
#endif // _DEBUG
    for(auto bone_iter = _bone_set.begin(); bone_iter != _bone_set.end(); bone_iter++)
    {
#ifdef _DEBUG
      count++;
      cout << "Move to count:" << count << "\t" << **bone_iter << "\n";
#endif // _DEBUG
      // Do not consider other full nodes.
      if((*bone_iter)->mod_sites == _full_node->mod_sites)
        continue;

      this->exploreEntryPoint(*bone_iter, _empty_node);

#ifdef _DEBUG
      cout << "\n***After insertion/appending:***\n";
      cout << "Current:" << **bone_iter << "\n";
      set<BackbonePtr> parents = (*bone_iter)->getParents();
      for(auto iter = parents.begin(); iter != parents.end(); iter++)
        cout << "--Parent:" << **iter << "\n";
      set<BackbonePtr> children = (*bone_iter)->getChildren();
      for(auto iter = children.begin(); iter != children.end(); iter++)
        cout << "--Children:" << **iter << "\n";
#endif // _DEBUG
    }
  }

  void FullMap::initialize()
  {
    _empty_node = boost::make_shared<Backbone>(ModificationSites(), 0);
    
    _full_node = boost::make_shared<Backbone>(seq->getModificationSitesBySymbol(_mod_type, 1), seq->getModificationConstraint(_mod_type));

    _empty_node->addParent(_full_node);
    _full_node->addChild(_empty_node);


#ifdef _DEBUG
    cout << "Empty:\n";
    cout << *_empty_node << "\n";
    cout << "Full:\n";
    cout << *_full_node << "\n";
#endif // _DEBUG
  }

  bool FullMap::exploreUpperNodes( BackbonePtr cur, BackbonePtr child_node, BackbonePtr parent_node )
  {

#ifdef _DEBUG
      cout << "\n******Before*****\n";
      cout << "Current:" << *cur << "\n";
      cout << "Checked child:" << *child_node << "\n";
      cout << "Checked parent:" << *parent_node << "\n";
#endif // _DEBUG
      // Get parent nodes.
      //set<BackbonePtr> parents = child_node->getParents();

      // Iterate over all parents. Decide if there is any chance to locate the position of the current node.
      bool sib = true;
      set<BackbonePtr> children = parent_node->getChildren();
      set<BackbonePtr> children_cur = cur->getChildren();
      set<BackbonePtr> parents_cur = cur->getParents();
      for(auto iter = children.begin(); iter != children.end(); iter++)
      {
#ifdef _DEBUG
        cout << "Current child:\t" << **iter << "\n";
#endif // _DEBUG

        if(*iter == cur) {
#ifdef _DEBUG
            cout << "Skip the check!\n";
#endif // _DEBUG
            if(sib)
              sib = false;

            continue;
        }

        // To the end.
        if((*iter)->isSmaller(cur)) {
#ifdef _DEBUG
          cout << "Found boundary!\n";
#endif // _DEBUG
          if(*iter == child_node) {
            // Successful appending.
#ifdef _DEBUG
            cout << "Unexpected situation\n";
#endif // _DEBUG
            child_node->addParent(cur);
            parent_node->addChild(cur);
            cur->addFamily(child_node, parent_node);
            if(sib)
              sib = false;
          } else {
            // Insertion.
#ifdef _DEBUG
            cout << "Insertion!\n";
#endif // _DEBUG
            (*iter)->replaceParent(parent_node, cur);
            parent_node->replaceChild(*iter, cur);
            cur->addFamily(*iter, parent_node);
          } 

        } else if((*iter)->isLarger(cur)) {
          // Dive deeper.
#ifdef _DEBUG
          cout << "Dive deeper!\n";
#endif // _DEBUG
          for(auto cp_iter = children.begin(); cp_iter != children.end(); cp_iter++)
          {
            this->exploreUpperNodes(cur, child_node, *cp_iter);
          }
          if(sib)
            sib = false;
        }
      }


      //if(cur->isSmaller(parent_node))
      //if(parent_node->isLarger(cur)){ // To the boundary.
      // /* child_node->replaceParent(parent_node, cur);
      //  parent_node->replaceChild(child_node, cur);*/

      //  // Append the node.
      //  //cur->addChild(child_node);
      //  //cur->addParent(parent_node);
      //  bool success_insertion = cur->addFamily(child_node, parent_node);
      //  
      //  // No need to check further information.
      //  if(success_insertion) {
      //    //child_node->addParent(cur);
      //    set<BackbonePtr> children_parent = parent_node->getChildren();
      //    //parent_node->addChild(cur);
      //    for(auto iter = children_parent.begin(); iter != children_parent.end(); iter++)
      //    {
      //      if((*iter)->isSmaller(cur)) {
      //        (*iter)->replaceParent(parent_node, cur);
      //        parent_node->replaceChild(*iter, cur);
      //        cur->addFamily(*iter, parent_node);
      //      }
      //    }
      //    
      //    child_node->addParent(cur);
      //    parent_node->addChild(cur);

      //  }

      //  if(sib)
      //    sib = false;
      //} 

      // Establish the appending case.
      if(sib) {
#ifdef _DEBUG
        cout << "Establish the appending case:\n";
#endif // _DEBUG
        child_node->addParent(cur);
        parent_node->addChild(cur);
        cur->addFamily(child_node, parent_node);
      }

      
      return sib;
      //else {
      //  set<BackbonePtr> parents_parent = parent_node->getParents();
      //  for(auto p_iter = parents_parent.begin(); p_iter != parents_parent.end(); p_iter++)
      //  {
      //    // Adjust the child node and parent node.
      //    this->exploreDeepNodes(cur, child_node, *p_iter);
      //  }

      //} 
  }

  void FullMap::exploreCompatibility( AssignmentPool& pool )
  {
   // TBD

  }

  bool FullMap::checkCompatibility( BackbonePtr small_bone, BackbonePtr large_bone )
  {
    // Maintain assignments with different mod numbers.
    bool flag = true;
    ModificationSites small_sites = small_bone->mod_sites;
    ModificationSites large_sites = large_bone->mod_sites;
    int diff_size = (int)getSiteDifference(large_sites, small_sites).size();

    if(_mod_type == "Ac") {
      
      set<int> small_num_set = small_bone->getModNumbers();
      
      set<int> large_num_set = large_bone->getModNumbers();
      
      for(auto it1 = small_num_set.begin(); it1 != small_num_set.end(); it1++)
      {
        for(auto it2 = large_num_set.begin(); it2 != large_num_set.end(); it2++)
        {
          if(this->checkCompatibility(small_bone, large_bone, *it1, *it2, diff_size) && flag)
            flag = false;
        }
      }

    } else if(_mod_type == "SO3") {
      int small_num = small_bone->getLargestModNumber();
      int large_num = large_bone->getLargestModNumber();

      flag = this->checkCompatibility(small_bone, large_bone, small_num, large_num, diff_size);
    }

    return flag;
  }

  bool FullMap::checkCompatibility( BackbonePtr small_bone, BackbonePtr large_bone, int small_num, int large_num, int diff_size )
  {
    int lower_num = small_num;
    int upper_num = diff_size + small_num;

    // Conflict!!!
    if(large_num < lower_num || large_num > upper_num)
    {
      set<AssignmentPtr> small_assign = small_bone->getAssignmentsByModNumber(small_num);
      set<AssignmentPtr> large_assign = large_bone->getAssignmentsByModNumber(large_num);
      for(auto assign_iter1 = small_assign.begin(); assign_iter1 != small_assign.end(); assign_iter1++)
        for(auto assign_iter2 = large_assign.begin(); assign_iter2 != large_assign.end(); assign_iter2++)
        {
          _status_table.insert(make_pair(*assign_iter1, *assign_iter2));
          _status_table.insert(make_pair(*assign_iter2, *assign_iter1));
        }
      
      return false;
    } else {
      return true;
    }
  }

  void FullMap::exploreEntryPoint( BackbonePtr cur, BackbonePtr child_node)
  {
#ifdef _DEBUG
    cout << "Entry Point:\n";
    cout << "Current:\n";
    cout << *cur << "\n";
    cout << "Child:\n";
    cout << *child_node << "\n";
#endif // _DEBUG

    bool sib = true;
    set<BackbonePtr> parents = child_node->getParents();
    set<BackbonePtr> children_cur = cur->getChildren();
    set<BackbonePtr> parents_cur = cur->getParents();
    for(auto iter = parents.begin(); iter != parents.end(); iter++)
    {
#ifdef _DEBUG
      cout << "Current parent:\t" << **iter << "\n";
#endif // _DEBUG
      // There is no need to check if the parent is cur itself or a recorded child of cur.
      if(*iter == cur || children_cur.find(*iter) != children_cur.end() || 
        parents_cur.find(*iter) != parents_cur.end()) {
#ifdef _DEBUG
          cout << "Skip the check!\n";
#endif // _DEBUG
          if(sib)
            sib = false;

          continue;
      }

      if((*iter)->isLarger(cur)){ // Insertion.
#ifdef _DEBUG
        cout << "For insertion!\n";
#endif // _DEBUG
        set<BackbonePtr> children_parent = (*iter)->getChildren();

        for(auto cp_iter = children_parent.begin(); cp_iter != children_parent.end(); cp_iter++)
        {
#ifdef _DEBUG
          cout << "Check the kids:" << **cp_iter << "\n";
#endif // _DEBUG
          if((*cp_iter)->isSmaller(cur)) {
            (*cp_iter)->replaceParent(*iter, cur);
            (*iter)->replaceChild(*cp_iter, cur);
            cur->addFamily(*cp_iter, *iter);
          }
        }

#ifdef _DEBUG
        cout << "Current:\n";
        cout << *cur << "\n";
        cout << "Child:\n";
        cout << *child_node << "\n";
        cout << "Parent:\n";
        cout << **iter << "\n";
#endif // _DEBUG
        
        if(sib)
          sib = false;

      } else if((*iter)->isSmaller(cur)) {   // Dive deeper for insertion.
#ifdef _DEBUG
        cout << "Dive deeper!\n";
#endif // _DEBUG
        this->exploreEntryPoint(cur, *iter);
        
        if(sib)
          sib = false;
      } 
    }

    // If all the paths are not acceptable. Creating a sibling.
    if(sib) { // Appending.
      // The child_node is fixed. Trying to find the upper bound.
      this->exploreUpperNodes(cur, child_node, _full_node);
      /*for(auto iter = parents.begin(); iter != parents.end(); iter++)
      {
        this->exploreUpperNodes(cur, child_node, *iter);
      }*/
    }
  }

  void FullMap::addAssignment(AssignmentPtr assignment)
  {
    // 0. Create a dummy backbone.
    // 1. If the backbone corresponding to the assignment is included in the containers.
    BackbonePtr single_bone = boost::make_shared<Backbone>(assignment->getBackboneModificationSites(_mod_type), _mod_type);
    string clv_type = assignment->getFragmentType();
    if(clv_type == "A" || clv_type == "B" || clv_type == "C"){
      // Try to insert the assignment into _nre_set.
      if(this->findNRECleavage(assignment)) {

      } else {
        this->addBackbone(single_bone);
      }
    } else if(clv_type == "X" || clv_type == "Y" || clv_type == "Z") {
      // Try to insert the assignment into _re_set.
      if(this->findRECleavage(assignment)) {

      } else {
        this->addBackbone(single_bone);
      }
    } else {
      // Do nothing.
    }


  }

  void FullMap::removeAssignment(AssignmentPtr assignment, BackbonePtr bone)
  {
    bone->removeAssignment(assignment);
    if(bone->isEmptyNode()) {
      this->removeBackbone(bone);
    }
  }

  void FullMap::removeBackbone(BackbonePtr bone)
  {
    set<BackbonePtr> parents = bone->getParents();
    set<BackbonePtr> children = bone->getChildren();

    for(auto it1 = parents.begin(); it1 != parents.end(); it1++)
    {
      for(auto it2 = children.begin(); it2 != children.end(); it2++)
      {
        (*it1)->replaceChild(bone, *it2);
        (*it2)->replaceParent(bone, *it1);
      }
    }

    if(bone->isNRECleavage())
      _nre_set.erase(bone);
    else if(bone->isRECleavage())
      _re_set.erase(bone);
    else {
      _nre_set.erase(bone);
      _re_set.erase(bone);
    }
  }

  void FullMap::addBackbone(BackbonePtr bone)
{
    this->exploreEntryPoint(bone, _empty_node);
    if(bone->isRECleavage())
      _re_set.insert(bone);
    else if(bone->isNRECleavage())
      _nre_set.insert(bone);
    else {
      _re_set.insert(bone);
      _nre_set.insert(bone);
    }

  }

  set<BackbonePtr> FullMap::getPotentialParents(AssignmentPtr assignment)
  {

  }

  set<BackbonePtr> FullMap::getPotentialChildren(AssignmentPtr assignment)
  {

  }

  ostream& operator<<(ostream& os, const FullMap& graph)
  {
      vector<const BackbonePtr> bone_vec;
      bone_vec.push_back(graph._empty_node);
      //std::copy(graph._bone_set.begin(), graph._bone_set.end(), std::back_inserter(bone_vec));
      for(auto iter = graph._bone_set.begin(); iter != graph._bone_set.end();iter++)
      {
        if((*iter)->mod_sites == graph._full_node->mod_sites)
          continue;
        else
          bone_vec.push_back(*iter);
      }
      bone_vec.push_back(graph._full_node);

      // Iterate over all backbones
      
      for(auto iter = bone_vec.begin(); iter != bone_vec.end(); iter++)
      {
        if((*iter)->mod_sites == graph._full_node->mod_sites)
          continue;

        os << "Backbone:" << **iter << "\n";
        set<BackbonePtr> parents = (*iter)->getParents();
        for(auto p_it = parents.begin(); p_it != parents.end(); p_it++)
          os << "Parent:" << **p_it << "\n";

        set<BackbonePtr> children = (*iter)->getChildren();
        for(auto c_it = children.begin(); c_it != children.end(); c_it++)
          os << "Children:" << **c_it << "\n";

      }

      return os;
  }


}