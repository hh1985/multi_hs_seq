#include <GAGPL/GAGLIBRARY/FullMap.h>

namespace gag
{

  void FullMap::connectBackboneSet()
  {
    
    for(auto bone_iter = _bone_set.begin(); bone_iter != _bone_set.end(); bone_iter++)
    {
      // Do not consider other full nodes.
      if((*bone_iter)->mod_sites == _full_node->mod_sites)
        continue;

      this->exploreEntryPoint(*bone_iter, _empty_node);
    }
  }

  void FullMap::initialize()
  {
    _empty_node = boost::make_shared<Backbone>(ModificationSites(), 0);
    
    _full_node = boost::make_shared<Backbone>(seq->getModificationSitesBySymbol(_mod_type, 1), seq->getModificationConstraint(_mod_type));
    //_full_node.addAssignment(end_assign);

    _empty_node->addParent(_full_node);
    _full_node->addChild(_empty_node);

    //_empty_node->addFamily(BackbonePtr(), _full_node);
    //_full_node->addFamily(_empty_node, BackbonePtr());

#ifdef _DEBUG
    cout << "Empty:\n";
    cout << *_empty_node << "\n";
    cout << "Full:\n";
    cout << *_full_node << "\n";
#endif // _DEBUG
    //_bone_set.insert(_empty_node);
    //_bone_set.insert(_full_node);
  }

  void FullMap::exploreDeepNodes(BackbonePtr cur, BackbonePtr child_node, BackbonePtr parent_node)
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
      
      if(parent_node->isLarger(cur)){ // To the boundary.
       /* child_node->replaceParent(parent_node, cur);
        parent_node->replaceChild(child_node, cur);*/

        // Append the node.
        //cur->addChild(child_node);
        //cur->addParent(parent_node);
        bool success_insertion = cur->addFamily(child_node, parent_node);
        
        // No need to check further information.
        if(success_insertion) {
          child_node->addParent(cur);
          parent_node->addChild(cur);
        }
      } else {
        set<BackbonePtr> parents_parent = parent_node->getParents();
        for(auto p_iter = parents_parent.begin(); p_iter != parents_parent.end(); p_iter++)
        {
          // Adjust the child node and parent node.
          this->exploreDeepNodes(cur, child_node, *p_iter);
        }

      } 
#ifdef _DEBUG
      cout << "\n******After*****\n"; 
      cout << "Current:" << *cur << "\n";
      cout << "Checked child:" << *child_node << "\n";
      cout << "Checked parent:" << *parent_node << "\n";
#endif // _DEBUG
  }

  void FullMap::exploreCompatibility( BackbonePtr bone)
  {
    set<BackbonePtr> parents = bone->getParents();

    set<BackbonePtr> children = bone->getChildren();

    for(auto iter = parents.begin(); iter != parents.end(); iter++)
    {
      this->checkCompatibility(bone, *iter);
      if((*iter)->isDummyNode())
        return;

      this->exploreCompatibility(*iter);
    }

    for(auto iter = children.begin(); iter != children.end(); iter++)
    {
      this->checkCompatibility(*iter, bone);
      this->exploreCompatibility(*iter);
    }

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
    for(auto iter = parents.begin(); iter != parents.end(); iter++)
    {
      if((*iter)->isLarger(cur)){ // Insertion.
        
        child_node->replaceParent(*iter, cur);
        (*iter)->replaceChild(child_node, cur);

        // Insert the node.
        //cur->addChild(child_node);
        //cur->addParent(*iter);
        cur->addFamily(child_node, *iter);

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
        this->exploreEntryPoint(cur, *iter);
        
        if(sib)
          sib = false;
      } 
    }

    // If all the paths are not acceptable. Creating a sibling.
    if(sib) { // Appending.
      // The child_node is fixed. Try to get the 
      for(auto iter = parents.begin(); iter != parents.end(); iter++)
      {
        this->exploreDeepNodes(cur, child_node, *iter);
      }
    }
  }

  ostream& operator<<(ostream& os, const FullMap& graph)
  {
      vector<const BackbonePtr> bone_vec;
      bone_vec.push_back(graph._empty_node);
      std::copy(graph._bone_set.begin(), graph._bone_set.end(), std::back_inserter(bone_vec));
      bone_vec.push_back(graph._full_node);

      // Iterate over all backbones
      
      for(auto iter = bone_vec.begin(); iter != bone_vec.end(); iter++)
      {
          os << "Backbone:" << **iter << "\n";
          set<BackbonePtr>& parents = (*iter)->getParents();
          for(auto p_it = parents.begin(); p_it != parents.end(); p_it++)
              os << "Parent:" << **p_it << "\n";

          set<BackbonePtr>& children = (*iter)->getChildren();
          for(auto c_it = children.begin(); c_it != children.end(); c_it++)
              os << "Children:" << **c_it << "\n";

      }
      return os;
  }


}