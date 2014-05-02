#include <GAGPL/GAGLIBRARY/FullMap.h>

namespace gag
{

  void FullMap::connectBackboneSet()
  {
    for(auto bone_iter = _bone_set.begin(); bone_iter != _bone_set.end(); bone_iter++)
    {
      /*set<BackbonePtr> parents = _empty_node->getParents();
      for(auto p_iter = parents.begin(); p_iter != parents.end(); p_iter++)
      this->exploreDeepNodes(*bone_iter, _empty_node, *p_iter);*/
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

    //_bone_set.insert(_empty_node);
    //_bone_set.insert(_full_node);
  }

  void FullMap::exploreDeepNodes(BackbonePtr cur, BackbonePtr child_node, BackbonePtr parent_node)
  {


#ifdef _DEBUG
      cout << "\n******Before*****\n"; 
      //cout << "Current:" << cur->mod_sites << "\n";
      //for(auto iter = cur->getParents().begin(); iter != cur->getParents().end(); iter++)
      //  cout << "--Checked parent:" << (*iter)->mod_sites << "\n";
      //for(auto iter = cur->getChildren().begin(); iter != cur->getChildren().end(); iter++)
      //  cout << "--Checked child:" << (*iter)->mod_sites << "\n";
      //for(auto iter = cur->getSiblings().begin(); iter != cur->getSiblings().end(); iter++)
      //  cout << "--Checked sibling:" << (*iter)->mod_sites << "\n";

      //cout << "Checked:" << child_node->mod_sites << "\n";
      //for(auto iter = child_node->getParents().begin(); iter != child_node->getParents().end(); iter++)
      //  cout << "--Checked parent:" << (*iter)->mod_sites << "\n";
      //for(auto iter = child_node->getChildren().begin(); iter != child_node->getChildren().end(); iter++)
      //  cout << "--Checked child:" << (*iter)->mod_sites << "\n";
      //for(auto iter = child_node->getSiblings().begin(); iter != child_node->getSiblings().end(); iter++)
      //  cout << "--Checked sibling:" << (*iter)->mod_sites << "\n";
      cout << "Current:" << *cur << "\n";
      cout << "Checked child:" << *child_node << "\n";
      cout << "Checked parent:" << *parent_node << "\n";

#endif // _DEBUG
      // Get parent nodes.
      set<BackbonePtr> parents = child_node->getParents();

      // Iterate over all parents. Decide if there is any chance to locate the position of the current node.
      
      if(parent_node->isLarger(cur)){ // To the boundary.
        // Removing the original connection between child:tree_node, parent:*iter
        child_node->replaceParent(parent_node, cur);
        parent_node->replaceChild(child_node, cur);

        // Insert the node.
        cur->addChild(child_node);
        cur->addParent(parent_node);
      } else {
        set<BackbonePtr> parents_parent = parent_node->getParents();
        if(parent_node->isSmaller(cur)) {
          for(auto p_iter = parents_parent.begin(); p_iter != parents_parent.end(); p_iter++)
          {
            // Adjust the child node and parent node.
            this->exploreDeepNodes(cur, parent_node, *p_iter);
          }
        } else {
          for(auto p_iter = parents_parent.begin(); p_iter != parents_parent.end(); p_iter++)
          {
            // Adjust the child node and parent node.
            this->exploreDeepNodes(cur, child_node, *p_iter);
          }

          // Check if parent_node is a sibling node.
          std::vector<BackbonePtr> common_parents;
          std::vector<BackbonePtr> common_children;

          set_intersection(parent_node->getParents().begin(), parent_node->getParents().end(), cur->getParents().begin(), cur->getParents().end(), std::back_inserter(common_parents));
          set_intersection(parent_node->getChildren().begin(), parent_node->getChildren().end(), cur->getChildren().begin(), cur->getChildren().end(), std::back_inserter(common_children));
          if(!common_parents.empty() && !common_children.empty())
          {
#ifdef _DEBUG
            cout << "Found a new sibling!\n";
#endif // _DEBUG

            cur->addSibling(parent_node);
            parent_node->addSibling(cur);

          }
        }

      } 


      //// Get parent nodes.
      //set<BackbonePtr> parents = child_node->getParents();

      //// Iterate over all parents. Decide if there is any chance to locate the position of the current node.
      //for(auto iter = parents.begin(); iter != parents.end(); iter++)
      //{
      //  if((*iter)->isLarger(cur)){ // To the boundary.
      //    // Removing the original connection between child:tree_node, parent:*iter
      //    child_node->replaceParent(*iter, cur);
      //    (*iter)->replaceChild(child_node, cur);

      //    // Insert the node.
      //    cur->addChild(child_node);
      //    cur->addParent(*iter);
      //  } else if((*iter)->isSmaller(cur)){
      //    set<BackbonePtr> parents_parents = (*iter)->getParents();
      //    for(auto p_iter = parents_parents.begin(); p_iter != parents_parents.end(); p_iter++)
      //    {
      //      // Adjust the child node and parent node.
      //      this->exploreDeepNodes(cur, *iter, *p_iter);
      //    }

      //  } else {
      //    set<BackbonePtr> parents_parents = (*iter)->getParents();
      //    for(auto p_iter = parents_parents.begin(); p_iter != parents_parents.end(); p_iter++)
      //    {
      //      // Adjust only the parent node.
      //      this->exploreDeepNodes(cur, child_node, *p_iter);
      //    }
      //  }
      //}

      // Check all the parents of the child node, decide if there is any sibling.

      

        //if(!cur->isSibling(*iter)){ // If recorded as sibling.

        //  // Add children.
        //  set<BackbonePtr> sib_children = (*iter)->getChildren();
        //  for(auto child_iter = sib_children.begin(); child_iter != sib_children.end(); child_iter++)
        //  {
        //    cur->addChild(*child_iter);
        //    (*child_iter)->addParent(cur);
        //  }

        //  // Add parents, all the parents of the sibling is also the parents of the 
        //  set<BackbonePtr> sib_parents = (*iter)->getParents();
        //  for(auto parent_iter = sib_parents.begin(); parent_iter != sib_parents.end(); parent_iter++)
        //  {
        //    cur->addParent(*parent_iter);
        //    (*parent_iter)->addChild(cur);
        //  }

        //  // Add sibling
        //  cur->addSibling(*iter);
        //  (*iter)->addSibling(cur);

        //  set<BackbonePtr> sib = (*iter)->getSiblings();
        //  for(auto sib_iter = sib.begin(); sib_iter != sib.end(); sib_iter++)
        //  {
        //    cur->addSibling(*sib_iter);
        //    (*sib_iter)->addSibling(cur);
        //  }


        //} else {
        //  // Should be checked already.
        //}
      //}
#ifdef _DEBUG
      cout << "\n******After*****\n"; 
      cout << "Current:" << *cur << "\n";
      cout << "Checked child:" << *child_node << "\n";
      cout << "Checked parent:" << *parent_node << "\n";
#endif // _DEBUG
  }

  void FullMap::exploreCompatibility( BackbonePtr bone)
  {
    set<BackbonePtr>& parents = bone->getParents();

    set<BackbonePtr>& children = bone->getChildren();

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

  bool FullMap::exploreEntryPoint( BackbonePtr cur, BackbonePtr child_node, set<BackbonePtr> parent_nodes)
  {
    bool sib = true;
    for(auto iter = parent_nodes.begin(); iter != parent_nodes.end(); iter++)
    {
      if((*iter)->isLarger(cur)){
        sib = false;
      } else if((*iter)->isSmaller(cur)) {
        sib = false;
      } 
    }

    if(sib) {
      // Create the connection between cur and child_node.
      cur->addChild(child_node);
      child_node->addParent(cur);

      set<BackbonePtr> new_parents;
      for(auto iter = parent_nodes.begin(); iter != parent_nodes.end(); iter++)
      {

      }
    }

    return sib;
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

          set<BackbonePtr>& siblings = (*iter)->getSiblings();
          for(auto s_it = children.begin(); s_it != children.end(); s_it++)
              os << "Siblings:" << **s_it << "\n";

      }
      return os;
  }


}