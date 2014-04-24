#include <GAGPL/GAGLIBRARY/CleavageTree.h>
#include <boost/make_shared.hpp>

namespace gag
{
	

	CleavageTree::CleavageTree( GlycanSequencePtr gs_ptr, std::string mod_symbol)
	{
		// Insert the ending node.
		ModificationSites complete_sites = gs_ptr->getModificationSitesBySymbol(mod_symbol, 1);
		ModificationSites empty_sites;

		CleavageItem begin_clv(empty_sites, 0);
		CleavageItem end_clv(complete_sites, gs_ptr->getModificationConstraint(mod_symbol));
		begin_clv.uniqueness_confidence = 1.0; end_clv.uniqueness_confidence = 1.0;

		begin_node = boost::make_shared<CleavageNode>(begin_clv);
		CleavageNodePtr end_node = boost::make_shared<CleavageNode>( end_clv);

		clv_tree.insert(std::make_pair(end_node->getModificationSites(), end_node));
		begin_node->addSuperSitesNode(end_node);
		end_node->addSubSitesNode(begin_node);

	}

	std::set<CleavageNodePtr> CleavageTree::getSubsets( CleavageItem& clv )
	{
		// Create the corresponding cleavage node with no parent and children information.
		if(clv.mod_sites.size() == 0) return std::set<CleavageNodePtr>();
				
		CleavageNode clv_node(clv);

		std::map<ModificationSites, CleavageNodePtr>::iterator node_iter = clv_tree.find(clv_node.getModificationSites());
    //std::set<CleavageNodePtr> clv_set;
		if(node_iter != clv_tree.end()) {
			return node_iter->second->last;
    } else {
      CleavageNodePtr temp_ptr = boost::make_shared<CleavageNode>(clv);
			CleavageEdgeSet clv_node_map = this->searchTree(temp_ptr);
			
      std::set<CleavageNodePtr> clv_set;
			for(CleavageEdgeSet::iterator iter = clv_node_map.begin();
				iter != clv_node_map.end(); iter++)
			{
				clv_set.insert(iter->first);
			}
			return clv_set;
		}
	}

	std::set<CleavageNodePtr> CleavageTree::getSubsets( CleavageNodePtr node )
	{
		if(node == nullptr) return std::set<CleavageNodePtr>();
		if(node->getModificationSites().size() == 0) return begin_node->last;

		std::map<ModificationSites, CleavageNodePtr>::iterator iter = clv_tree.find(node->getModificationSites());
		if(iter != clv_tree.end())
			return iter->second->last;
		else {
			CleavageEdgeSet clv_node_map = this->searchTree(node);
			std::set<CleavageNodePtr> clv_set;
			for(CleavageEdgeSet::iterator iter = clv_node_map.begin();
				iter != clv_node_map.end(); iter++)
			{
				clv_set.insert(iter->first);
			}
			return clv_set;
		}

	}

	std::set<CleavageNodePtr> CleavageTree::getSupersets( CleavageItem& clv )
	{
		// Create the corresponding cleavage node with no parent and children information
		if(clv.mod_sites.size() == 0) return begin_node->next;

		CleavageNode clv_node(clv);

		std::map<ModificationSites, CleavageNodePtr>::iterator iter = clv_tree.find(clv_node.getModificationSites());
		if(iter != clv_tree.end())
			return iter->second->next;
		else {

			CleavageEdgeSet clv_node_map = this->searchTree(boost::make_shared<CleavageNode>(clv));
			std::set<CleavageNodePtr> clv_set;
			for(CleavageEdgeSet::iterator iter = clv_node_map.begin();
				iter != clv_node_map.end(); iter++)
			{
				clv_set.insert(iter->second);
			}
			return clv_set;
		}
	}

	std::set<CleavageNodePtr> CleavageTree::getSupersets( CleavageNodePtr node )
	{
		if(node == nullptr) return std::set<CleavageNodePtr>();
		if(node->getModificationSites().size() == 0) return begin_node->next;

		std::map<ModificationSites, CleavageNodePtr>::iterator iter = clv_tree.find(node->getModificationSites());
		if(iter != clv_tree.end())
			return iter->second->next;
		else {
			CleavageEdgeSet clv_node_map = this->searchTree(node);
			std::set<CleavageNodePtr> clv_set;
			for(CleavageEdgeSet::iterator iter = clv_node_map.begin();
				iter != clv_node_map.end(); iter++)
			{
				clv_set.insert(iter->second);
			}
			return clv_set;
		}
			
	}

	bool CleavageTree::insertNode(CleavageItem& clv )
	{
    // If the node has been recorded, simply ignore the new node.
    CleavageNodePtr locate_node = this->locateModificationSites(clv.mod_sites);
    if(locate_node != nullptr) {
      locate_node->addIsoNode(clv);
      return true;
    }

		// Create new cleavage node.
		CleavageNodePtr clv_node = boost::make_shared<CleavageNode>(clv);

		// Locate the position on the tree.
		// this->searchChildren(begin_node, begin_node, clv_node);
		CleavageEdgeSet clv_node_map = this->searchTree(clv_node);

		std::set<CleavageNodePtr> clv_node_set;
		
    // Add the set into clv_tree, and set up the relationship.
		for(CleavageEdgeSet::iterator iter = clv_node_map.begin(); iter != clv_node_map.end(); iter++)
		{
			clv_node_set.insert(iter->first);
		}

		for(std::set<CleavageNodePtr>::iterator iter = clv_node_set.begin(); iter != clv_node_set.end(); iter++)
		{
			std::pair<CleavageEdgeSet::iterator, CleavageEdgeSet::iterator> p = clv_node_map.equal_range(*iter);
			while(p.first != p.second) {
				clv_node->addSubSitesNode(p.first->first);
				clv_node->addSuperSitesNode(p.first->second);
				p.first->first->replaceSuperNodeLink(p.first->second, clv_node);
				p.first->second->replaceSubNodeLink(p.first->first, clv_node);
				++p.first;
			}
		}

		clv_tree.insert(std::make_pair(clv_node->getModificationSites(), clv_node));

		return true;
	}

	void CleavageTree::searchChildren( CleavageNodePtr parent, CleavageNodePtr current, CleavageNodePtr new_node, CleavageEdgeSet& temp_store )
	{
		std::set<CleavageNodePtr> children = current->getSuperModificationSites();
		
		// Iterate over all parents, try to find the locate the new_node. 
		bool pass = false;
		for(std::set<CleavageNodePtr>::iterator iter = children.begin(); iter != children.end(); iter++) {
			if(new_node->hasSuperSetNode(*iter)) {
				
				// Bingo !!! Update status of the new_node.
				// 1. The child node is now new_node's child. The node current is node new_node's parent.
				//new_node->addSuperSitesNode(*iter);
				//new_node->addSubSitesNode(parent);
				
				// 2. new_node is the child node's parent.
				//iter->replaceSuperNodeLink(current, new_node);
				
				// 3. new_node is current's child.
				//current->replaceSubNodeLink(*iter, new_node);
				
				// 4. Store the node.
				//clv_tree.insert(new_node);
				temp_store.insert(std::make_pair(parent, *iter));

				if(!pass) pass = true;
			} else if(new_node->hasSubSetNode(*iter)) {
				searchChildren(*iter, *iter, new_node, temp_store);
				return;
			} else {
				continue;
			}
		}

		// If none of them works, append the new node directly to current.
		if(pass) return;

		for(std::set<CleavageNodePtr>::iterator iter = children.begin(); iter != children.end(); iter++) {
			searchChildren(parent, *iter, new_node, temp_store);
		}
		
	}

	CleavageEdgeSet CleavageTree::searchTree( CleavageNodePtr new_node )
	{

		CleavageEdgeSet temp_store;

		this->searchChildren(begin_node, begin_node, new_node, temp_store);
		
		return temp_store;
	}

  CleavageNodePtr CleavageTree::locateModificationSites( const ModificationSites& mod_sites )
  {
		// Create the object directly from modification sites.
		// CleavageItem clv_item(mod_sites);
		// CleavageNode matched_node(clv_item);

		std::map<ModificationSites, CleavageNodePtr>::iterator iter = clv_tree.find(mod_sites);
    
    if(iter == clv_tree.end())
      return boost::shared_ptr<CleavageNode>();
    else
      return iter->second;
	}

	void CleavageNode::replaceSubNodeLink( CleavageNodePtr old_node, CleavageNodePtr new_node )
	{
		NodeSet::iterator iter = last.find(old_node);
		
		if(iter != last.end()) last.erase(iter);

		this->addSubSitesNode(new_node);
	}

	void CleavageNode::replaceSuperNodeLink( CleavageNodePtr old_node, CleavageNodePtr new_node )
	{
		NodeSet::iterator iter = next.find(old_node);
		if(iter != next.end()) { // If found, erase the old link.
			//*iter = new_node;
			next.erase(iter);
		}
		this->addSuperSitesNode(new_node);
	}

}

