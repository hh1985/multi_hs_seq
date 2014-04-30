/********************************************************************
	created:	2013/05/10
	created:	10:5:2013   13:29
	filename: 	CleavageTree.h
	file path:	GAGPL\GAGLIBRARY
	file base:	CleavageTree
	file ext:	h
	author:		Han Hu
	
	purpose:	Organize the relationship between different modification 
						sites.
*********************************************************************/
#ifndef GAG_CLEAVAGETREE_H
#define GAG_CLEAVAGETREE_H

#include <GAGPL/CHEMISTRY/Modification.h>
#include <GAGPL/GAGLIBRARY/GeneralCleavage.h>
#include <boost/shared_ptr.hpp>
//#include <boost/ptr_container/ptr_set.hpp>
#include <boost/ptr_container/ptr_map.hpp>

namespace gag
{
	// If m1 contain subset m2.
	//bool containSubset(const ModificationSites& m1, const ModificationSites& m2);

  class CleavageNode;
  //typedef boost::ptr_set<CleavageNode> NodeSet;
  
  typedef boost::shared_ptr<CleavageNode> CleavageNodePtr;

  typedef std::set<CleavageNodePtr> NodeSet;

	class CleavageNode
	{
	public:
		CleavageNode(CleavageItem& clv)
			: mod_sites(clv.mod_sites)
		{
      clv_vec.push_back(clv);
    }

    ~CleavageNode()
    {}

		inline bool operator<(const CleavageNode& clv_node) const
		{
			return mod_sites < clv_node.mod_sites;
		}

		inline void addSubSitesNode(CleavageNodePtr node)
		{
			last.insert(node);
		}
		
		inline void addSuperSitesNode(CleavageNodePtr node)
		{
			next.insert(node);
		}
    
    inline void addIsoNode(CleavageItem& clv)
    {
      clv_vec.push_back(clv);
    }

    NodeSet getSubModificationSites()
		{
			return last;
		}

		NodeSet getSuperModificationSites()
		{
			return next;
		}

		inline bool isEndingNode() const
		{
			return next.size() == 0;
		}

		inline bool hasSubSetNode(CleavageNodePtr clv_node)
		{
			return std::includes(mod_sites.begin(), mod_sites.end(), clv_node->mod_sites.begin(), clv_node->mod_sites.end());
		}
		inline bool hasSuperSetNode(CleavageNodePtr clv_node)
		{
			return std::includes(clv_node->mod_sites.begin(), clv_node->mod_sites.end(), mod_sites.begin(), mod_sites.end());
		}
    inline bool isIsoNode(CleavageNodePtr clv_node)
    {
      return mod_sites == clv_node->mod_sites;
    }

		void replaceSubNodeLink(CleavageNodePtr old_node, CleavageNodePtr new_node);
		void replaceSuperNodeLink(CleavageNodePtr old_node, CleavageNodePtr new_node);

    const std::vector<CleavageItem>& getCleavageItems() const
    {
      return clv_vec;
    }

    const ModificationSites& getModificationSites() const
    {
      return mod_sites;
    }

	public:
		// Super set. Has to be non-redundant.
		//std::set<CleavageNode*> last;
    NodeSet last;
		// Subset. Has to be non-redundant.
		//std::set<CleavageNode*> next;
    NodeSet next;

  private:
    // CleavageItem objects with the same mod type are grouped together.
    std::vector<CleavageItem> clv_vec;
    
    ModificationSites mod_sites;

	};

	typedef std::multimap<CleavageNodePtr, CleavageNodePtr> CleavageEdgeSet;
	typedef std::pair<CleavageNodePtr, CleavageNodePtr> CleavageBoundary;

	class CleavageTree
	{
	public:
		// Constructor. Insert the starting node and ending node.
		CleavageTree(GlycanSequencePtr gs_ptr, std::string mod_symbol);
		// Destructor.
		~CleavageTree() {
			//delete begin_node; begin_node = nullptr;
			//delete end_node; end_node = NULL;
		}

		// Find the right position for the new modification sites. If the node has been recorded in the tree, do nothing for that.
		// The return value indicates the results of insersion 
		// events.
		bool insertNode(CleavageItem& clv);

		// Search the children of the current node to see if there are matches. When searching for parents, no longer cut the position.
		void searchChildren(CleavageNodePtr parent, CleavageNodePtr current, CleavageNodePtr new_node, CleavageEdgeSet& temp_store);

		// Wrapper of searchChilden() function. This function returns a vector of pair of cleavage nodes. The first element of the pair is parent and the second is a child.
		CleavageEdgeSet searchTree(CleavageNodePtr new_node);
		
		std::set<CleavageNodePtr> getSubsets(CleavageItem& clv);
		std::set<CleavageNodePtr> getSubsets(CleavageNodePtr node);

		std::set<CleavageNodePtr> getSupersets(CleavageItem& clv);
		std::set<CleavageNodePtr> getSupersets(CleavageNodePtr node);

		inline CleavageNodePtr getBeginNode()
		{
			return begin_node;
		}

		// Locate the node from modification sites.
		// If no node has been found, return NULL instead.
    CleavageNodePtr locateModificationSites(const ModificationSites& mod_sites);

	private:
		// Reason for using ptr container:
		// 1. Good for destroy the tree.
		// 2. Good for find function.
		//NodeSet clv_tree;
    std::map<ModificationSites, CleavageNodePtr> clv_tree;

		CleavageNodePtr begin_node;
		
	};
}

#endif /* GAG_CLEAVAGETREE_H */