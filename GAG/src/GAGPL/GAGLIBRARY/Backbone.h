/********************************************************************
	created:	2014/04/30
	created:	30:4:2014   13:46
	filename: 	Backbone.h
	file path:	GAGPL\GAGLIBRARY
	file base:	Backbone
	file ext:	h
	author:		Han Hu
	
	purpose:	Backbone is a container for assignments with the same modification sites. Note that the assignments may actually come from different fragment types and may be virtual.
  1. Backbone connects to other backbone based on their modification sites.
  2. The compatibility between assignments is determined by both the modification sites information and the mod number information, which is monitored by FullMap object.
*********************************************************************/

#ifndef GAG_BACKBONE_H
#define GAG_BACKBONE_H

#include <GAGPL/GAGLIBRARY/Assignment.h>
#include <boost/shared_ptr.hpp>
//#include <boost/bimap/bimap.hpp>
//#include <boost/bimap/multiset_of.hpp>
//#include <boost/bimap/unordered_multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/enable_shared_from_this.hpp>

            
namespace gag
{
  class Backbone;
  typedef boost::shared_ptr<Backbone> BackbonePtr;

  class Backbone : public ModificationSites
  {
  public:
    // Assignments which share the same backbone structure.
    // The mod number is no longer important since it has been separated.
  
  public:
    Backbone(const ModificationSites& mod_sites, int top_n)
      : ModificationSites(mod_sites), top_num(top_n) {}

    //Backbone(const ModificationSites& mod_sites, int mod_num)
    //  : mod_sites(mod_sites) {
    //    members.insert(make_pair(mod_num, set<AssignmentPtr>()));
    //}

    // Constructor.  Either create a native backbone or a complementary backbone.
    //Backbone(AssignmentPtr assignment, const string& mod_symbol, bool comp=false)
    //{
    //  this->addAssignment(assignment, mod_symbol);
    //}

    BackbonePtr self()
    {
     return shared_from_this();
    }
    boost::shared_ptr<Backbone const> self() const
    {
     return shared_from_this();
    }

    /* Operation of internal members. */
    void addAssignment(AssignmentPtr assignment);
    //void addAssignment(AssignmentPtr assignment, const string& mod_symbol);
    void removeAssignment(AssignmentPtr assignment);
    //void removeAssignment(AssignmentPtr assignment, const string& mod_symbol);

    /* Operation of other backbones */
    //bool addFamily(BackbonePtr child, BackbonePtr parent);
    // In order to set up the relationship from parent node's side, work on the replace*** method from next_node object.
    void replaceParent(BackbonePtr last_node, BackbonePtr next_node);
    void replaceChild(BackbonePtr last_node, BackbonePtr next_node);

    set<BackbonePtr>& getParents();
    set<BackbonePtr>& getChildren();

    // The operation of adding parents and children is costly.  All the children will be used for making up the pairs. Notice that if the child set is empty. The null pointer will be added. This is useful in the case of appending.
    void addParent(BackbonePtr node);
    void addChild(BackbonePtr node);
    void removeParent(BackbonePtr node);
    void removeChild(BackbonePtr node);

    bool isSmaller(BackbonePtr cur)
    {
      return containSubset(cur, *this) && (this->size() < cur->size());
    }
    bool isLarger(BackbonePtr cur)
    {
      return containSubset(*this, cur) && (this->size() > cur->size());
    }

    /* Decide status. */
    inline bool isEmpty() const
    {
      return members.size() == 0;
    }

    
    friend ostream& operator<<(ostream& os, const Backbone& bone);
  
  private:

    int top_num;

    // neighbors are used for recording the context of internal cleavages (nominal)
    //BackboneNeighbor _neighbors;

    set<BackbonePtr> _parents;
    set<BackbonePtr> _children;

    set<AssignmentPtr> members;
    
  };
}

#endif /* GAG_BACKBONE_H */