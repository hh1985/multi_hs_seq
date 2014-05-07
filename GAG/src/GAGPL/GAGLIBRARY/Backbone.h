/********************************************************************
	created:	2014/04/30
	created:	30:4:2014   13:46
	filename: 	Backbone.h
	file path:	GAGPL\GAGLIBRARY
	file base:	Backbone
	file ext:	h
	author:		Han Hu
	
	purpose:	Backbone groups the assignments by their modification 
            sites. 
*********************************************************************/

#ifndef GAG_BACKBONE_H
#define GAG_BACKBONE_H

#include <GAGPL/GAGLIBRARY/Assignment.h>
#include <boost/shared_ptr.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/unordered_multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

namespace gag
{
  class Backbone;
  typedef boost::shared_ptr<Backbone> BackbonePtr;
  
  struct child {};
  struct parent {};

  using namespace boost::bimaps;
  typedef bimap<
    unordered_multiset_of<tagged<BackbonePtr, child>>, 
    unordered_multiset_of<tagged<BackbonePtr, parent>>,
    set_of_relation<>
  > BackboneNeighbor;

  typedef BackboneNeighbor::value_type NeighborType;

  class Backbone
  {
  public:

    ModificationSites mod_sites;
    map<int, set<AssignmentPtr>> members;
  
  public:
    Backbone(const ModificationSites& mod_sites, int mod_num)
      : mod_sites(mod_sites) {
        members.insert(make_pair(mod_num, set<AssignmentPtr>()));
    }

    // Constructor.
    Backbone(AssignmentPtr assignment, const string& mod_symbol)
    {
      this->addAssignment(assignment, mod_symbol);
    }

    set<AssignmentPtr> getAssignmentsByModNumber(int mod_num);

    set<int> getModNumbers() const;

    int getLargestModNumber() const;

    void addAssignment(AssignmentPtr assignment, const string& mod_symbol);

    inline bool isDummyNode() const
    {
      return members.begin()->second.size() == 0;
    }

    inline bool operator<(const Backbone& bone)
    {
      return mod_sites < bone.mod_sites;
    }

    bool addFamily(BackbonePtr child, BackbonePtr parent);

    void replaceParent(BackbonePtr last_node, BackbonePtr next_node);

    void replaceChild(BackbonePtr last_node, BackbonePtr next_node);

    set<BackbonePtr> getParents();

    const set<BackbonePtr> getParents() const;

    set<BackbonePtr> getChildren();

    const set<BackbonePtr> getChildren() const;

    // The operation of adding parents and children is costly.  All the children will be used for making up the pairs. Notice that if the child set is empty. The null pointer will be added. This is useful in the case of appending.
    void addParent(BackbonePtr node);

    void addChild(BackbonePtr node);

    bool isSmaller(BackbonePtr cur)
    {
      return containSubset(cur->mod_sites, mod_sites) && (mod_sites.size() < cur->mod_sites.size());
    }

    bool isLarger(BackbonePtr cur)
    {
      return containSubset(mod_sites, cur->mod_sites) && (mod_sites.size() > cur->mod_sites.size());
    }

    friend ostream& operator<<(ostream& os, const Backbone& bone);
  
  private:
    // parents, children and siblings are only for terminal cleavages (nominal).
    //set<BackbonePtr> _parents;
    //set<BackbonePtr> _children;
    //set<BackbonePtr> _siblings;

    // neighbors are used for recording the context of internal cleavages (nominal)
    BackboneNeighbor _neighbors;
    
  };
}

#endif /* GAG_BACKBONE_H */