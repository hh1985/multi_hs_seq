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

namespace gag
{
  class Backbone;
  typedef boost::shared_ptr<Backbone> BackbonePtr;

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

    inline size_t getParentsNumber() const
    {
      return _parents.size();
    }

    inline size_t getChildenNumber() const
    {
      return _children.size();
    }

    inline size_t getSiblingsNumber() const
    {
      return _siblings.size();
    }

    inline set<BackbonePtr>& getParents()
    {
      return _parents;
    }

    inline set<BackbonePtr>& getChildren()
    {
      return _children;
    }

    inline set<BackbonePtr>& getSiblings()
    {
      return _siblings;
    }

    inline void addParent(BackbonePtr parent)
    {
      _parents.insert(parent);
    }

    inline void addChild(BackbonePtr child)
    {
      _children.insert(child);
    }

    inline void addSibling(BackbonePtr sib)
    {
      _siblings.insert(sib);
    }

    void replaceParent(BackbonePtr old_bone, BackbonePtr new_bone)
    {
      _parents.erase(old_bone);
      _parents.insert(new_bone);
    }

    void replaceChild(BackbonePtr old_bone, BackbonePtr new_bone)
    {
      _children.erase(old_bone);
      _children.insert(new_bone);
    }

    bool isSibling(BackbonePtr cur)
    {
      return _siblings.find(cur) != _siblings.end();
    }

    bool isSmaller(BackbonePtr cur)
    {
      return containSubset(cur->mod_sites, mod_sites);
    }

    bool isLarger(BackbonePtr cur)
    {
      return containSubset(mod_sites, cur->mod_sites);
    }

    friend ostream& operator<<(ostream& os, const Backbone& bone);
  private:
    set<BackbonePtr> _parents;
    set<BackbonePtr> _children;
    set<BackbonePtr> _siblings;
  };
}

#endif /* GAG_BACKBONE_H */