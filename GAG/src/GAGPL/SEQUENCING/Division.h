/*!
 * \file Division.h
 *
 * \author Han
 * \date June 2014
 *
 */

#ifndef GAG_DIVISION_H
#define GAG_DIVISION_H

#include <GAGPL/GAGLIBRARY/Assignment.h>
#include <boost/enable_shared_from_this.hpp>

namespace gag
{
/*!
 * \class Division
 *
 * \brief Division describes the pair of modification sites and modification number, which serves to divide the sequence into two regions.
 *
 * \author Han
 * \date June 2014
 */

  class Division;
  typedef boost::shared_ptr<Division> DivisionPtr;

  class Division: public boost::enable_shared_from_this<Division>
  {
  public:
    Division()
      : correction(0)
    {}
    Division(const ModificationSites& sites, int num)
      : mod_sites(sites), mod_num(num), correction(0)
    {}
    DivisionPtr self()
    {
      return shared_from_this();
    }

    boost::shared_ptr<Division const> self() const
    {
      return shared_from_this();
    }

    inline int getModificationNumber() const
    {
      return mod_num;
    }
    inline int getCorrectionNumber() const
    {
      return correction;
    }
    inline int getCorrectedModificationNumber() const
    {
      return mod_num + correction;
}

    inline int getModificationSitesNumber() const
    {
      return (int)mod_sites.size();
    }

    inline string getGeneralType() const
    {
      return _type;
    }

    // A - B
    int getDiffSitesNumber(const DivisionPtr div_node) const;
    int getInterSitesNumber(const DivisionPtr div_node) const;

    inline const ModificationSites& getModificationSites() const
    {
      return mod_sites;
    }

    inline void updateCorrectionNumber(int cor)
    {
      correction = cor;
    }

    inline bool isCorrected() const
    {
      return correction == 0;
    }

    void addAssignment(AssignmentPtr assign);

    inline const set<AssignmentPtr>& getSupportAssignments() const
    {
      return assign_support;
    }

    // Provide relationship examination. If A is B's parent, then B is automatically A's child.
    void addParent(DivisionPtr p_div);
    void addChild(DivisionPtr c_div);
    void addNode(DivisionPtr c_div, DivisionPtr p_div);

    void replaceParent(DivisionPtr old_div, DivisionPtr new_div);
    void replaceChild(DivisionPtr old_div, DivisionPtr new_div);

    set<DivisionPtr>& getParents();
    set<DivisionPtr>& getChildren();


    int getInDegree() const;
    int getOutDegree() const;

    bool isLargerThan(DivisionPtr div_ptr) const;
    bool isSmallerThan(DivisionPtr div_ptr)const ;
    bool isSibling(DivisionPtr div_ptr) const;
    bool isCompatible(DivisionPtr div_ptr) const;

    bool isClosestParent(DivisionPtr c_div) const ;
    bool isClosestChild(DivisionPtr) const;

    friend ostream& operator<<(ostream& os, const Division& div);
  
  private:
    // TBD: check the overall confidence score based on the assignments information.
    double getComprehensiveConfidence() const;
    double shiftCost() const;

    // Guarantee the new added node is the "best" node.
    // No qualification check.
    void updateParent(DivisionPtr new_ptr);
    void updateChild(DivisionPtr new_ptr);

    // Low level operation of child and parent set with no protection.
    void insertParent(DivisionPtr p_div);
    void insertChild(DivisionPtr c_div);
    void insertNode(DivisionPtr c_div, DivisionPtr p_div);

  private:
    /* 
     * Structure
     * For terminal cleavage, only the NRE part is stored. For internal cleavage, keep the original mod_sites and mod_num.
     */
    ModificationSites mod_sites;
    int mod_num;
    
    // "correction" specifies the nominal modification number. It is useful when considering about the case of complete sulfate loss, correction will record the manually corrected sulfate shift. 
    int correction;
    // Specify the cleavage type of the Division object. Possible values include: "C"(cross-ring), "G"(glycosidic-bond), "I"(internal). Note that the type is only determined by the modification sites, and has nothing to do with the original fragment type.
    string _type;

    /* 
     * Support 
     */
    set<AssignmentPtr> assign_support;

    /* 
     * Environment 
     */
    set<DivisionPtr> parent_set;
    set<DivisionPtr> child_set;

  };

}

#endif /* GAG_DIVISION_H */