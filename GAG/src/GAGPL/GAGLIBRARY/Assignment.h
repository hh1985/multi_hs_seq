/********************************************************************
	created:	2014/04/25
	created:	25:4:2014   16:04
	filename: 	Assignment.h
	file path:	GAG\src\GAGPL\GAGLIBRARY
	file base:	Assignment
	file ext:	h
	author:		Han Hu
	
	purpose:	Assignment describes the structural information regarding sulfate distribution: candidate modification sites and modification number.
*********************************************************************/

#ifndef GAG_ASSIGNMENT_H
#define GAG_ASSIGNMENT_H

//#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/CHEMISTRY/Modification.h>
#include <GAGPL/SPECTRUM/MonoPeak.h>
#include <boost/shared_ptr.hpp>

namespace gag
{
  // Only concern about the neutral loss, 
  // It has nothing to do with the sulfate and acetate number.

  using namespace std;
  using namespace ::boost;
  
  //using namespace ::boost::multi_index;
  

  class Assignment: public Unit {

  public:
    Assignment(MonoPeakPtr pk, NodeItem& node)
      : pk(pk), fg(node.getFragment()), Unit(node.getComposition()) 
    {
        const Satellite& sate = node.getCompositionShiftList();

        auto iter = sate.begin();
        while(iter != sate.end())
        {
            if(iter->first == "C2H2O")
                mod_count.insert(make_pair("Ac", abs(iter->second)));
            else if(iter->first == "O3S")
                mod_count.insert(make_pair("SO3", abs(iter->second)));
            else
                neu_loss.insert(*iter);

            iter++;
        }
    }
    
    // Used for creating dummy assisgnment.
    Assignment() {}

    // copy constructor.
    Assignment(const Assignment&);
    
    // operator=
    Assignment& operator=(const Assignment&);

    int getNeutralLossNumber(const std::string& loss) const
    {
      auto iter = neu_loss.find(loss);
      return iter == neu_loss.end() ? 0 : iter->second;
    }

    int getModificationNumber(const std::string& mod) const
    {
      auto iter = mod_count.find(mod);
      return iter == mod_count.end() ? 0 : iter->second;
    }

    inline void setModificationNumber(const std::string& mod, int num)
    {
      mod_count[mod] = num;
    }

    inline ModificationSites getBackboneModificationSites(const std::string& mod) const
    {
        return fg->getModificationSitesBySymbol(mod, 1);
    }

    inline ModificationSites getAcetateSites() const
    {
        return fg->getModificationSitesBySymbol("Ac", 1);
    }
    inline ModificationSites getSulfateSites() const
    {
        return fg->getModificationSitesBySymbol("SO3", 1);
    }

    inline double getMass() const
    {
        return fg->getMass();
    }

    /*inline void addParent(AssignmentPtr parent)
    {
        _parents.insert(parent);
    }

    inline void addChild(AssignmentPtr child)
    {
        _children.insert(child);
    }

    inline void addSibling(AssignmentPtr sibling)
    {
        _siblings.insert(sibling);
    }*/

    string getCleavageType()
    {
        return fg->getGeneralType();
    }

    /*inline size_t getChildNumber()
    {
        return _parents.size();
    }

    inline size_t getParentNumber()
    {
        return _children.size();
    }

    inline size_t getSiblingNumber()
    {
        return _siblings.size();
    }*/

    pair<ModificationSites, int> describeStructure(const string&);

    friend ostream& operator<<(ostream& os, const Assignment& assignment);

  private:
    FragmentPtr fg;
    map<string, int> neu_loss;

    // For GAG, the modification contains only sulfate and acetate group.
    map<string, int> mod_count;

    MonoPeakPtr pk;

  };

  typedef boost::shared_ptr<Assignment> AssignmentPtr;
  

}


#endif /* GAG_ASSIGNMENT_H */