/********************************************************************
	created:	2014/04/25
	created:	25:4:2014   16:04
	filename: 	Assignment.h
	file path:	GAG\src\GAGPL\GAGLIBRARY
	file base:	Assignment
	file ext:	h
	author:		Han Hu
	
	purpose:	Assignment describes the structural information regarding sulfate distribution: candidate modification sites and
 modification number.
*********************************************************************/

#ifndef GAG_ASSIGNMENT_H
#define GAG_ASSIGNMENT_H

#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <GAGPL/CHEMISTRY/Modification.h>
#include <GAGPL/SPECTRUM/MonoPeak.h>

namespace gag
{
  // Only concern about the neutral loss, 
  // It has nothing to do with the sulfate and acetate number.

  //typedef std::map<std::string, int> NeutralLoss;
  using namespace std;

  class Assignment;
  typedef boost::shared_ptr<Assignment> AssignmentPtr;
  typedef set<AssignmentPtr> 

  class Assignment: public Unit {

  public:
    Assignment(FragmentPtr fg)
      : fg(fg), Unit(fg->getComposition()) 
    {}
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
      auto iter = mod_num.find(mod);
      return iter == mod_num.end() ? 0 : iter->second;
    }

    inline void setModificationNumber(const std::string& mod, int num)
    {
      mod[mod] = num;
    }

    // 
    ModificationSites getBackboneModificationSites() const;

  private:
    FragmentPtr fg;
    std::map<std::string, int> neu_loss;

    // For GAG, the modification 
    std::map<std::string, int> mod_num;

    
  };
}


#endif /* ASSIGNMENT_H */