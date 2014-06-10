/*!
 * \file DivisionCluster.h
 *
 * \author Han
 * \date June 2014
 *
 */

#include <GAGPL/SEQUENCING/Division.h>

namespace gag
{

  class DivisionCluster 
  {

  public:
    //DivisionCluster() {}
    DivisionCluster(const string& symbol, const ModificationSites& mod_sites, int mod_num)
      : _symbol(symbol), full_sites(mod_sites), full_num(mod_num), _size(0)  
    {}

    void addAssignment(AssignmentPtr assign);
    void addDivision(DivisionPtr div);

    // Only select division with the highest number for each modification sites tag.
    set<DivisionPtr> getSelectedDivision(bool highest);
    set<const DivisionPtr> getSelectedDivision(bool highest) const;

    set<const DivisionPtr> getAllDivision() const;

    size_t size() const;

    bool isQualifiedDivision(DivisionPtr div_ptr) const;
  
  private:
    DivisionPtr locateDivision(const ModificationSites& sites, int num);

  private:

    string _symbol;
    ModificationSites full_sites;
    int full_num;
    
    map<ModificationSites, map<int, DivisionPtr>> _backbone;
    size_t _size;
  };
}