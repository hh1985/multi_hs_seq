/********************************************************************
	created:	2014/04/28
	created:	28:4:2014   20:12
	filename: 	E:\Project\multi_hs_seq\GAG\src\GAGPL\GAGLIBRARY\AssignmentPool.h
	file path:	E:\Project\multi_hs_seq\GAG\src\GAGPL\GAGLIBRARY
	file base:	AssignmentPool
	file ext:	h
	author:		Han Hu
	
	purpose:	Container of assignments. It provides initial sorting and hashing
            functions for assignments.
*********************************************************************/

#ifndef GAG_ASSIGNMENTPOOL_H
#define GAG_ASSIGNMENTPOOL_H

//#include <GAGPL/GAGLIBRARY/Assignment.h>
//#include <GAGPL/GAGLIBRARY/Backbone.h>
#include <GAGPL/SEQUENCING/Division.h>
//#include <GAGPL/SEQUENCING/DivisionCluster.h>
#include <boost/multi_index_container.hpp>

namespace gag
{
  using namespace std;
  using namespace ::boost;
  using namespace ::boost::multi_index;

  //struct theo_mass {};
  struct ac_sites {};
  struct sulfate_sites {};
  struct general_type {};

  typedef multi_index_container<
    AssignmentPtr,
    indexed_by<
      ordered_non_unique<
      tag<theo_mass>, const_mem_fun<Assignment, double, &Assignment::getMass>
      >,
      ordered_non_unique<
      tag<ac_sites>, const_mem_fun<Assignment, ModificationSites, &Assignment::getAcetateSites>
      >,
      ordered_non_unique<
      tag<sulfate_sites>, const_mem_fun<Assignment, ModificationSites, &Assignment::getSulfateSites>
      >
    >
  > AssignmentContainer;

  typedef AssignmentContainer::index<theo_mass>::type AssignmentsByMass;
  typedef AssignmentContainer::index<ac_sites>::type AssignmentsByACSites;
  typedef AssignmentContainer::index<sulfate_sites>::type AssignmentsBySulfateSites;


  // Wrapper of assignment container.
  class AssignmentPool
  {
    // FullMap needs to access the assignment container.
    //friend class FullMap;

  public:
    //AssignmentPool(const ModificationSites& mod_sites, int mod_num)
    //  : full_sites(mod_sites), full_num(mod_num)
    //{}
    AssignmentPool() {}
    
    void addAssignment(AssignmentPtr);

    // 1. Getting isomeric assignment structure.
    set<AssignmentPtr> getAssignmentsByMass(double mass);

    // 2. Getting assignments with alternative structure configuration.
    set<AssignmentPtr> getAssignmentsByModificationSites(const string& mod_symbol, const ModificationSites&);

    // 3. Select qualified assignments for constructing assignment graph.
    // For each mass, only the terminal assignments can be selected.
    // For assignments with the the same modification sites, only the one with the largest modification number can be selected.
    //set<BackbonePtr> selectQualifiedAssignments(const string& mod_symbol);

    // Only terminal assignments can be selected!
    // The grouping procedure is left to the DivMap step.
    set<DivisionPtr> selectAssignments(const string& mod_symbol, const ModificationSites& mod_sites, int mod_num);

    // Get the uniqueness value by assignment.
    double getUniquenessValue(AssignmentPtr assignment);

  private:
    AssignmentContainer _pool;

  };
}
#endif /* GAG_ASSIGNMENTPOOL_H */