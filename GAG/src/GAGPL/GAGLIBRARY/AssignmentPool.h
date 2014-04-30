/********************************************************************
	created:	2014/04/28
	created:	28:4:2014   20:12
	filename: 	E:\Project\multi_hs_seq\GAG\src\GAGPL\GAGLIBRARY\AssignmentPool.h
	file path:	E:\Project\multi_hs_seq\GAG\src\GAGPL\GAGLIBRARY
	file base:	AssignmentPool
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef GAG_ASSIGNMENTPOOL_H
#define GAG_ASSIGNMENTPOOL_H

#include <GAGPL/GAGLIBRARY/Assignment.h>
//#include <GAGPL/GAGLIBRARY/FullMap.h>
#include <boost/multi_index_container.hpp>
//#include <boost/shared_ptr.hpp>


namespace gag
{
    using namespace ::boost;
    using namespace ::boost::multi_index;

    struct theo_mass {};
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
            hashed_non_unique<
            tag<general_type>, const_mem_fun<Assignment, string, &Assignment::getCleavageType>
            >
        >
    > AssignmentContainer;

    typedef AssignmentContainer::index<theo_mass>::type AssignmentsByMass;
    typedef AssignmentContainer::index<ac_sites>::type AssignmentByACSites;
    typedef AssignmentContainer::index<sulfate_sites>::type AssignmentBySulfateSites;
    typedef AssignmentContainer::index<general_type>::type AssignmentsByCleavageType;


    // Wrapper of assignment container.
    class AssignmentPool
    {
    public:
        // 1. Getting isomeric assignment structure.
        set<AssignmentPtr> getAssignmentsByMass(double mass);

        // 2. Getting assignments with alternative structure configuration.
        set<AssignmentPtr> getAssignmentsByModificationSites(const std::string&, const ModificationSites&);

        // 3. The definition of uniqueness is based on ambiguous modificaiton sites.  If an assignment does not cause any ambiguity, it is unique.
        map<int, set<AssignmentPtr>> sortAssignmentsByUniqueness(ModificationSites&);

        void addAssignment(AssignmentPtr);
        
        // FullMap needs to access assignment container.
        friend class FullMap;

        // Calculate the ambiguity of the assignment set. Note that unqualified assignments (e.g. cross-ring/cross-ring assignment) are not maintained.
        pair<int, set<AssignmentPtr>> calculateUniquenessValue(set<AssignmentPtr>);


    private:
        AssignmentContainer _pool;
    };
}
#endif /* GAG_ASSIGNMENTPOOL_H */