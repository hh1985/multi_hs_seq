/********************************************************************
	created:	2014/04/29
	created:	29:4:2014   0:04
	filename: 	FullMap.h
	file path:	GAG\src\GAGPL\GAGLIBRARY
	file base:	FullMap
	file ext:	h
	author:		Han Hu
	
	purpose:	Connecting different assignments into a graph.  Note that 
*********************************************************************/

#ifndef GAG_FULLMAP_H
#define GAG_FULLMAP_H

//#include <GAGPL/GAGLIBRARY/Assignment.h>
#include <GAGPL/GAGLIBRARY/AssignmentPool.h>

namespace gag
{
    struct StatusItem
    {
        AssignmentPtr target;
        AssignmentPtr neighbor;
        int direction;
    };

    class FullMap
    {
    public:
        // Constructor.
        FullMap(const std::string& mod)
            : _mod_type(mod){}
        
        // Each time when an new assignment is inserted into the graph, check the graph topology and try to set up the connection.
        // Not all assignments are inserted into the graph.
        void insertAssignment(AssignmentPtr);

        // This function does not update the status of the graph.
        bool checkCompatibility(AssignmentPtr);

        // The function is responsible for optimizing the graph. 1. Switching to the internal assignments (more assignments will be used to for insertion into the graph) recorded in the pool
        void resolveConflicts(AssignmentPool& pool);

        // Removing the assignment from the graph and the status table.
        void removeAssignment(AssignmentPtr);

        // Append the assignment
        // 1. insert the assignment into the pool;
        // 2. modify the assignment status.
        void addAssignment(AssignmentPtr, AssignmentPool&);

        // Print the map.
        friend ostream& operator<<(ostream&, const FullMap&);

        // TBD: the output of the modification distribution.

    private:
        AssignmentPtr _empty_node;
        AssignmentPtr _full_node;

        std::string _mod_type;

        // The status table records the conflicting assignments.
        std::vector<StatusItem> _status_table;
        
    };
}
#endif /* GAG_FULLMAP_H */