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

#include <GAGPL/GAGLIBRARY/AssignmentPool.h>
#include <GAGPL/GLYCAN/GlycanSequence.h>

namespace gag
{
    using namespace std;

    class FullMap
    {
    public:
        // Constructor.
        FullMap(GlycanSequencePtr seq, const std::string& mod, set<BackbonePtr>& backbone_set)
            : seq(seq), _mod_type(mod), _bone_set(backbone_set)
        {
          this->initialize();
          this->connectBackboneSet();
          //this->exploreCompatibility(_empty_node);
        }

        // The function is responsible for optimizing the graph. 1. Switching to the internal assignments (more assignments will be used to for insertion into the graph) recorded in the pool
        void resolveConflicts(AssignmentPool& pool);

        Backbone& getBackbone(AssignmentPtr assignment);
        set<Backbone> getBackboneNeighbor(const Backbone& bone);
        set<Backbone> getBackboneNeighbor(AssignmentPtr assignment);

        bool exploreUpperNodes(BackbonePtr cur, BackbonePtr child_node, BackbonePtr parent_node);
        // Decide if any of the path is OK for insertion.
        void exploreEntryPoint(BackbonePtr cur, BackbonePtr child_node);

        // Check the compatibility between assignments of each backbone pair.
        void exploreCompatibility(AssignmentPool& pool);

        bool checkCompatibility(BackbonePtr small_bone, BackbonePtr large_bone);
        bool checkCompatibility(BackbonePtr small_bone, BackbonePtr large_bone, int small_num, int large_num, int diff_size);

        // Print the map.
        friend ostream& operator<<(ostream&, const FullMap&);

        bool isNRECleavage(BackbonePtr cur) const;
        bool isRECleavage(BackbonePtr cur) const;

        // TBD: the output of the modification distribution.
    private:

      // Iterate over all bones, update the connections between bones
      void connectBackboneSet();
      void initialize();

    private:
        BackbonePtr _empty_node;
        BackbonePtr _full_node;

        std::string _mod_type;

        set<BackbonePtr>& _bone_set;

        GlycanSequencePtr seq;
        // The status table records the conflicting assignments.
        // the key is the target, and the value is the neighbor.
        multimap<AssignmentPtr, AssignmentPtr> _status_table;
        
    };
}
#endif /* GAG_FULLMAP_H */