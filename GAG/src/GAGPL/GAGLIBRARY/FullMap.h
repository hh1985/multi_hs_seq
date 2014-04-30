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

namespace gag
{

    class FullMap
    {
    public:
        // Constructor.
        FullMap(GlycanSequencePtr seq, const std::string& mod, set<Backbone>& backbone_set)
            : seq(seq), _mod_type(mod), _bone_set(backbone_set)
        {
          this->initialize();
          this->connectBackboneSet();
        }

        // The function is responsible for optimizing the graph. 1. Switching to the internal assignments (more assignments will be used to for insertion into the graph) recorded in the pool
        void resolveConflicts(AssignmentPool& pool);

        Backbone& getBackbone(AssignmentPtr assignment);
        set<Backbone> getBackboneNeighbor(const Backbone& bone);
        set<Backbone> getBackboneNeighbor(AssignmentPtr assignment);

        void exploreDeepNodes(BackbonePtr cur, BackbonePtr tree_node);

        // Check the compatibility between assignments of each backbone pair.
        void exploreCompatibility(BackbonePtr bone);

        bool checkCompatibility(BackbonePtr small_bone, BackbonePtr large_bone);
        bool checkCompatibility(BackbonePtr small_bone, BackbonePtr large_bone, int small_num, int large_num, int diff_size);

        // Print the map.
        friend ostream& operator<<(ostream&, const FullMap&);


        // TBD: the output of the modification distribution.
    private:
      // Iterate over all bones, update the connections between bones
      void connectBackboneSet(set<Backbone>& backbone_set);
      void initialize();

    private:
        Backbone _empty_node;
        Backbone _full_node;

        std::string _mod_type;

        set<Backbone> _bone_set;

        GlycanSequencePtr seq;
        // The status table records the conflicting assignments.
        // the key is the target, and the value is the neighbor.
        multimap<AssignmentPtr, AssignmentPtr> _status_table;
        
    };
}
#endif /* GAG_FULLMAP_H */