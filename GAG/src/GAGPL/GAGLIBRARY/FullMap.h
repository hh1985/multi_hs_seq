/********************************************************************
	created:	2014/04/29
	created:	29:4:2014   0:04
	filename: 	FullMap.h
	file path:	GAG\src\GAGPL\GAGLIBRARY
	file base:	FullMap
	file ext:	h
	author:		Han Hu
	
	purpose:	Manage backbones into a graph. 
*********************************************************************/

#ifndef GAG_FULLMAP_H
#define GAG_FULLMAP_H

#include <GAGPL/GAGLIBRARY/AssignmentPool.h>
#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/GAGLIBRARY/Backbone.h>

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
      }

    public:
      /*  Operation of assignments */

      // Insert the assignment into the graph. If there is backbone corresponding to the assignment in the graph, the assignment will be added directly into the backbone. Otherwise, a new backbone will be created.
      void addAssignment(AssignmentPtr assignment);

      // Remove the assignment from the graph, if this causes the removal of backbone, the topology of the graph will be adjusted.
      void removeAssignment(AssignmentPtr assignment, BackbonePtr bone);


      /* Operation of backbones */

      // Insert the backbone into the graph.
      void addBackbone(BackbonePtr bone);

      // Remove the backbone from the graph.
      void removeBackbone(BackbonePtr bone);

      set<BackbonePtr> getPotentialParents(AssignmentPtr assignment);
      set<BackbonePtr> getPotentialParents(BackbonePtr bone);
      
      set<BackbonePtr> getPotentialChildren(AssignmentPtr assignment);
      set<BackbonePtr> getPotentialChildren(BackbonePtr bone);

      // Get the child-parent pair, Useful only for internal cleavage. 
      set<pair<BackbonePtr, BackbonePtr>> getPotentialNeighbors(AssignmentPtr assignment);
      set<pair<BackbonePtr, BackbonePtr>> getPotentialNeighbors(BackbonePtr bone);


      /* Operation of the graph */

      bool exploreUpperNodes(BackbonePtr cur, BackbonePtr child_node, BackbonePtr parent_node);
      // Decide if any of the path is OK for insertion.
      void exploreEntryPoint(BackbonePtr cur, BackbonePtr child_node);
      
      /* Check the compatibility between assignments */
      // Iterate over the backbone container, check the compatibility between child and parent
      void screenBackbones();

      // Check the compatibility between assignments of each backbone pair.
      void exploreCompatibility(AssignmentPool& pool);

      bool checkCompatibility(BackbonePtr small_bone, BackbonePtr large_bone);
      bool checkCompatibility(BackbonePtr small_bone, BackbonePtr large_bone, int small_num, int large_num, int diff_size);
      
      void checkCompatibility();
      int calculatePressure(AssignmentPtr);

      /* Solve the incompatibility */
      set<AssignmentPtr> findSolution(AssignmentPtr assignment, AssignmentPool& pool);
      
      /* Check the status */  
      bool findNRECleavage(AssignmentPtr assignment) const;
      bool findRECleavage(AssignmentPtr assignment) const;

      /* Operator overload */
      friend ostream& operator<<(ostream&, const FullMap&);

    // TBD: the output of the modification distribution.
    private:

      // Iterate over all bones, update the connections between bones
      void connectBackboneSet();

      // Initialize the backbones for the empty_node and full_node.
      void initialize();

    private:
        BackbonePtr _empty_node;
        BackbonePtr _full_node;

        std::string _mod_type;

        set<BackbonePtr>& _bone_set;

        // Store the backbones from XYZ fragments.
        set<BackbonePtr> _re_set;

        // Store the backbones from ABC fragments.
        set<BackbonePtr> _nre_set;

        // Support the query of modification sites
        // Not good for finding neighbors.
        map<ModificationSites, BackbonePtr> _bone_map;

        GlycanSequencePtr seq;

        // The status table records the conflicting assignments.
        // the key is the target, and the value is the neighbor.
        multimap<AssignmentPtr, AssignmentPtr> _status_table;
        
    };
}
#endif /* GAG_FULLMAP_H */