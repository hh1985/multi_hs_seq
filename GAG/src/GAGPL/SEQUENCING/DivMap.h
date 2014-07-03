/*!
 * \file DivMap.h
 *
 * \author Han
 * \date June 2014
 *
 */

#ifndef GAG_DIVMAP_H
#define GAG_DIVMAP_H

//#include <GAGPL/SEQUENCING/Division.h>
#include <GAGPL/SEQUENCING/DivisionCluster.h>
#include <GAGPL/SEQUENCING/AssignmentPool.h>

namespace gag
{

  struct DivPath
  {
    // Each vector represents a path.
    set<vector<DivisionPtr>> path_set;

    // For reporting the results of DivMap.
    friend ostream& operator<<(ostream& os, const DivPath& path);
  };

/*!
 * \class DivMap
 *
 * \brief The class is responsible for organizing Division objects.
 *
 * \author Han
 * \date June 2014
 */

  class DivMap
  {
  public:
    DivMap(GlycanSequencePtr gs, string mod_symbol);
    DivMap(const ModificationSites& mod_sites, int mod_num, string mod_symbol);

    //void addAssignment(AssignmentPtr assign);
    
    // Locate and insert the Division object.
    void addDivisionNode(DivisionPtr div_ptr);

    int getTotalPathNumber() const;
    int getMinimumPathNumber() const;

    DivPath getMaximumPath() const;
    DivPath getMinimumPath() const;

    // Organize the graph based on the assignment information.
    void cleanMap(const AssignmentPool& assign_dict);
 
    // Report the adjacency list.
    friend ostream& operator<<(ostream& os, const DivMap& map);

  private:
    // Create the network with only two dummy nodes.
    void initilize();

    // Check the type of the assignment
    bool qualityCheck(AssignmentPtr assign)const ;

    // Check if the div object has been stored in the map.
    //DivisionPtr checkRecord(const ModificationSites& sites, int num) const;
    bool checkCompatibility(DivisionPtr div1, DivisionPtr div2);

    // For the situation where there is connection between the child and the parent.
    void insertNode(DivisionPtr child_node, DivisionPtr parent_node, DivisionPtr div_nod);
    // For the situation where there is no connection between the child and parent.
    void connectNode(DivisionPtr child_node, DivisionPtr parent_node, DivisionPtr div_node);
    void exploreMap(DivisionPtr child_node, DivisionPtr check_node, DivisionPtr div_node);
     
  private:
    // Dummy nodes
    DivisionPtr _empty_node;
    DivisionPtr _full_node;

    // The total information.
    ModificationSites full_sites;
    int full_num;
    string mod_symbol;

    DivisionCluster _map;

  };
}
#endif /* GAG_DIVMAP_H */