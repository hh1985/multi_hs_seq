/********************************************************************
	created:	2013/04/24
	created:	24:4:2013   10:25
	filename: 	SequencePrediction.h
	file path:	GAGPL\GAGLIBRARY
	file base:	SequencePrediction
	file ext:	h
	author:		Han Hu
	
	purpose:	Predict the original sequence with multiple methods.
*********************************************************************/

#ifndef GAG_SEQUENCEPREDICTION_H
#define GAG_SEQUENCEPREDICTION_H

#include <GAGPL/GAGLIBRARY/GeneralCleavage.h>
#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/SPECTRUM/MonoPeak.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/GAGLIBRARY/CleavageTree.h>

namespace gag
{
	// The modification sites, merged mod number and cleavage type.
	typedef std::map<ModificationPosition, double> ModificationDistribution;

  void printModificationDistribution(const ModificationDistribution& dist, const std::string& mod_name);


	// The modification symbol and corresponding distribution.
	typedef std::map<std::string, ModificationDistribution> ModificationPattern;

  //typedef std::map<std::string, ModificationSites> ModificationSequence;

 

	// The space for calculating uniqueness of fragments.
	typedef CleavageContainer CleavageSpace;
	typedef CleavageContainer CleavageElite;
	
	// The represents for building the distribution.
  //typedef std::multimap<MonoPeakPtr, NodeItem> AssignmentSpace;
  //struct AssignmentElite {
  //  AssignmentSpace elite;
  //  double weight;
  //  
  //};
  struct AssignmentMember
  {
    MonoPeakPtr mono_pk_ptr;
    NodeItem node;
    double weight;

    AssignmentMember()
      : weight(1.0) {}
    AssignmentMember(MonoPeakPtr pk, NodeItem& node_item, double cluster_size)
      : mono_pk_ptr(pk), node(node_item), weight(cluster_size)
    {}

    inline bool operator<(const AssignmentMember& member) const 
    {
      return mono_pk_ptr < member.mono_pk_ptr;
    }
  };

  typedef std::multiset<AssignmentMember> AssignmentGroup;
  typedef std::multimap<MonoPeakPtr, NodeItem> AssignmentSpace;
  typedef std::map<std::string, std::set<std::pair<ModificationSites, int>>> AssignmentClassification;
	//typedef std::map<std::pair<MonoPeakPtr, NodeItem>, double> AssignmentElite;
	

	class SequencePrediction
	{
	public:
		// Constructors.
		SequencePrediction(std::multimap<MonoPeakPtr, NodeItem>& matched_results, GlycanSequencePtr seq) : gs(seq), param(Param::Instance())
		{
			this->build(matched_results);
		}

		// Estimate the distribution of the modification across each site. If missing is true, which means the modification site might be missing; otherwise, it will fix to the backbone.
		// If missing is true, consider about cleaning the fragments with missing modification away. (grouping)
		// For identification of Ac, consider only fragments from glycosidic-bond.
		// For identification of SO3, consider the terminal cleavages first to draw the framework and then the internal cleavages to tune the distribution of sulfate groups.
		// New: The interface for estimating the modification distribution.

		// Convert the data into a unique mass list.
		//std::set<double> getExperimentalMassSet() const;

		// Get uniqueness value by mass
		//size_t getNumberOfAssignments(const double& mass) const;
		//size_t getNumberOfGeneralCleavages(const double& mass) const;

		// The function can get complementary modification sites which will be used for searching recorded number of modification events.
		
		// Get the background modification number.
		//double getModificationNum(const ModificationSites& mod_sites);

    // Print the internal structure of the distribution.
    void printModificationDistribution(std::string mod_symbol);
	  void printModificationDistribution(const ModificationDistribution& dist) const;
    ModificationDistribution getModificationDistribution(std::string mod_symbol)
    {
      return mod_pattern[mod_symbol];
    }

    // Naive method 1 for calculating the coverage of mass values for each candidate sequence.
    //double calculateCoverage(std::multimap<MonoPeakPtr, NodeItem>& filtered_data, int cleavage_num = 2);

     // Naive method 2 for calculating the number of golden pair. Notice to remove redundant pairs.  Basically, if the modification sites of two assignments are complementary to each other, we define this as a pair.
    //double calculateGoldenPairNum(std::multimap<MonoPeakPtr, NodeItem>& filtered_data);
    
    // Statistical method for calculating the merged cost of HS-SEQ.
    double calculateCost(const ModificationSequence& seq);

    double calculateMergedCost(const ModificationSequence& seq);

    // This function will get subset of the data based on specified sequence code.  It should be run before any sequencing methods.
    // status specifies if the data should be dealt for merged sequence.
  	//std::multimap<MonoPeakPtr, NodeItem> filterData(std::multimap<MonoPeakPtr, NodeItem>& matched_results, const ModificationSequence& seq, bool status);

	// New version of functions which match exactly the description on the manuscript.
	private:
    // Organize the procedure of calculating the modification distribution.
    void build(std::multimap<MonoPeakPtr, NodeItem>& matched_results);

    ModificationSites getComplementaryModificationSites(const ModificationSites& ms);

		// Low level processing of assignments.  Assignments differing only in neutral loss, hydrogen shift and sulfate difference should be grouped.
		//AssignmentElite groupAssignments(const AssignmentSpace& matched_results);

    AssignmentGroup groupAssignments(const AssignmentSpace& matched_results, bool loss = true);
    
    // Convert the data into non-redundant observation.
    AssignmentSpace groupAssignmentsByFragmentType(const AssignmentSpace& matched_results, bool missing = true);

    // For Ac, the function accumulates evidence for its locations;
    // For SO3, the function gets the one with top sulfur number.
    // The cleavage items for grouping should have the same uniqueness values.
    AssignmentSpace groupAssignmentsByModificationSites(AssignmentSpace& evidence, bool missing = true);

		// Selection of assignments fitting specific fragment type requirements

    // Decide if the assignment is a qualifed type.
    // Only glycosidic-bond fragment should be counted.
    bool isQualifiedAssignment(const Fragment& frag) const;

		// Uniqueness value of an assignment in the given pool.
		void estimateUniquenessValue(CleavageElite& clv_elite, CleavageSpace& clv_pool);
    //double estimateUniquenessValue(const CleavageItem& clv_item, CleavageSpace& clv_pool);

    // Uniqueness of assignment in the pool
    double estimateUniquenessValue(AssignmentSpace::iterator ele, const AssignmentSpace& pool);

    

    double calculateSiteUniqueness(std::string type, const AssignmentClassification& assignments);
    // 1. Remove false assignments based on known distribution.
    // 2. Update modification status.
		void filterAssignments(const std::string& mod_symbol, AssignmentSpace& assignment_pool);

		// Constructing fragment graph v1. 
        // The modification distribution is modified during the graph 
        // construction process.
		//void constructFragmentGraph(CleavageElite& fragment_rep);
        // The modification distribution is modified based on accumulation of evidences from different assignments.
    //void accumulateModificationEvidence(CleavageElite& clv_elite, std::string mod_symbol = "Ac");

		// Sulfate intactness of an assignment.
		//void estimateSulfateIntactness(CleavageItem& clv, CleavageTree& clv_tree);

    // Initialize the modification distribution assuming all sites are equally distributed. Only works for the specified mod_symbol.
    void initialize();

    // Since the identification is on residue level, the distribution must reflect this trend.
    void initializeModificationDistribution(ModificationDistribution& mod_dist);

    // The functional groups have equal chance to exist across all sites.
    void equalInitialization(ModificationDistribution& mod_dist);

		//void updateCleavageTree(CleavageItem& clv, CleavageTree& clv_tree);

		// Convert NodeItem to abstract cleavage type.
		CleavageItem convertNodeItem(const MonoPeak& mono_pk, const NodeItem& node);

    void averageBySite(ModificationDistribution& mod_dist, const ModificationSites& mod_sites, double total_num);
    void averageByResidue(ModificationDistribution& mod_dist, const ModificationSites& mod_sites, double total_num);
		/* 
		// TBD: Constructing fragment graph v2.
		CleavageTree constructFragmentGraph(AssignmentElite& fragment_rep);
		ModificationDistribution generateDistribution(CleavageTree& clv_tree); 
		*/

    // New version of estimating modification distribution.
    // If separate is false, cross-ring cleavage and gb-cleavage are treated equally.
    bool generateModificationDistribution(std::string mod_symbol, const AssignmentSpace& assignment_pool, bool separate = true);

    // The latest and oldest version of the function.
    bool generateModificationDistributionV2(std::string symbol, const AssignmentSpace& pool, bool missing = false);

    // Based on graph.
    bool generateModificationDistributionV3(std::string symbol, const AssignmentSpace& pool, bool missing = false);

    // Cleavage type and corresponding assignments (abstract form)
    //typedef std::map<std::string, std::set<std::pair<ModificationSites, int>>> AssignmentClassification;

    // The unique assignments are used for calculating uniqueness as well as estimate the background.
    AssignmentClassification generateUniqueAssignments(const std::pair<AssignmentSpace::const_iterator, AssignmentSpace::const_iterator>& assign_range);

    // dist is the background distribution.
    ModificationDistribution updateDistributionV2(const AssignmentClassification& uni_set, const ModificationDistribution& dist);

    // Calculate p(x|D) -- likelihood.  Estimate the possibility of each choice based on prior distribution.
    ModificationDistribution updateDistributionV3(const AssignmentClassification& uni_set, const ModificationDistribution& dist);

    // For given parent and child node, update the modification distribution defined between them.
    // The update is directly applicable to mod_pattern, and mod_symbol is used to retrieve site information. 
    void updateLocalDistribution(ModificationDistribution& mod_dist, CleavageItem& clv, CleavageNodePtr child, CleavageNodePtr parent);

    void updateLocalDistribution(CleavageItem& clv, CleavageTree& clv_tree);

    // The method is used to distinguish the site difference which assist further evaluation. Normalization to 0 and 1.
    void scalingModificationDistribution();

    double getAccumulatedDensity(const ModificationDistribution& mod_dist, const ModificationSites& child_sites, const ModificationSites& parent_sites);

    ModificationSites getDiffSites(const ModificationSites& child_sites, const ModificationSites& parent_sites);

    // The universal solution for uniqueness-based distribution estimation.
    void updateGlobalDistribution(CleavageElite& clv_rep, CleavageTree& last_tree);

    // Check the cleavage tree for locating the assignment, estimate the distance from the background number.
    void estimateBackgroundConfidence(std::pair<CleavageByUniquenessConfidence::iterator, CleavageByUniquenessConfidence::iterator> uni_pair, CleavageTree& clv_tree);

    void estimatePeerConfidence(std::pair<CleavageByUniquenessConfidence::iterator, CleavageByUniquenessConfidence::iterator> uni_pair);

    // Internal function for estimating the confidence value from the graph structure. The graph can be current graph including clv, or it can be the background graph.
    // If complementary == 1, the complementary node of clv will be checked. Otherwise it will be simply skipped.
    double inferFromGraph(CleavageItem& clv, CleavageTree& clv_tree, bool complemenatry = true);
       
    // This module will correct the sulfation likelihood within residue using biosynthetic rule: 2-N > 6-O > 3-O.
    void correctSulfateDistribution();
    
    

	private:
		// 1. Filter for assignments.
		// 2. Calculation of uniqueness values for qualified assignments.
		//CleavageContainer generateCleavageContainer(bool missing = true);

		// This function returns a bool value deciding if the modification site is a type of NRE-intact cleavage. If the sites include the sites on the first ring, return yes, otherwise, return false.
    bool isNRECleavage(const ModificationSites& ms) const;
		// RE-intact cleavage. If the sites include the sites on the last ring, return yes, otherwise, return false.
    bool isRECleavage(const ModificationSites& ms) const;

		// Terminal cleavage in a general sense.
    inline bool isGeneralTerminalCleavage(const ModificationSites& ms) const
		{
			return (isNRECleavage(ms) || isRECleavage(ms)) && !isIntactCleavage(ms);
		}
		inline bool isTerminalCleavage(const CleavageItem& clv) const
		{
			return clv.clv_num == 1;
		}
		inline bool isInternalCleavage(const CleavageItem& clv) const
		{
			return clv.clv_num == 2;
		}

		inline bool isIntactCleavage(const ModificationSites& ms) const
		{
			return isNRECleavage(ms) && isRECleavage(ms);
		}

    

		inline bool containSubset(const ModificationSites& ms1, const ModificationSites& ms2) const
		{
			return std::includes(ms1.begin(), ms1.end(), ms2.begin(), ms2.end());
		}


	private:
		//std::multimap<double, NodeItem> data;
		GlycanSequencePtr gs;
		
    Param& param;
    // Stores current modification symbol.
		std::string mod_symbol;

    // The cleavage type: G or C.
    std::string mod_status;

		// Storing the overall modification distribution.
		ModificationPattern mod_pattern;

    std::vector<std::string> type_vec;

    std::map<std::string, double> type_score;
		
	};
}

#endif /* GAG_SEQUENCEPREDICTION_H */