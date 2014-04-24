/********************************************************************
	created:	2013/04/23
	created:	23:4:2013   21:45
	filename: 	GeneralCleavage.h
	file path:	GAGPL\GAGLIBRARY
	file base:	GeneralCleavage
	file ext:	h
	author:		Han Hu
	
	purpose:	General cleavage is for the position of modification.
						It only cares about the neighbor modification sites of 
						the cleavage.
*********************************************************************/

#ifndef GAG_GENERALCLEAVAGE_H
#define GAG_GENERALCLEAVAGE_H

#include <boost/shared_ptr.hpp>
#include <GAGPL/CHEMISTRY/Modification.h>
#include <GAGPL/FRAGMENTATION/Fragmentation.h>

namespace gag
{
	// For general cleavage, only the neighbor modification site(s) matter.

	// Basic functions: 
	// 1. sort all the nodes by uniqueness (uniqueness of general cleavages);
	// 1.a exp_mass => multiple nodes => set of general cleavages.
	// 1.b It seems that no need to group by one of the general cleavages.

	struct CleavageItem
	{
		// Constructor. The information stored in node item will be used for updating the cleavage_set.
		CleavageItem(const double& mass, const ModificationSites& mod_sites, const int& num, const int& clv_num, const std::string& type)
			: experimental_mass(mass), mod_sites(mod_sites), mod_num(num), clv_num(clv_num), general_type(type), estimated_mod_num((double)mod_num), uniqueness_confidence(0.0), intactness_confidence(1.0), background_confidence(1.0), repeat_num(1)
		{}
		
		CleavageItem(const ModificationSites& mod_sites)
			: mod_sites(mod_sites), uniqueness_confidence(0.0), intactness_confidence(1.0), background_confidence(1.0), repeat_num(1)
		{}
		// New version of assignment constructor.
		CleavageItem(const ModificationSites& mod_sites, const int& mod_num)
			: mod_sites(mod_sites), mod_num(mod_num), estimated_mod_num((double)mod_num), uniqueness_confidence(0.0), intactness_confidence(1.0), background_confidence(1.0), repeat_num(1)
		{}

		ModificationSites mod_sites;

		double experimental_mass;

    // The mod number specified by the assignment.
		int mod_num;
    // The mod number adjusted by the confidence.
    double estimated_mod_num;

		//int uniqueness_value;

		// New version of confidence values. Both are assumed to be 0.
		double uniqueness_confidence;
		double intactness_confidence;
    double background_confidence;

    // New version: the number of members in the same cluster.
    int repeat_num;

		int clv_num;

		std::string general_type; // Cross-ring or glycosidic-bond cleavage: it should be considered as the different observations dependent on the cleavage type.

		// New version of fragment types. Extended to support internal fragments.
		//std::vector<std::string> type_vec;
		
		//double confidence;
		//double weight;

		inline bool operator<(const CleavageItem& clv) const
		{
			return mod_sites < clv.mod_sites || (mod_sites == clv.mod_sites && general_type < clv.general_type);
		}

		inline bool isGlycosidicBondCleavage() const
		{
			return general_type == "G";
		}

		inline bool isCrossringCleavage() const
		{
			return general_type == "C"; 
		}

		inline bool isTerminalCleavage() const
		{
			return general_type == "G" || general_type == "C";
		}

    inline bool isInternalCleavage() const
    {
      return general_type == "GG" || general_type == "GC" || general_type == "CG" || general_type == "CC" || general_type == "I";
    }

    inline bool isSpecifiedType(std::string frag_type) const
    {
      return general_type == frag_type;
    }

    inline std::string getGeneralType() const
    {
      return general_type;
    }

        // New version: stop using the internal variable for confidence.
		inline double getConfidenceValue() const
		{
			return background_confidence;
		}

    inline void addRepeatNumber(int num)
    {
      repeat_num = num;
    }
    inline int getRepeatNumber() const
    {
      return repeat_num;
    }

    inline void setUniquenessValue(double uni_value)
    {
      uniqueness_confidence = uni_value; background_confidence *= uni_value;
    }

    inline double getAdjustedModificationNumber() const
    {
      return this->getConfidenceValue() * (double)mod_num;
    }

	};

	typedef boost::shared_ptr<CleavageItem> CleavagePtr;

	struct general_cleavage{};
	struct fragment_type{};
	struct exp_mass{};
	struct mod_num{};
	struct unique_num{};
  struct uniqueness_confidence{};
  struct intactness_confidence{};

	typedef multi_index_container<
		CleavageItem,
		indexed_by<
			ordered_non_unique<
				tag<general_cleavage>, member<CleavageItem, ModificationSites, &CleavageItem::mod_sites>
			>,
			ordered_non_unique<  // Sort the explanations by mass.
				tag<exp_mass>, member<CleavageItem, double, &CleavageItem::experimental_mass>
			>,
      ordered_non_unique<
      tag<uniqueness_confidence>, member<CleavageItem, double, &CleavageItem::uniqueness_confidence>,std::greater<double>
      >,
      ordered_non_unique<
          tag<intactness_confidence>, member<CleavageItem, double, &CleavageItem::intactness_confidence>, std::greater<double>
      >,
			ordered_non_unique<
				tag<fragment_type>, member<CleavageItem, std::string, &CleavageItem::general_type>
			>
		>
	> CleavageContainer;

	typedef CleavageContainer::index<general_cleavage>::type CleavageByModificationSites;
	typedef CleavageContainer::index<exp_mass>::type CleavageByExperimentalMass;
  typedef CleavageContainer::index<uniqueness_confidence>::type CleavageByUniquenessConfidence;
  typedef CleavageContainer::index<intactness_confidence>::type CleavageByIntactnessConfidence;
	typedef CleavageContainer::index<fragment_type>::type CleavagebyType;

}


#endif /* GAG_GENERALCLEAVAGE_H */