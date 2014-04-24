/*
 * =====================================================================================
 *
 *       Filename:  Fragmentation.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  4/23/2012 4:39:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef  GAG_FRAGMENTATION_H_INC
#define  GAG_FRAGMENTATION_H_INC

#include <GAGPL/FRAGMENTATION/FragmentationTable.h>
#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/MISC/Param.h>
#include <boost/shared_ptr.hpp>

using namespace param;

namespace gag
{
  //typedef std::vector<std::vector<Branch> > Fragments;
	enum re_cleavage {RE_INTACT, B, C, A, RE_N};
	enum nre_cleavage {NRE_INTACT, Y, Z, X, NRE_N};

	struct FragmentPosition
	{
		// The branch ID should follow the leaf-to-root style.
		//size_t branch_id;
		//size_t mono_id;

		MacroPosition mac_pos;
		// Ring ID. Please be careful of the conversion from Carbon ID to
		// Ring ID.
		size_t xring_first;
		size_t xring_second;

		inline size_t getBranchID() const { return mac_pos.branch_id; } 
		inline size_t getMonosaccharideID() const { return mac_pos.mono_id; }
		inline void setBranchID(const size_t& b_id) 
		{
			mac_pos.branch_id = b_id;
		}
		inline void setMonosaccharideID(const size_t& m_id)
		{
			mac_pos.mono_id = m_id;
		}

    inline bool operator<(const FragmentPosition& frag_pos) const
    {
      return mac_pos < frag_pos.mac_pos || (mac_pos == frag_pos.mac_pos && xring_first < frag_pos.xring_first) || (mac_pos == frag_pos.mac_pos && xring_first == frag_pos.xring_first && xring_second < frag_pos.xring_second);
    }
		FragmentPosition(const size_t& b_id, const size_t& m_id, const size_t& x1, const size_t& x2)
			: mac_pos(b_id, m_id), xring_first(x1), xring_second(x2) {}
	};

	// Cleavage type and corresponding fragmentation position.
	// A fragment may consist of a couple of cleavage types and their corresponding fragmentation positions.
	typedef std::multimap<std::string, FragmentPosition> CleavageCollection;
	
	

	class Fragment: public Unit, public Modifier
  {

    private:
			

			// The item should be the same as the one stored in the glycan sequence.
			// with the difference of their status.
			//ModificationMap sub_mod_map;

			//FragmentationTable& ft;
			Param& param;

			// Actions related to modification should be processed by mod_assistant.
			// Modifier mod_assistant;
	//public:
			// All calculation is based on position.
			GlycanSequencePtr glyco_seq;
			// Cleavage type and site.
			CleavageCollection cleavage_sites;
		
		protected:
			// A generic function which can be called by generate[X]Type() ([X] = X,Y,Z,A,B,C here). Using information from cleavage_sites.
			Composition getFragmentComposition(const FragmentPosition& fp, bool cw = true/* clockwise */);
					
		public:
			// Constructor.
			Fragment(GlycanSequencePtr seq)
				: glyco_seq(seq), /*ft(FragmentationTable::Instance()),*/ param(Param::Instance()), cleavage_sites(), Unit(seq->getComposition()), Modifier(seq->getModificationMap())
			{}

			Fragment(const Fragment& fg);

			Fragment() 
				: param(Param::Instance())
			{}

			inline void addSequence(GlycanSequencePtr seq_ptr)
			{
				glyco_seq = seq_ptr;
			}
			void printFragment() const;

			// Check if the cleavage is a type of internal cleavage.
			bool isInternalCleavage(const FragmentPosition& fp);

			bool isCrossringCleavage() const;
			bool isGlycosidicBondCleavage() const;
      inline bool isTerminalCleavage() const
      {
        return cleavage_sites.size()==1;
      }

      // If cross-ring cleavage, glycosidic-bond cleavage.
			std::string getGeneralType() const;

			// Repeatedly store all the fragmentation information. The modification information should be updated correspondingly.
			// type: the cleavage type.
			// fp: the branch id, mono id and cross-ring cleavages.
			// Notice that this is different from 
			void setFragmentation(const std::string& type, const FragmentPosition& fp);
			// A more friendly version, which will finish the conversion from 
			// nomenclature to internal implementation.
			void setFragmentation(const std::string& type, const size_t& index, const std::string leaf, const size_t& xring_first, const size_t& xring_second);
      // Correct version.
      void setFragmentation(const std::string& type, const size_t& index, const size_t& xring_first, const size_t& xring_second);

			void updateCleavage(const CleavageCollection& cleavage_col);

			inline const CleavageCollection& getCleavages() const
			{
				return cleavage_sites;
			}

			Fragment& operator=(const Fragment&);

			// Convert the cleavage collection information to string.
      // Notice the string contains string information.
			std::string getCleavageType() const;

      std::string getFragmentType() const;

			inline size_t getCleavageNum() const
			{
				return cleavage_sites.size();
			}

			// A generic function for fragmentation.
			//void updateFragmentByType(const std::string& type);	
			// The update process can be divided into several steps:
			// 1. Locate the cleavage position using information from Fragmentation.
			// 2. Update the mass and composition of that Branch.
			// 3. Clear the composition and mass of all related branch.
			// 4. Recalculate the overall composition and mass.
			// void update();

			// The value of flag: 1 -- available; 
			//                   -1 -- unavailable;
			//                    0 -- both.
			// By default, the function will return all the available sites.

			void cleanModificationStatus(const FragmentPosition& site, bool cw = true);

			void cleanSubTreeModificationStatus(const size_t& branch_id);
			void cleanTreeModificationStatus(const size_t& branch_id);

			void cleanChainModificationStatus(const size_t& branch_id, const size_t& m1, const size_t& m2);
			//void cleanUnitModificationStatus(const Monosaccharide& mono);
			
			// If s1 < s2, clean the site ids from s1 to s2.
			// If s1 >= s2, clean the site ids except s1 to s2.
			// Convert ring id to carbon id first and do the clean job.
			void cleanSubUnitModificationStatusByRingID(const MacroPosition& mac_pos, const size_t& s1, const size_t& s2);

			void cleanSubUnitModificationStatusByCarbonID(const MacroPosition& mac_pos, const size_t& c1, const size_t& c2);

			inline ModificationByPosition& getModificationSetByPosition()
			{
				return this->getModificationContainer().get<0>();
			}

			inline GlycanSequencePtr getGlycanSequence()
			{
				return glyco_seq;
			}
			inline CleavageCollection& getCleavageCollection()
			{
				return cleavage_sites;
			}

      // Terminal cleavage only, otherwise return -1.
      int getMonosaccharidePosition() const;
      size_t getCleavageIndex() const;

      // This function will expand current cross-ring cleavage.
      ModificationSites getExpandedModificationSites(std::string mod_symbol, int flag = 1);
      ModificationSites getReducedModificationSites(std::string mod_symbol, int flag = 1);

  };

	typedef boost::shared_ptr<Fragment> FragmentPtr;


}





#endif   /* ----- #ifndef GAG_FRAGMENTATION_H_INC  ----- */
