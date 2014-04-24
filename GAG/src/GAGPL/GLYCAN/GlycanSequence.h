/*
 * =====================================================================================
 *
 *       Filename:  GlycanSequence.h
 *
 *    Description:  The class for building glycan chain.
 *
 *        Version:  1.0
 *        Created:  05/ 3/2012  3:49:50 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef  GAG_GLYCANSEQUENCE_INC
#define  GAG_GLYCANSEQUENCE_INC

#include <algorithm>
#include <GAGPL/GLYCAN/Branch.h>
#include <GAGPL/CHEMISTRY/Modifier.h>
#include <GAGPL/GLYCAN/GlycanComposition.h>
#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/foreach.hpp>

namespace gag
{
	typedef boost::bimap< boost::bimaps::set_of<size_t, std::less<size_t> >, boost::bimaps::multiset_of<size_t> > BranchMap;
	typedef BranchMap::left_map::const_iterator left_const_iterator;
	typedef BranchMap::right_map::const_iterator right_const_iterator;

	// start and end of the connecting branch ids.
	typedef std::map<size_t, size_t> BranchDescendants;
	// Current branch id => first branch id on the subtree;
	typedef std::map<size_t, size_t> BranchOrigin;
	typedef std::vector<std::pair<std::string, size_t> > SequenceCode;
	

	class GlycanSequence: public Unit
	{	

		private:
			// The basic idea for designing the structure of glycan sequence is that the backbone structure should be separated from the modification part.
			std::vector<Branch> glycan_sequence;
			
			// The modification type and corresponding positions and status (true means available). Notice that the positions might be dynamically modified to fit the fragment. The content should be controlled by class Modifier.
			// ModificationMap mod_pos_map;

			// The total number of modification which should be used as a constraint for mod_position.
			std::map<std::string, int> mod_constraint;

			Modifier mod_assistant;

			// The left value should be the NRE branch and the right value should be the RE branch.
			BranchMap branch_links;
			// The key is the target branch id, the value is lowest children id, the range of the children ids is [lowest, branch_id -1].
			BranchDescendants children;
			BranchOrigin originmap;

		public:
	
			GlycanSequence() {}
			// The sequence can either be built by composition or by code.
			// When built by composition, the modification position is not specified.
			void buildByGAGComposition(const GlycanComposition& glycan_compo, std::string = "GlcA");
			// When built by code, the position is specified.
			void buildByGAGCode(const std::string code);
			

			// When a new branch is added, the masses of the corresponding 
			// monosaccharide unit will be modified.
			void addBranch(Branch&);

			void addBranchLink(const size_t& re_id, const size_t& nre_id);
			
			void addBranchLinks(const std::set<size_t>& re_ids, const size_t& nre_id);

			void updateChildrenIDs(const size_t branch_id);

			inline std::vector<Branch>& getBranches()
			{
				return glycan_sequence;
			}

			inline size_t getBranchSize()
			{
				return glycan_sequence.size();
			}

			inline Branch& getBranchByID(const size_t id)
			{
				return glycan_sequence.at(id);
			}

			inline BranchMap& getBranchLinks()
			{
				return branch_links;
			}
			
			// The site id is based on carbon id.
			inline Site& getInternalSite(const ModificationPosition& mp)
			{
				return this->getBranchByID(mp.getBranchID()).getUnitByID(mp.getMonosaccharideID()).getInternalSiteByCarbonID(mp.site_id);
			}

			void update();

			// The second value is used to evaluate if this is a leaf branch.
			std::pair<size_t, size_t> getDescendantBranchIDs(const size_t& branch_id);
			size_t getOriginBranchID(const size_t& branch_id);
		
			Composition getSubComposition(const size_t& id1, const size_t& id2);

			Composition getSubTreeComposition(const size_t& branch_id);
			Composition getTreeComposition(const size_t& branch_id);
			
			// This method only works for GAG.
			SequenceCode getGAGSequenceCode();
			std::string getGAGSequenceCodeString();

			// The function is used in fragment calculation.
			// The change of ModificationMap has to go through class Modifier.
			//inline const ModificationMap& getModificationMap() const
			//{
			//	return mod_assistant.mod_map;
			//}
			//// Notice that this is a copy.
			inline Modifier getModificationMap()
			{
				return mod_assistant;
			}
			inline ModificationByPosition& getModificationSetByPosition()
			{
				return mod_assistant.getModificationContainer().get<0>();
			}
			ModificationSites getModificationSitesBySymbol(const std::string& mod_symbol);
			ModificationSites getModificationSitesBySymbol(const std::string& mod_symbol, int status);

			// Functions for modification.
			friend class GlycanComposition;

			void initializeModificationSites(const std::string& mod_symbol);
			
			// The code of status: 
			//  1 -- available
			//  0 -- occupied
			// -1 -- unavailable (occupied or connected by other factors)
			void modifyModificationStatus(const std::string& mod_symbol, const ModificationPosition& pos, int status = 0);
			// Locate the site position. It is used only when the site information is specified. Notice that the modification status table is also updated correspondingly. If the modification is not relavent to any of the status table. Use the function addModification.
			void updateModification(const std::string& mod_symbol, ModificationPosition& pos);
			// Locate the functional group. When the position is specified, the position of modification on the ring will be used irrespective of the information stored in pos. Use it cautiously.
			void updateModification(const std::string& mod_symbol, ModificationPosition& pos, size_t site_id);

			// Useful for derivatization.
			void addModification(const std::string& mod_symbol, ModificationPosition& pos);

			std::set<std::string> getModificationTypes();

			int getModificationConstraint(const std::string& mod_symbol);
			void addModificationConstraint(const std::string& mod_symbol, int number);

			void printStructure();

   ModificationSites getComplementaryModificationSites(const ModificationSites& ms, std::string mod_symbol);
			
	};

	typedef boost::shared_ptr<GlycanSequence> GlycanSequencePtr;
}

#endif   /* ----- #ifndef GAG_GLYCANSEQUENCE_INC  ----- */
