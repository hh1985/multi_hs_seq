#ifndef GAG_LIBRARYFOREST_H
#define GAG_LIBRARYFOREST_H

#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <algorithm>

namespace gag
{
	typedef std::vector<ModificationPosition> ModificationSiteList;
	typedef std::vector<ModificationSiteList> ModificationPatternList;
	typedef std::pair<ModificationSiteList, ModificationSiteList> GlobalModPattern;

	// The class LibraryForest is created to allow exhaustive enumeration of modification sites.
	class LibraryForest
	{
	private:
		GlycanSequence& gag_seq;
		std::map<std::string, size_t> mod_list;
		std::vector<GlobalModPattern> mod_pattern_list;

	protected:
		// Generate all theoretical candidate modification patterns.
		ModificationPatternList guessModificationPositions(ModificationSiteList& sites, size_t num);
		
	public:
		LibraryForest(GlycanSequence& backbone) :
			gag_seq(backbone), mod_pattern_list()
		{
		}

		// Update the gag_sequence;
		LibraryTree applyModification(GlobalModPattern& global_pattern);
		//inline void addLibraryTree(LibraryTree& tree)
		//{
		//	forest.push_back(tree);
		//}
		inline size_t getTreeNum()
		{
			return mod_pattern_list.size();
		}
		inline void addModList(const std::string& mod, const size_t num)
		{
			mod_list.insert(std::make_pair(mod, num));
		}
		//LibraryTree getLibraryTreeByModPattern(GlobalModPattern& global_pattern);
		//inline std::vector<LibraryTree>& getLibraryTrees()
		//{
		//	return forest;
		//}
		// Generating all chemically possible sites for acetation.
		ModificationSiteList getCandidateAcetationSites();
		// Generating all chemically possible sites for sulfation. 
		ModificationSiteList getCandidateSulfationSites();
		
		//void printGAGSequence(GlycanSequence);
		// This algorithm stores the information into mod_pattern_list;
		void generateTreeByComposition(const size_t& num_s, const size_t& num_ac);

		inline std::vector<GlobalModPattern>& getGlobalModification()
		{
			return mod_pattern_list;
		}

	};
}

#endif /* GAG_LIBRARYFOREST_H */