#ifndef  GAG_BRANCH_INC
#define  GAG_BRANCH_INC

#include <algorithm>
#include <GAGPL/GLYCAN/Linkage.h>
#include <GAGPL/GLYCAN/Monosaccharide.h>

namespace gag
{
	class Branch : public Unit
	{
		private:
			size_t branch_id;
			
			std::vector<Monosaccharide> mono_chain;
			std::vector<Composition> re_extension;

			// For glycan, there will be only one reducing end.
			
			// The mono_id and the linkage.
			std::vector<Linkage> links;

		public:
			Branch(const size_t id)
				: branch_id(id), mono_chain(), links()
			{}
			Branch(){}

			void addUnit(Monosaccharide& mono_unit);
			
			void addLinkage(const Linkage& link);

			void addExtension(Composition& compo);

			std::vector<Linkage> getNeighborLinks(const size_t mono_id);

			inline size_t getBranchID() const
			{
				return branch_id;
			}

			inline std::vector<Monosaccharide>& getGlycanChainUnits()
			{
				return mono_chain;
			}
			inline const std::vector<Monosaccharide>& getGlycanChainUnits() const
			{
				return mono_chain;
			}

			inline size_t getUnitNum()
			{
				return mono_chain.size();
			}
			inline std::vector<Linkage>& getLinkages()
			{
				return links;
			}
			inline size_t getExtensionNum()
			{
				return re_extension.size();
			}

			Composition getExtensionComposition();

			inline Monosaccharide& getUnitByID(const size_t id)
			{
				return mono_chain.at(id);
			}
			
			// Calculate the mass and composition. The correction of the terminal mass will be considered.
			void update();

			// Useful for calculating fragment mass.
			Composition getSubComposition(const size_t start, const size_t end);
		
			void printStructure();
	};
}

#endif   /* ----- #ifndef GAG_BRANCH_INC  ----- */
