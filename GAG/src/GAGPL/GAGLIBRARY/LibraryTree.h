#ifndef GAG_LIBRARYTREE_H
#define GAG_LIBRARYTREE_H

#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <GAGPL/SPECTRUM/MonoPeak.h>
#include <map>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/member.hpp>

namespace gag
{
	
	// A combination of neutral loss. 
	typedef std::map<std::string, int> Satellite;

	//bool operator<(const Satellite& sate1, const Satellite& sate2)
	//{
	//	return sate1.size() < sate2.size();
	//}

	using namespace ::boost;
	using namespace ::boost::multi_index;

	class NodeItem: public Unit
	{
	private:
		FragmentPtr fg; // Consider about using shared_ptr.
		Satellite sate;
		
	public:
		NodeItem(FragmentPtr fg)
			: fg(fg), Unit(fg->getComposition())
		{}
    NodeItem() {}
		NodeItem(const NodeItem& node);

		NodeItem& operator=(const NodeItem& node);

		// process the neutral loss and number information.
		bool appendMassShift(const std::string& loss_type, int num);

		inline Satellite& getCompositionShiftList()
		{
			return sate;
		}
		inline const Satellite& getCompositionShiftList() const
		{
			return sate;
		}
		inline FragmentPtr getFragment()
		{
			return fg;
		}
		inline const FragmentPtr getFragment() const
		{
			return fg;
		}

		int getModificationNum(const std::string& compo_string) const;

		inline std::string getCleavageType() const
		{
			return fg->getCleavageType();
		}
    // Stupid function.
    inline std::string getCleavageClass() const
    {
      return fg->getCleavageCollection().begin()->first;
    }
    
    inline FragmentPosition getFragmentPosition() const
    {
      return fg->getCleavageCollection().begin()->second;
    }

    inline std::string getGeneralType() const
    {
      return fg->getGeneralType();
    }

		inline size_t getCleavageNum() const
		{
			return fg->getCleavageNum();
		}

		std::string getCompositionShift() const;

		// if format is "regular", the results are separate by multiple lines with title added. If it is "compact", the results are compressed into a single line.
		void printNodeItem(const std::string format = "regular", const size_t clv_num = 2) const;
	
	};

	typedef boost::shared_ptr<NodeItem> NodePtr;

	struct theo_mass{};
	struct node_frag{};
	typedef multi_index_container<
		NodeItem,
		indexed_by<
			ordered_non_unique<  // Sorted by theoretical mass.
			tag<theo_mass>, const_mem_fun<Unit, double, &NodeItem::getMass>
			>,
			ordered_non_unique<  // Sorted by fragment type.
			tag<node_frag>, const_mem_fun<NodeItem, std::string, &NodeItem::getCleavageType>
			>
		>
	> TreeContainer;

	typedef TreeContainer::index<theo_mass>::type TreeByMass;
	typedef TreeContainer::index<node_frag>::type TreeByCleavage;

	// The class is designed for gag. To be implemented for complicate structure.
	class FragmentTree
	{
	public:

		FragmentTree(GlycanSequencePtr seq)
			: glycan_seq(seq), frag_table(FragmentationTable::Instance()), param(Param::Instance())
		{
			this->build();
		}

		inline void addNewNode(NodeItem& node)
		{
			node_pool.insert(node);
		}
		// Get items within error range.
		void getNodeItemsByMass(MonoPeakPtr pk, std::multimap<MonoPeakPtr, NodeItem>& node_map);
		std::multimap<MonoPeakPtr, NodeItem> searchLibrary(const std::set<MonoPeakPtr>& pk_list);

		//inline size_t calculateUniqueness(double mass)
		//{
		//	return this->getNodeItemsByMass(mass).size();
		//}
		// The input can be a specified sulfation pattern.

		inline size_t getNodeSize() const
		{
			return node_pool.size();
		}
		
		void printLibrary();

		void exportLibrary(const std::string& filename);

	private:
		// Generate all the fragments.
		void build();
		void addCleavage();
		// The partially processed fragment and location.
		void addRestCleavage(const CleavageCollection& cc, const size_t& cur);

		// The function will result in the expansion of node_pool.
		void processFragment(const CleavageCollection& cc, MassLossWindow mlw);
		void processFragment(FragmentPtr fg, MassLossWindow& mlw, size_t degree, Satellite sate);

		// Map the fragments into sequence and identify the possibility of each modification sites. Notice that the way of estimating acetate group is different from estimating sulfate groups.
		//std::map<ModificationPosition, double> estimateAcetateGroupPosition(const std::vector<NodeItem>& node_vec);

	private:
		Param& param;
		FragmentationTable& frag_table;
		GlycanSequencePtr glycan_seq;
		TreeContainer node_pool;
	};

	void loadMonoPeakList(const std::string& filename, std::set<MonoPeakPtr>& peak_list);
}





#endif /* GAG_LIBRARYTREE_H */