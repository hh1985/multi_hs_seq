#include <GAGPL/GAGLIBRARY/FullMap.h>
#include <GAGPL/IO/SequenceReader.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>

using namespace std;
using namespace gag;
using namespace gag::io;

int main(int argc, char **argv)
{
  try
  {
    // 1. Load the sequence.
    SequenceReader seq_reader;
    string seqfile = "./data/structure/hex6.gl";
    //string pk_list = "hex6_5_NETD";

    GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
    seq_reader.readSequenceStructure(seqfile, gs);

    string mono_file = "./data/output/hex6_5_NETD_mono_list.txt";

    set<MonoPeakPtr> mono_set;
    seq_reader.readMonoPeakList(mono_file, mono_set);

    // 2. Create the library

    FragmentTree glyco_tree(gs);
    std::cout << "Get data size: " << glyco_tree.getNodeSize() << "\n";
    std::multimap<MonoPeakPtr, NodeItem> node_map = glyco_tree.searchLibrary(mono_set);

    std::cout << "Assignment data size: " << node_map.size() << "\n";

    // 3. Create the map.
    AssignmentPool assignment_pool;
    // 3.1 Convert node item into assignment.
    for(auto iter = node_map.begin(); iter != node_map.end(); iter++)
    {
        AssignmentPtr pk_ptr = boost::make_shared<Assignment>(iter->first, iter->second);
        assignment_pool.addAssignment(pk_ptr);
    }

    // 3.2 Deal with acetate group.
    string mod_symbol = "Ac";
    set<BackbonePtr> bone_set = assignment_pool.selectQualifiedAssignments(mod_symbol);

    for(auto iter = bone_set.begin(); iter != bone_set.end(); iter++)
    {
        set<int> mod_num_set = (*iter)->getModNumbers();
        for(auto num_it = mod_num_set.begin(); num_it != mod_num_set.end(); num_it++)
        {
            set<AssignmentPtr> assign_set = (*iter)->getAssignmentsByModNumber(*num_it);
            cout << "Assignment size:" << assign_set.size();
        }

    }

    FullMap ass_graph(gs, mod_symbol, bone_set);
    cout << ass_graph << "\n";

    // 3.3 Deal with sulfate group
    mod_symbol = "SO3";
    set<BackbonePtr> bone_set2 = assignment_pool.selectQualifiedAssignments(mod_symbol);

    FullMap sulfate_graph(gs, mod_symbol, bone_set2);

    cout << sulfate_graph << "\n";
  }
  catch (std::exception& e)
  {
    cout << "Exception: " << e.what() << endl;
  }
}