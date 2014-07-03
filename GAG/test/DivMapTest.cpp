#include <GAGPL/SEQUENCING/DivMap.h>
#include <GAGPL/IO/SequenceReader.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/SEQUENCING/AssignmentPool.h>

using namespace std;
using namespace gag;
using namespace gag::io;

int main(int argc, char *argv[])
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

    vector<string> symbol_vec;
    symbol_vec.push_back("Ac");
    symbol_vec.push_back("SO3");

    // 3.2 Predict the Ac position.
   
    for(size_t i = 0; i < symbol_vec.size(); i++) 
    {
      string mod_symbol = symbol_vec[i];

      ModificationSites full_sites = gs->getModificationSitesBySymbol(mod_symbol, 1);
      int full_num = gs->getModificationConstraint(mod_symbol);

      cout << "Mod type: " << mod_symbol << "\n";
      cout << "Full sites: " << full_sites << "\n";
      cout << "Mod number: " << full_num << "\n";

      set<DivisionPtr> div_set = assignment_pool.selectAssignments(mod_symbol, full_sites, full_num);
      cout << "Selected division size: " << div_set.size() << "\n";

      DivMap div_map(gs, mod_symbol);
      for(auto iter = div_set.begin(); iter != div_set.end(); iter++) {
        cout << "\nAdd division node: " << **iter << "\n";
        div_map.addDivisionNode(*iter);

        cout << "\nUpdated map:\n";
        cout << div_map << "\n";
      }

      cout << "End of " << mod_symbol << "\n";


      // TBD: Filter assignments once the positions has been identified. 
      // 1. Convert path to modification distribution.
      // 2. Update assignment_pool based on identified modification list.
    }

  } catch(std::exception& e) {
    cout << "Exception: " << e.what() << endl;
  }
}