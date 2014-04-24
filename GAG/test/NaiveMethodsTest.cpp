#include <GAGPL/GAGLIBRARY/NaiveMethods.h>
#include <GAGPL/IO/SequenceReader.h>

using namespace std;
using namespace gag;
using namespace gag::io;

int main(int argc, char **argv)
{
  try
  {
    SequenceReader seq_reader;

    // The test file is hex6.gl
    std::string name = "hex6.gl";
    std::string pk_list = "hex6_5_NETD";

    // 1. Load the sequence from file.
    std::string structure_filename = "./data/structure/";
    structure_filename.append(name);

    GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
    seq_reader.readSequenceStructure(structure_filename, gs);

    // 2. Load the mono_isotopic peak list from file. The precursor information should be with the spectrum file.
    std::string mono_filename = "./data/output/";
    mono_filename.append(pk_list);
    mono_filename.append("_mono_list.txt");


    std::set<MonoPeakPtr> mono_set;
    seq_reader.readMonoPeakList(mono_filename, mono_set);

    FragmentTree glyco_tree(gs);
    std::cout << "Get data size: " << glyco_tree.getNodeSize() << "\n";
    std::multimap<MonoPeakPtr, NodeItem> node_map = glyco_tree.searchLibrary(mono_set);

    std::cout << "Assignment data size: " << node_map.size() << "\n";

    // Construct a mod sequence.
    ModificationSequence mod_seq;
    ModificationSites ac_sites;
    //ModificationPosition mod_pos1(0, 3);
    ac_sites.insert(ModificationPosition(0, 3, 2));
    mod_seq.insert(std::make_pair("Ac", ac_sites));

    ModificationSites sulfate_sites;
    sulfate_sites.insert(ModificationPosition(0, 0, 2));
    sulfate_sites.insert(ModificationPosition(0, 1, 2));
    sulfate_sites.insert(ModificationPosition(0, 1, 6));
    sulfate_sites.insert(ModificationPosition(0, 3, 6));
    sulfate_sites.insert(ModificationPosition(0, 5, 2));
    sulfate_sites.insert(ModificationPosition(0, 5, 6));
    mod_seq.insert(std::make_pair("SO3", sulfate_sites));

    ModificationSites alternative_sites;
    ModificationSequence alternative_seq;
    alternative_seq.insert(std::make_pair("Ac", ac_sites));
    alternative_sites.insert(ModificationPosition(0, 0, 2));
    alternative_sites.insert(ModificationPosition(0, 1, 6));
    alternative_sites.insert(ModificationPosition(0, 1, 3));
    alternative_sites.insert(ModificationPosition(0, 3, 6));
    alternative_sites.insert(ModificationPosition(0, 5, 6));
    alternative_sites.insert(ModificationPosition(0, 5, 3));
    alternative_seq.insert(std::make_pair("SO3", sulfate_sites));


    cout << "*************************************" << endl;
    cout << "*   Test 0: Filter data             *" << endl;
    cout << "*************************************" << endl;
    
    NaiveMethods nm(node_map, gs);
    /*std::multimap<MonoPeakPtr, NodeItem> new_data = nm.filterData(node_map);
    cout << "Old data: " << node_map.size() << "\tNew data: " << new_data.size() << "\n";*/
    
    cout << "*************************************" << endl;
    cout << "*   Test 1: Coverage methods        *" << endl;
    cout << "*************************************" << endl;

    cout << "Merged version:\n";
    
    cout << nm.calculateCoverage(mod_seq) << "\n";
    cout << nm.calculateCoverage(alternative_seq) << "\n";

    cout << "Full version:\n";
    nm.setStatus(0);
    cout << nm.calculateCoverage(mod_seq) << "\n";
    cout << nm.calculateCoverage(alternative_seq) << "\n";

    cout << "*************************************" << endl;
    cout << "*   Test 2: Golden pair             *" << endl;
    cout << "*************************************" << endl;

    cout << "Merged version:\n";
    nm.setStatus(1);
    cout << nm.calculateGoldenPairNum(mod_seq) << "\n";
    cout << nm.calculateGoldenPairNum(alternative_seq) << "\n";

    cout << "Full version:\n";
    nm.setStatus(0);
    cout << nm.calculateGoldenPairNum(mod_seq) << "\n";
    cout << nm.calculateGoldenPairNum(alternative_seq) << "\n";
  }
  catch (std::exception& e)
  {
    cout << "Exception: " << e.what() << endl;
  }
}
