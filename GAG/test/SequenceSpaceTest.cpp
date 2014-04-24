#include <GAGPL/GAGLIBRARY/SequenceSpace.h>
#include <GAGPL/GAGLIBRARY/SequencePrediction.h>
#include <GAGPL/GAGLIBRARY/NaiveMethods.h>
#include <GAGPL/IO/SequenceReader.h>

using namespace std;
using namespace gag;
using namespace gag::io;

int main(int argc, char **argv)
{
  try
  {
    cout << "*************************************" << endl;
    cout << "*   Test 0: Enumerate sequence      *" << endl;
    cout << "*************************************" << endl;

    SequenceReader seq_reader;

    // 1. Load the sequence from file.
    std::string structure_filename = "./data/structure/";
    structure_filename.append(argv[1]);

    GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
    // The sequence used for comparison.
    GlycanSequencePtr compare_gs = boost::make_shared<GlycanSequence>();
    seq_reader.readSequenceStructure(structure_filename, gs);
    seq_reader.readSequenceStructure(structure_filename, compare_gs);

    SequenceSpace seq_space(gs);
    SequenceSpace compare_space(compare_gs);

#ifdef _DEBUG
    //SequenceList seq_list = seq_space.generateModificationCombination();
    //SequenceList seq_list;
    set<string> merged_string;

    int seq_count(0);

    while(1)
    {
      ModificationSequence mod_seq = seq_space.nextModificationSequence();
      
      if(mod_seq.size() == 0) break;
     
      seq_count++;

      //seq_list.push_back(mod_seq);
      std::string m_str = seq_space.getMergedCodeString(mod_seq);
      cout << seq_space.getSequenceString(mod_seq) << " " << m_str << "\n";
      merged_string.insert(m_str);

    }

    cout << "Filename: " << structure_filename << "\n";

    cout << "There are " << seq_count << " sequences.\n";

    cout << "There are " << merged_string.size() << " merged sequences." << endl;

#endif // _DEBUG

    //cout << "*************************************" << endl;
    //cout << "*   Test 1: Calculate score         *" << endl;
    //cout << "*************************************" << endl;
    ////Construct the library tree, which is used as the complete set.
    //FragmentTree glyco_tree(gs);
    //FragmentTree compare_tree(compare_gs);

    //std::string mono_filename = "./data/output/";
    //mono_filename.append(argv[2]);
    //mono_filename.append("_mono_list.txt");

    //std::set<MonoPeakPtr> mono_set;
    //seq_reader.readMonoPeakList(mono_filename, mono_set);

    //std::multimap<MonoPeakPtr, NodeItem> node_map = glyco_tree.searchLibrary(mono_set);
    //std::multimap<MonoPeakPtr, NodeItem> compare_map = compare_tree.searchLibrary(mono_set);

    //SequencePrediction seq_predictor(node_map, gs);
    //NaiveMethods compare_predictor(compare_map, compare_gs);

    //std::set<std::string> mod_types = gs->getModificationTypes();
    //std::cout << "Position\tModification\tIntensity\n";
    //for(auto iter = mod_types.begin(); iter != mod_types.end(); iter++)
    //{
    //  ModificationDistribution mod_dist = seq_predictor.getModificationDistribution(*iter);
    //  printModificationDistribution(mod_dist, *iter);
    //  
    //}

    //std::string results_file = "./data/output/";
    //results_file.append(argv[2]);
    //results_file.append("_score.txt");

    //std::ofstream outfile(results_file.c_str());
    //if(!outfile.is_open()) throw std::runtime_error("Unable to open the file!");
    //
    //// Print head.
    //outfile << "SEQ\tM_SEQ\tCoverage\tGP\tCost\tM_Coverage\tM_GP\tM_Cost\n";

    //while(1) {
    //  ModificationSequence mod_seq = seq_space.nextModificationSequence();
    //  if(mod_seq.size() == 0) break;
    //  
    //  outfile << seq_space.getSequenceString(mod_seq) << "\t" << seq_space.getMergedCodeString(mod_seq);

    //  for(int i = 0; i < 2; i++) {
    //    bool flag = (i == 0 ? false : true);

    //    //NaiveMethods compare_predictor(compare_gs, mod_seq);
    //    compare_predictor.setStatus(flag);

    //    outfile << "\t" << compare_predictor.calculateCoverage(mod_seq);
    //    outfile << "\t" << compare_predictor.calculateGoldenPairNum(mod_seq);

    //    if(!flag)
    //      outfile << "\t" << seq_predictor.calculateCost(mod_seq);
    //    else
    //      outfile << "\t" << seq_predictor.calculateMergedCost(mod_seq) << "\n";          
    //  }
    //  //seq_reader.writePredictionResults(results_file, seq_scores);
    //}

    //outfile.close();


  }

  catch (std::exception& e)
  {
    cout << "Exception: " << e.what() << endl;
  }
	
}