/********************************************************************
	created:	2014/01/15
	created:	15:1:2014   13:24
	filename: 	HHSequencing.cpp
	file path:	GAGPL\GAGLIBRARY
	file base:	HHSequencing
	file ext:	cpp
	author:		Han Hu
	
	purpose:	The main function for HS-SEQ algorithm
*********************************************************************/
//#include <iostream>
#include <GAGPL/GAGLIBRARY/SequencePrediction.h>
#include <GAGPL/IO/SequenceReader.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <time.h>

using namespace std;
using namespace gag;
using namespace gag::io;

//void printModificationDistribution(const ModificationDistribution& dist, const std::string& mod_name)
//{
//  // Print title.
//  cout << "Position\tModification\tIntensity" << endl;
//
//  for(ModificationDistribution::const_iterator iter = dist.begin(); 
//    iter != dist.end(); iter++)
//  {
//    cout << iter->first.printString() << "\t" << mod_name << "\t" << iter->second << endl;
//  }
//}
// argv[1] -- structure file.
// argv[2] -- spectrum file
int main(int argc, char **argv)
{
  try
  {
    clock_t t;
    t = clock();
    SequenceReader seq_reader;

    // 1. Load the sequence from file.
    std::string structure_filename = "./data/structure/";
    structure_filename.append(argv[1]);

    GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
    seq_reader.readSequenceStructure(structure_filename, gs);

    //SequenceReader seq_reader;

    //// 1. Load the sequence from file.
    //std::string structure_filename = "../data/structure/";
    //structure_filename.append(argv[1]);

    //GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
    //seq_reader.readSequenceStructure(structure_filename, gs);

    // 2. Load the mono_isotopic peak list from file. The precursor information should be with the spectrum file.
    std::string mono_filename = "./data/output/";
    mono_filename.append(argv[2]);
    mono_filename.append("_mono_list.txt");

    std::set<MonoPeakPtr> mono_set;
    seq_reader.readMonoPeakList(mono_filename, mono_set);

    FragmentTree glyco_tree(gs);
    std::cout << "Get data size: " << glyco_tree.getNodeSize() << "\n";
    std::multimap<MonoPeakPtr, NodeItem> node_map = glyco_tree.searchLibrary(mono_set);

    std::cout << "Assignment data size: " << node_map.size() << "\n";

    // Predict the structure from the backbone and mass list.
    SequencePrediction seq_predictor(node_map, gs);

    std::set<std::string> mod_types = gs->getModificationTypes();

    std::string results_file = "./data/output/";
    results_file.append(argv[2]);
    results_file.append("_dist.txt");
    std::ofstream outfile(results_file.c_str());

    if(outfile.is_open()) throw std::runtime_error("Unable to open the file!");

    outfile << "Position\tModification\tIntensity\n";
    for(auto iter = mod_types.begin(); iter != mod_types.end(); iter++)
    {
      ModificationDistribution mod_dist = seq_predictor.getModificationDistribution(*iter);
      for(ModificationDistribution::const_iterator mod_iter = mod_dist.begin(); 
        mod_iter != mod_dist.end(); mod_iter++)
      {
        outfile << mod_iter->first.printString() << "\t" << *iter << "\t" << mod_iter->second << "\n";
      }
      
    }
    outfile.close();

    ////int ac_num = gs->getModificationConstraint("Ac");
    //ModificationDistribution ac_dist = seq_predictor.getModificationDistribution("Ac");
    //printModificationDistribution(ac_dist, "Ac");

    //cout << "*************************************" << endl;
    //cout << "*   Test 2: Distribution of SO3     *" << endl;
    //cout << "*************************************" << endl;

    //ModificationDistribution s_dist = seq_predictor.getModificationDistribution("SO3");
    //cout << "Final distribution of SO3: " << std::endl;
    ////for(ModificationDistribution::iterator iter = s_dist.begin(); 
    ////	iter != s_dist.end(); iter++)
    ////{
    ////	std::cout << "Modification Position: " << iter->first.getBranchID() << " " << iter->first.getMonosaccharideID() << " " << iter->first.site_id << std::endl;
    ////	std::cout << "mod intensity: " << iter->second << std::endl;
    ////}

    //printModificationDistribution(s_dist, "SO3");

    t = clock() - t;
    //std::cout << "Done!" << std::endl;
    cout << "Time cost: " << (double)t/CLOCKS_PER_SEC << "s\n";
    //cin.get();
  }
  catch(std::exception& e)
  {
    cout << "Exception: " << e.what() << endl;
  }
}