/********************************************************************
	created:	2014/01/15
	created:	15:1:2014   13:26
	filename: 	SequencePredictionTest.cpp
	file path:	test
	file base:	SequencePredictionTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	Test file for HS-SEQ algorithm
*********************************************************************/
#include <GAGPL/GAGLIBRARY/SequencePrediction.h>
#include <GAGPL/IO/SequenceReader.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <time.h>

using namespace std;
using namespace gag;
using namespace gag::io;

//void printModificationDistribution(const ModificationDistribution& dist, const std::string& mod_name)
//{
//	// Print title.
//	cout << "Position\tModification\tIntensity" << endl;
//
//	for(ModificationDistribution::const_iterator iter = dist.begin(); 
//		iter != dist.end(); iter++)
//	{
//		cout << iter->first.printString() << "\t" << mod_name << "\t" << iter->second << endl;
//	}
//}
// argv[1] -- structure file.
// argv[2] -- spectrum file
int main(int argc, char **argv)
{
  clock_t t;
  t = clock();
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

	cout << "*************************************" << endl;
	cout << "*   Test 1: Distribution of Ac      *" << endl;
	cout << "*************************************" << endl;

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

	// Predict the structure from the backbone and mass list.
	SequencePrediction seq_predictor(node_map, gs);

	//int ac_num = gs->getModificationConstraint("Ac");
  ModificationDistribution ac_dist = seq_predictor.getModificationDistribution("Ac");
  printModificationDistribution(ac_dist, "Ac");

	cout << "*************************************" << endl;
	cout << "*   Test 2: Distribution of SO3     *" << endl;
	cout << "*************************************" << endl;
	
	ModificationDistribution s_dist = seq_predictor.getModificationDistribution("SO3");
	cout << "Final distribution of SO3: " << std::endl;
  printModificationDistribution(s_dist, "SO3");
  
  cout << "****************************************" << endl;
  cout << "*   Test 3: Test of sequence score     *" << endl;
  cout << "****************************************" << endl;

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

  std::cout << "Cost: " << seq_predictor.calculateCost(mod_seq) << "\n";
  std::cout << "Merged cost: " << seq_predictor.calculateMergedCost(mod_seq) << "\n";
  //std::cout << "Coverage: " << seq_predictor.calculateCoverage(node_map) << "\n";
  //std::cout << "Golden pair: " << seq_predictor.calculateGoldenPairNum(node_map) << "\n";

  cout << "****************************************" << endl;
  cout << "*   Test 4: Compare of isomers         *" << endl;
  cout << "****************************************" << endl;
  ModificationSequence mod_seq1;
  ModificationSites new_ac_sites;

  new_ac_sites.insert(ModificationPosition(0,5,2));
  mod_seq1.insert(std::make_pair("Ac", new_ac_sites));

  ModificationSites sulfate_sites1;
  sulfate_sites1.insert(ModificationPosition(0, 0, 2));
  sulfate_sites1.insert(ModificationPosition(0, 1, 2));
  sulfate_sites1.insert(ModificationPosition(0, 1, 3));
  sulfate_sites1.insert(ModificationPosition(0, 1, 6));
  sulfate_sites1.insert(ModificationPosition(0, 3, 2));
  sulfate_sites1.insert(ModificationPosition(0, 3, 3));

  mod_seq1.insert(std::make_pair("SO3", sulfate_sites1));

  ModificationSequence mod_seq2;
  mod_seq2.insert(std::make_pair("Ac", new_ac_sites));

  ModificationSites sulfate_sites2;
  sulfate_sites2.insert(ModificationPosition(0, 0, 2));
  sulfate_sites2.insert(ModificationPosition(0, 1, 2));
  sulfate_sites2.insert(ModificationPosition(0, 1, 3));
  sulfate_sites2.insert(ModificationPosition(0, 1, 6));
  sulfate_sites2.insert(ModificationPosition(0, 3, 3));
  sulfate_sites2.insert(ModificationPosition(0, 3, 6));

  mod_seq2.insert(std::make_pair("SO3", sulfate_sites2));


  std::cout << "SEQ1 Cost: " << seq_predictor.calculateCost(mod_seq1) << " ";
  std::cout << "Merged cost: " << seq_predictor.calculateMergedCost(mod_seq1) << "\n";

  std::cout << "SEQ2 Cost: " << seq_predictor.calculateCost(mod_seq2) << " ";
  std::cout << "Merged cost: " << seq_predictor.calculateMergedCost(mod_seq2) << "\n";





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