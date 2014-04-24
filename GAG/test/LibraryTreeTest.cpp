#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/FRAGMENTATION/FragmentationTable.h>
#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <GAGPL/CHEMISTRY/Modifier.h>
#include <GAGPL/MATH/MassConversion.h>
#include <GAGPL/IO/SequenceReader.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <utility>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <time.h>


using namespace std;
using namespace gag;
using namespace msmath;
using namespace gag::io;

int main(int argc, char* argv[])
{
	//// match peak list against tree_vec. 

	std::cout << "**************************************" << std::endl;
	std::cout << "         Test 1: Node Item            " << std::endl;
	std::cout << "**************************************" << std::endl;

	GlycanComposition g_compo;
	g_compo.addGlycanComposition("DeltaGlcA", 1);
	g_compo.addGlycanComposition("GlcA", 2);
	g_compo.addGlycanComposition("GlcN", 3);
	g_compo.addModificationComposition("Ac", 1);
	g_compo.addModificationComposition("SO3", 5);

	GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
	gs->buildByGAGComposition(g_compo);

	// Create the fragment.
	FragmentPtr fg_ptr(new Fragment(gs));
	// Set the fragment type: Y5 A5(1,5)
	FragmentPosition fp1(0,0,0,0);
	//FragmentPosition fp2(0,4,1,5);
	fg_ptr->setFragmentation("Y", fp1);
	//fg->setFragmentation("A", fp2);
	fg_ptr->printFragment();

	//// Copy the fragment
	FragmentPtr fgx(new Fragment(*fg_ptr));

	FragmentPosition fp2(0,4,1,5);
	fgx->setFragmentation("A", fp2);
	fgx->printFragment();

	std::cout << "Original fragment" << std::endl;
	fg_ptr->printFragment();

	NodeItem node(fg_ptr);
	node.appendMassShift("H2O", 1);
	node.appendMassShift("H", -1);
	//node.printNodeItem();

	std::cout << "**************************************" << std::endl;
	std::cout << "         Test 2: Node Item            " << std::endl;
	std::cout << "**************************************" << std::endl;
	
	std::string filename = "./data/output/";
	filename.append(argv[2]);

	std::string mono_file = filename;
	mono_file.append("_mono_list.txt");

	// Load the data and search for matching.
	std::set<MonoPeakPtr> peak_list;
	loadMonoPeakList(mono_file, peak_list);

	SequenceReader seq_reader;

	// 1. Load the sequence from file.
	std::string structure_filename = "./data/structure/";
	structure_filename.append(argv[1]);

	GlycanSequencePtr glyco_seq = boost::make_shared<GlycanSequence>();
	seq_reader.readSequenceStructure(structure_filename, glyco_seq);


	// Start the clock.
	clock_t t = clock();
	FragmentTree glyco_tree(glyco_seq);

	// Update the time value.
	t = clock() - t;
	std::cout << "It took " << ((double)t)/CLOCKS_PER_SEC << " seconds to construct the library." << std::endl;
	std::string library_name = filename;
	library_name.append("_library.txt");
	glyco_tree.exportLibrary(library_name);

	std::cout << "**************************************" << std::endl;
	std::cout << "       Test 3: Library matching       " << std::endl;
	std::cout << "**************************************" << std::endl;

	std::multimap<MonoPeakPtr, NodeItem> node_map = glyco_tree.searchLibrary(peak_list);
	std::cout << "There are " << node_map.size() << " matches!" << std::endl; 

	std::string matched_pk_file = filename;
	matched_pk_file.append("_assignments.txt");
	std::ofstream outfile(matched_pk_file.c_str());
	outfile.precision(5);

	if(outfile.is_open()) {
		// print title.
		outfile << "MZ\tCharge\tEXP_Mass\tTheo_Mass\tError\tCleavage_Type\tCleavage_Num\tComposition_Shift\tComposition\tAc\tSO3\n";
		for(std::multimap<MonoPeakPtr, NodeItem>::iterator iter = node_map.begin(); iter != node_map.end(); iter++)
		{
			double mass = msmath::calculateMass(iter->first->mz, -1 * iter->first->z);
      double ppm = (iter->second.getMass() - mass)/mass * 1e6;

			outfile << std::fixed << iter->first->mz << "\t" << iter->first->z << "\t" << mass << "\t" << iter->second.getMass() << "\t" << ppm << "\t" << iter->second.getCleavageType() << "\t" << iter->second.getCleavageNum() << "\t" << iter->second.getCompositionShift() << "\t" << iter->second.getCompositionString() << "\t" << iter->second.getModificationNum("Ac") << "\t" << iter->second.getModificationNum("SO3")  << "\n";
		}
		outfile.close();
	} else {
		std::cout << "Unable to open file!\n" << std::endl;
	}

}