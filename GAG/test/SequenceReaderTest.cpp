#include <iostream>
#include <GAGPL/IO/SequenceReader.h>
#include <GAGPL/SPECTRUM/SimpleFinder.h>


using namespace std;
using namespace gag::io;
using namespace gag;

// argv[1] -- filename
int main(int argc, char *argv[])
{
	
	SequenceReader seq_reader;
	
	cout << "*************************************" << endl;
	cout << "*   Test 1: Code-based Sequence     *" << endl;
	cout << "*************************************" << endl;

	// 1. Specify the folder for storing the structure .gl file.
	std::string filename = "../data/structure/";
	filename.append("dp15.gl");

	// 2. Read the sequence.
	GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
	seq_reader.readSequenceStructure(filename, gs);

	// 3. Explore the sequence information.
	gs->printStructure();

	cout << "*************************************" << endl;
	cout << "*   Test 2: Spectrum                *" << endl;
	cout << "*************************************" << endl;

	std::string spec_file = "../data/spectrum/";
	spec_file.append("hex6_5_NETD.txt");

	RichList spectrum;
	seq_reader.readSpectrum(spec_file, spectrum);
	cout << "Spectrum size: " << spectrum.getSize() << endl;

	cout << "*************************************" << endl;
	cout << "*   Test 3: Monoisotopic peaks      *" << endl;
	cout << "*************************************" << endl;

	std::string mono_file = "../data/output/";
	mono_file.append("hex6_5_NETD_mono_list.txt");

	std::set<MonoPeakPtr> mono_set;
	seq_reader.readMonoPeakList(mono_file, mono_set);
	cout << "Mono peak list size: " << mono_set.size() << endl;

	cin.get();
	
}