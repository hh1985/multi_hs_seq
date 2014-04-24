/********************************************************************
	created:	2012/11/16
	created:	16:11:2012   16:59
	filename: 	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\test\EnvelopFinderTest.cpp
	file path:	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\test
	file base:	EnvelopFinderTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include <GAGPL/SPECTRUM/EnvelopFinder.h>
#include <GAGPL/SPECTRUM/PeakList.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <boost/algorithm/string.hpp>

using namespace gag;
using namespace std;

RichList readSpectrum()
{
	std::cout << "Load spectrum..." << std::endl;

	std::ifstream infile("../data/spectrum.txt");

	std::string line;

	RichList spec;
	if(infile.is_open())
	{
		while(std::getline(infile, line))
		{
			std::istringstream is;
			is.str(line);
			double k1, k2, k3, k4; 
			is >> k1 >> k2 >> k3 >> k4;
			//gag::Peak pk(k1, k2);
			RichPeakPtr pk = createRichPeak(k1, k2, k3, k4);
			spec.addPeak(pk);
		}
	}
	infile.close();
	infile.clear();
	return spec;

}

int main(int argc, char *argv[])
{
	// Test 1
	cout << "Test 1: Manual interpretation" << endl;
	RichList spec1;
	RichPeakPtr pk0 = createRichPeak(376.49227,2125405,5184,3.7);
	RichPeakPtr pk1 = createRichPeak(376.99506,51087228,311536,130);
	RichPeakPtr pk2 = createRichPeak(377.32956,12896527,68383,31.5);
	RichPeakPtr pk3 = createRichPeak(377.65959,3833504,4942,8.1);
	RichPeakPtr pk4 = createRichPeak(377.66405,4087202,5864,8.7);
	RichPeakPtr pk5 = createRichPeak(377.99828,1626549,5004,2.4);
	spec1.addPeak(pk0);
	spec1.addPeak(pk1);
	spec1.addPeak(pk2);
	spec1.addPeak(pk3);
	spec1.addPeak(pk4);
	spec1.addPeak(pk5);

	EnvelopFinder env_finder1(spec1);

	std::vector<EnvelopPtr> env_pool1 = env_finder1.getEnvelops();
	BOOST_FOREACH(EnvelopPtr env, env_pool1)
	{
		env_finder1.printEnvelop(env);
	}
	// Test 2
	//cout << endl;
	//cout << "Test 2: Real spectrum" << endl;
	//// 1. TBD: Read the data and examine the format of the data.
	//
	//RichList spectrum = readSpectrum();
	//spectrum.printPeakList<peak_mz>();

	//EnvelopFinder env_finder2(spectrum);

	//// 2. print the identified envelops.
	//std::vector<EnvelopPtr> env_pool = env_finder2.getEnvelops();

	//BOOST_FOREACH(EnvelopPtr env, env_pool)
	//{
	//	env_finder2.printEnvelop(env);
	//}

}