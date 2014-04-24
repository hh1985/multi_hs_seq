/********************************************************************
	created:	2012/11/16
	created:	16:11:2012   16:59
	filename: 	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\test\EnvelopFinderTest.cpp
	file path:	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\test
	file base:	SimpleFinderTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/SPECTRUM/SimpleFinder.h"
#include "GAGPL/SPECTRUM/PeakList.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <boost/algorithm/string.hpp>
#include <stdlib.h>

using namespace gag;
using namespace std;
using namespace param;

RichList readSpectrum(const std::string& filename)
{
	std::cout << "Load spectrum..." << std::endl;

	std::ifstream infile(filename.c_str());

	std::string line;

	RichList spec;
	if(infile.is_open())
	{
		while(std::getline(infile, line))
		{
			std::istringstream is;
			is.str(line);
			double k1, k2, k3, k4; 
			// The four columns are: m/z, intensity, resolution, s/n
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

// argv[1] -- filename without suffix.
// argv[2] -- precursor m/z
// argv[3] -- precursor charge

int main(int argc, char *argv[])
{
	// Test 1
	//cout << "Test 1: Manual interpretation" << endl;

	//EnvelopFinder env_finder1(spec1);

	//std::vector<EnvelopPtr> env_pool1 = env_finder1.getEnvelops();
	//BOOST_FOREACH(EnvelopPtr env, env_pool1)
	//{
	//	env_finder1.printEnvelop(env);
	//}
	// Test 2
	//cout << endl;
	try {
		cout << "Test 2: Real spectrum" << endl;
		// 1. TBD: Read the data and examine the format of the data.

		std::string filename = "./data/spectrum/";
		filename.append(argv[1]);
		filename.append(".txt");

		RichList spectrum = readSpectrum(filename);

		Param& param = Param::Instance();
		param.setParameter<double>("precursor_mz", atof(argv[2]));
		param.setParameter<int>("precursor_charge", atoi(argv[3]));

		// Threshold for most abundant peak in an envelop (mono peak in current case.
		//param.setParameter<double>("signal_noise", 15.0);
		param.setParameter<std::string>("instrument_type", "ft");

		// threshold for observed peak.
		param.setParameter<double>("lower_bound_sn",5.0);
		param.setParameter<double>("learning_rate", 0.01);
		param.setParameter<double>("max_intensity_shift", 0.3);
		param.setParameter<int>("max_sulfur_num",6);
		param.setParameter<int>("num_missing_peaks",1);

		SimpleFinder env_finder2(spectrum);

		// 2. print the identified envelops.
		std::multimap<RichPeakPtr, EnvelopPtr> env_pool = env_finder2.getEnvelops();

		cout << "Output the results into file." << endl;

		std::string path = "./data/output/";
		path.append(argv[1]);
		path.append("_mono_list.txt");

		ofstream outfile(path.c_str());

		if(outfile.is_open()) {
			outfile << "MZ\tIntensity\tCharge\n"; 
			outfile.precision(5);
			cout.precision(5);

			for(std::multimap<RichPeakPtr, EnvelopPtr>::iterator iter = env_pool.begin(); iter != env_pool.end(); iter++)
			{
				env_finder2.printEnvelop(iter->second);

				outfile << fixed << iter->first->mz << "\t" << iter->first->intensity << "\t" << iter->second->charge_state << "\n"; 		
			}

			outfile.close();
		} else {
			std::cout << "Unable to open file!\n" << std::endl;
		}

		std::string und_pks = "./data/output/";
		und_pks.append(argv[1]);
		und_pks.append("_undetermined_list.txt");
		ofstream outfile2(und_pks.c_str());

		if(outfile2.is_open()) {
			outfile2 << "MZ\tIntensity\tCharge\n";
			outfile2.precision(5);

			std::multimap<RichPeakPtr, int>& pks = env_finder2.getUndeterminedPeaks();
			for(std::multimap<RichPeakPtr, int>::iterator iter = pks.begin(); iter!= pks.end(); iter++) 
			{
				if(iter->first->pk_type == "NOISE")
					continue;
				outfile2 << fixed << iter->first->mz << "\t" << iter->first->intensity << "\t" << iter->second << "\n";
			}
			outfile2.close();

		} else {
			std::cout << "Unable to open file " << und_pks << std::endl;
		}
	} catch(std::exception& e) {
		cout << "Exception: " << e.what() << endl;
	}

	return 0;

}