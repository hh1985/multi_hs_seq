#include <GAGPL/SPECTRUM/PseudoclusterFinder.h>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>


int main(int argc, char *argv[]) 
{
	// 1. Read test data.
	std::ifstream infile (argv[1]);

	std::vector<gag::Peak> peaks;

	//vector< iterator_range<string::iterator> > vec;

	std::string line;

	if(infile.is_open())
	{
		while(std::getline(infile, line))
		{
			//std::vector<std::string> peak;
			//boost::split(peak, line, boost::is_any_of("\t"));	
			std::istringstream is;
			is.str(line);
			double k1; float k2;
			is >> k1 >> k2;
			gag::Peak pk(k1, k2);
			peaks.push_back(pk);
		}
	}
	infile.close();

	gag::PseudoclusterFinder pf;
	gag::Param& param = gag::Param::instance();
	param.setParameter("max_charge",4);
	param.setParameter("min_charge",1);
	param.setParameter<double>("tolerance",5e-6);
	param.setParameter<double>("A2tolerance", 1e-5);
	param.setParameter("min_num_peaks",2);
	param.setParameter("max_continuous_peaks",2);
	param.setParameter("gap_distance",1.1);

	std::map<std::string, gag::PseudoclusterList> cluster;
	pf.run(peaks,cluster);
	std::cout << "Cluster_ID\tPeak_ID\tCharge\tMax\tMonoMass" << std::endl; 
	std::map<std::string, gag::PseudoclusterList>::iterator it;
	for(it=cluster.begin(); it !=cluster.end(); it++)
	{
	//	//std::cout << "Charge: " << (*it).first << std::endl;
		gag::PseudoclusterList& pl = (*it).second;
		for(std::vector<gag::Pseudocluster>::iterator ip = pl.begin(); ip != pl.end(); ip++)
		{
			std::cout << (*ip).getClusterID() << "\t" << (*ip).getHeadID() << "\t" << (*ip).getChargeState() 
				<< "\t" << (*ip).getMaxPeaks() << "\t" << (*ip).getPeakList().at(0).getMZ() << std::endl;

			//gag::IsopeakList& xy = (*ip).getPeakList();
	//		for(size_t xx=0; xx < xy.size(); xx++)
	//		{
	//			std::cout << xy[xx].getMZ() << " ";
	//		}
	//		std::cout << std::endl;
	//		//std::cout << (*ip).getChargeState() << "\t" << xy[0].getMZ() << std::endl;
		}
	}

	return 0;
}
