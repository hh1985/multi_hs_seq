#include "GAGPL/SPECTRUM/IsotopicDistribution.h"
#include <time.h>
#include <boost/filesystem.hpp>

#pragma warning(disable : 4996)

using namespace std;
using namespace gag;
using namespace boost::filesystem;

int main(int argc, char* argv[])
{
	try {

		vector<pair<Composition, int> > compo_vec;
		//Composition compo1("C50H71N12O13"); 
    Composition compo1("C44H64N3O12S4");
    Composition compo2("C23832H37816N6528O7031S170");
		compo_vec.push_back(std::make_pair(compo1, 11));
		compo_vec.push_back(std::make_pair(compo2, 1325));

		clock_t t;

		// Parse the input path and generate output file name.
		for(vector<pair<Composition, int> >::iterator iter = compo_vec.begin(); iter != compo_vec.end(); iter++)
		{
			// Initial time.
			t = clock();

			IsotopicDistribution iso_dist(iter->first, iter->second);
			AggregatedIsotopicVariants peakset = iso_dist.getAggregatedIsotopicVariants();

			// Ending time.
			t = clock()-t;
			cout << "It took " << ((double)t)/CLOCKS_PER_SEC << " seconds to calculate " << iter->first.getCompositionString() << endl;

			cout << "Composition: " << iter->first.getCompositionString() << endl;
			cout << "Monoisotopic mass: " << iter->first.getMass() << endl;
			cout << "Average mass: " << iter->first.getAverageMass() << endl;
			cout << "Estimated average mass: " << iso_dist.getAverageMass() << endl;
			peakset.printPeakList<peak_mz>(cout);

		}
	}
	catch(std::exception& e)
	{
		cout << "Exception: " << e.what() << endl;
	}

	cout << "The End!" << endl;
	//cin.get();

	return 0;


}