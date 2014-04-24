/********************************************************************
	created:	2012/11/13
	created:	13:11:2012   9:12
	filename: 	PeakListTest.cpp
	file path:	GAG\test
	file base:	PeakListTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include <GAGPL/SPECTRUM/PeakList.h>
#include <boost/make_shared.hpp>

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace gag;

	RichList peaklist;

	// 1. Add peaks into PeakList
	RichPeakPtr pk1 = createRichPeak(475.2, 28.0, 30.0, 4.2);
	peaklist.addPeak(pk1);

	RichPeakPtr pk2 = createRichPeak(476.2, 80.0, 123.0, 8.7);
	peaklist.addPeak(pk2);

	RichPeakPtr pk3 = createRichPeak(477.2, 10.0, 15.0, 3.2);
	pk3->pk_type = "ISO";
	peaklist.addPeak(pk3);
	
	RichPeakPtr pk4 = createRichPeak(478.2, 5.0, 5.0, 2.2);
	pk4->pk_type = "ISO";
	peaklist.addPeak(pk4);

	// print out the size of the peaklist.
	cout << "Size: " << peaklist.getSize() << endl;

	RichPeakListByResolution& i = peaklist.getPeakListByType<peak_resolution>();

	peaklist.print<peak_resolution>(i);

	cout << "Get peaks list by type." << endl;
	pk4->pk_status = false;
	RichPeakListByID& pks_id = peaklist.getPeakListByType<peak_id>();
	RichPeakListByID::iterator it = pks_id.find(3);
	RichPeakListByType::iterator type_it = peaklist.getPeakContainer().project<peak_type>(it);

	RichPeakListByType& pks_type = peaklist.getPeakListByType<peak_type>();
	peaklist.modifyType(type_it, "ENV");
	peaklist.modifyType(type_it, "ISO");

	std::pair<RichPeakListByType::iterator, RichPeakListByType::iterator> p = pks_type.equal_range("ISO");

	size_t count = 0;
	for(RichPeakListByType::iterator iter = p.first; iter != p.second; iter++)
		count++;
	cout << "There are " << count << " ISO peaks" << endl;

	// 2. printPeakList with given order.
	cout << "Sort by ID" << endl;
	peaklist.printPeakList<peak_id>();

	cout << "Sort by MZ" << endl;
	peaklist.printPeakList<peak_mz>();

	cout << "Sort by Area" << endl;
	peaklist.printPeakList<peak_resolution>();

	// 3. get the specified peak.
	//int base_id = peaklist.getBasePeakID();
	cout << "Base RichPeak: " << endl;
	cout << peaklist.getBasePeak<peak_resolution>() << endl;

	cout << "A+1: " << endl;
	RichPeakPtr& pk5 = peaklist.getPeakByShift<peak_resolution>(1);
	cout << pk5 << endl;

	cout << "A+2: " << endl;
	RichPeakPtr& pk6 = peaklist.getPeakByShift<peak_resolution>(2);
	cout << pk6 << endl;

	cout << "A-1: " << endl;
	RichPeakPtr& pk7 = peaklist.getPeakByShift<peak_resolution>(-1);
	cout << pk7 << endl;

	cout << "A-2: " << endl;
	RichPeakPtr& pk8 = peaklist.getPeakByShift<peak_resolution>(-2);
	cout << pk8 << endl;

	// 4. Test Peak.
	//PeakList<Peak, PeakContainer, peak_intensity> peakset;
	RawList peakset;

	PeakPtr pk9 = boost::make_shared<Peak>(471.2, 1320.5);
	peakset.addPeak(pk9);

	PeakPtr pk10 = boost::make_shared<Peak>(472.2, 170.6);
	peakset.addPeak(pk10);

	PeakPtr pk11 = boost::make_shared<Peak>(473.2, 222222.7);
	peakset.addPeak(pk11);

	PeakPtr pk12 = boost::make_shared<Peak>(474.2, 26.4);
	peakset.addPeak(pk12);

	PeakPtr pk13 = boost::make_shared<Peak>(475.2, 32.5);
	peakset.addPeak(pk13);

	// 2. printPeakList with given order.

	cout << "Sort by MZ" << endl;
	peakset.printPeakList<peak_mz>();

	cout << "Sort by Intensity" << endl;
	peakset.printPeakList<peak_intensity>();

	cout << "Base Peak: " << endl;
	cout << peakset.getBasePeak<peak_intensity>() << endl;

	cout << "A+1: " << endl;
	cout << peakset.getPeakByShift<peak_intensity>(1) << endl;

	cout << "A+2: " << endl;
	cout << peakset.getPeakByShift<peak_intensity>(2) << endl;

	cout << "A+3: " << endl;
	cout << peakset.getPeakByShift<peak_intensity>(3) << endl;

	cout << "A-1: " << endl;
	cout << peakset.getPeakByShift<peak_intensity>(-1) << endl;

	cout << "A-2: " << endl;
	cout << peakset.getPeakByShift<peak_intensity>(-2) << endl;

	cout << "A-3: " << endl;
	cout << peakset.getPeakByShift<peak_intensity>(-3) << endl;

}