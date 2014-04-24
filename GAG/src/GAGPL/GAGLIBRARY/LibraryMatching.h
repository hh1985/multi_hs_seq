/********************************************************************
	created:	2013/04/05
	created:	5:4:2013   16:25
	filename: 	LibraryMatching.h
	file path:	GAGPL/GAGLIBRARY
	file base:	LibraryMatching
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef GAG_LIBRARYMATCHING_H
#define GAG_LIBRARYMATCHING_H

#include <GAGPL/SPECTRUM/InternalCalibration.h>
#include <GAGPL/SPECTRUM/PeakList.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/MISC/Param.h>

namespace gag
{
	class LibraryMatching
	{
	public:
		//LibraryMatching() {}
		LibraryMatching(LibraryTree& libtree)
			: sr_table(libtree.getReverseTable()), param(Param::Instance()) 
		{}

	public:
		// load the deconvoluted results. The input is the file, and the output is the richlist.
		// The monoisotopic peak list can either be loaded from file or read from RichList object created by SimpleFinder.
		void load(const std::string& filename, MonoList& mono_list);

		// There are several factors to consider: the number of possible sulfate and acetate groups.
		SatelliteReverse matchLibraryByMass(double mass, double error_range);
		void matchLibraryByMassList(MonoList& mono_list, std::map<MonoPeakPtr, double>& peak_error);
		SatelliteReverse matchLibraryByAdjustedMZ(std::map<double, MonoPeakPtr>& adjusted_mz);

	private:
		SatelliteReverse& sr_table;
		Param& param;
		
		// This function includes two part: 
		// 1. try to find specified types of ions.
		// 2. adjust the shift degree to minimize the error.
		SatelliteReverse matchLibrary(MonoList& mono_list);

		// The cost of modifying multi-index container's key value is too high.
		std::map<double, MonoPeakPtr> adjustMZ(MonoList& mono_list, std::map<MonoPeakPtr, double>& peak_error);

	};
}



#endif /* GAG_LIBRARYMATCHING_H */

