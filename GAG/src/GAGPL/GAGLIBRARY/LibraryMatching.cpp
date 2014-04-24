#include <GAGPL/GAGLIBRARY/LibraryMatching.h>
#include <GAGPL/SPECTRUM/InternalCalibration.h>
#include <iostream>
#include <fstream>
#include <sstream>

namespace gag
{

	void LibraryMatching::load( const std::string& filename, MonoList& mono_list )
	{
		ifstream infile (filename.c_str());
		std::string line;
		if(infile.is_open()) {
			// Remove title.
			getline(infile, line);
			while(getline(infile, line))
			{
				istringstream is;
				is.str(line);
				// The basic format is: m/z	z	intensity.
				double k1; int k2; double k3;
				is >> k1 >> k2 >> k3;
				MonoPeak mp(k1, k2, k3);
				mono_list.addPeak(mp);
			}
		}
		infile.close();
		
	}

	SatelliteReverse LibraryMatching::matchLibrary( MonoList& mono_list )
	{
		// Create the mapping relationship between MonoPeak and Assignments and error.
		MonoPeakListByIntensity& mono_pk_int = mono_list.getPeakListByType<peak_intensity>();

		std::map<MonoPeakPtr, double> peak_error;

		this->matchLibraryByMassList(mono_list, peak_error);

		// Adjust the mz value.
		std::map<double, MonoPeakPtr> adjusted_results = adjustMZ(mono_list, peak_error);

		// Redo the matching process.
		return this->matchLibraryByAdjustedMZ(adjusted_results);
		
	}

	SatelliteReverse LibraryMatching::matchLibraryByMass( double mass, double error_range )
	{
		SatelliteReverse matched_results;

		// consider the possibilities of different modifications.

		SatelliteReverse::iterator sr_iter = sr_table.lower_bound(mass);

		// Some key values are larger than mass.
		if(sr_iter != sr_table.end())
		{
			// Move the iterator upwards or downwards to search for closest element.
			for(SatelliteReverse::iterator temp_iter = sr_iter; temp_iter != sr_table.end(); temp_iter++) {
				double distance = abs(mass - temp_iter->first)/mass;
				if(distance < error_range) {
					matched_results.insert(*temp_iter);
					continue;
				} 

				break;
			}

			for(SatelliteReverse::reverse_iterator rev_temp_iter(sr_iter); rev_temp_iter != sr_table.rend(); rev_temp_iter++)
			{
				double distance = abs(mass - rev_temp_iter->first)/mass;
				if(distance < error_range) {
					matched_results.insert(*rev_temp_iter);
					continue;
				}

				break;
			}
		}

		return matched_results;
	}

	void LibraryMatching::matchLibraryByMassList( MonoList& mono_list, std::map<MonoPeakPtr, double>& peak_error)
	{
		if(peak_error.size() != 0)
			peak_error.clear();

		MonoPeakListByIntensity& mono_pk_int = mono_list.getPeakListByType<peak_intensity>();

		for(MonoPeakListByIntensity::iterator mono_iter = mono_pk_int.begin(); mono_iter != mono_pk_int.end(); mono_iter++)
		{
			MonoPeakPtr pk = *mono_iter;
			double mz = pk->mz; int z = pk->z;
			double mass = calculateMass(mz, z);

			SatelliteReverse results = matchLibraryByMass(mass, 1e-5);
			
			if(results.size() == 0)
				continue;
			else {
				SatelliteReverse::iterator iter = results.begin();
				for(; iter != results.end(); iter++)
				{
					TreeNodePtr tn_ptr = iter->second.first;
					const SatelliteReverse& sm = tn_ptr->satellites.at(iter->second.second);
					double theo_mz = calculateMZ(sm.second, z);
					double exp_ppm = (mz - theo_mz)/theo_mz * 1e6;
					peak_error.insert(std::make_pair(MonoPeakPtr, exp_ppm));
				}
			}

		}
	}
 
	std::map<double, MonoPeakPtr> LibraryMatching::adjustMZ( MonoList& mono_list, std::map<MonoPeakPtr, double>& peak_error )
	{
		using namespace alglib;

		real_2d_array xy;
		int npoints = peak_error.size(); int nvars = 1;
		xy.setlength(npoints, nvars+1);

		for(int i=0; i<npoints; i++)
		{
			// x should be the 
			xy(i,0) = peak_error[i]->first->mz; xy(i,1) = peak_error[i]->second;
		}

		// 1. Construct model
		linearmodel lm;
		int info_code = 0; lrreport ar;
		lrbuild(xy, npoints, nvars, info_code, lm, ar);

		// 2. Correct the m/z of MonoList.
		std::map<double, MonoPeakPtr> corrected_results;
		MonoPeakListByMZ& mono_list_mz = mono_list.getPeakListByType<peak_mz>();
		for(MonoPeakListByMZ::iterator iter = mono_list_mz.begin(); iter != mono_list_mz.end(); iter++)
		{
			MonoPeakPtr pk = *iter;
			real_1d_array x;
			x.setlength(1); x = pk->mz;
			double pre_error = lrprocess(lm, x);
			corrected_results.insert(pk->mz - pre_error, MonoPeakPtr);
		}

		return corrected_results;
	}

	SatelliteReverse LibraryMatching::matchLibraryByAdjustedMZ( std::map<double, MonoPeakPtr>& adjusted_mz )
	{
		SatelliteReverse matched_results;

		for(std::map<double, MonoPeakPtr>::iterator iter = adjusted_mz.begin(); iter != adjusted_mz.end(); iter++)
		{
			MonoPeakPtr pk = *iter;
			double mz = pk->mz; int z = pk->z;
			double mass = calculateMass(mz, z);

			SatelliteReverse results = matchLibraryByMass(mass, 1e-6);

			if(results.size() == 0)
				continue;
			else {
				SatelliteReverse::iterator iter = results.begin();
				for(; iter != results.end(); iter++)
				{
					TreeNodePtr tn_ptr = iter->second.first;
					const SatelliteReverse& sm = tn_ptr->satellites.at(iter->second.second);
					// Append the information into the overall pool.
					SatelliteReverse::const_iterator const_iter = sm.begin();
					for(; const_iter != sm.end(); const_iter++)
					{
						matched_results.insert(*const_iter);
					}
				}
			}
		}

		return matched_results;
	}



}