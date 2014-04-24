/*
 * =====================================================================================
 *
 *       Filename:  Spectrum.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  9/6/2012 2:14:17 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu(hh1985@bu.edu), 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef  GAG_SPECTRUM_H_INC
#define  GAG_SPECTRUM_H_INC

#include <GAGPL/SPECTRUM/PeakList.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/shared_ptr.hpp>

namespace gag
{
	//enum Flag {USED, NEW};
  
	//class Envelop;

	class RawPeak
  {
	private:
		size_t id;
    Peak pk;
		float area;
		float relative_area;
		// Signal over noise.
    float s_n;
		// Store the charge states which has been estimated.
		std::set<int> charge_states;
		//Flag status;
		//boost::ptr_vector<Envelop> env_ptr_vec;
	public:
		inline size_t getID() const
		{	return id;	}

		inline double getMZ() const
		{	return pk.mz;	}

		inline double getIntensity() const
		{	return pk.intensity; }

		inline float getSignalOverNoise() const
		{	return s_n; }

		inline float getPeakArea() const
		{	return area;	}

		inline float getRelativePeakArea() const
		{	return relative_area; }

		inline void setRelativePeakArea(float area_percentage)
		{	relative_area = area_percentage;	}
		
		// Decide if the peak has been considered as a candidate base peak.
		bool estimated(int charge);

		static inline bool mzSmaller(const RawPeak& left, const RawPeak& right)
		{
			return Peak::mzSmaller(left.pk, right.pk);
		}
		static inline bool mzLarger(const RawPeak& left, const RawPeak& right)
		{
			return Peak::mzLarger(left.pk, right.pk);
		}
		static bool intensitySmaller(const RawPeak& left, const RawPeak& right)
		{
			return Peak::intensitySmaller(left.pk, right.pk);
		}
		static bool intensityLarger(const RawPeak& left, const RawPeak& right)
		{
			return Peak::intensityLarger(left.pk, right.pk);
		}
		static bool SNSmaller(const RawPeak& left, const RawPeak& right)
		{
			return (left.s_n <= right.s_n);
		}
		static bool SNLarger(const RawPeak& left, const RawPeak& right)
		{
			return (left.s_n > right.s_n);
		}
  };

	typedef multi_index_container<
		boost::shared_ptr<RawPeak>,
		indexed_by<
		hashed_index<mem_fun<RawPeak, size_t, &RawPeak::getID> >,
		// sort by less<double> on m/z.
		ordered_unique<member<RawPeak, double, &RawPeak::getMZ> >,
		// sort by greater<double> on intensity.
		ordered_unique<member<RawPeak, double, &RawPeak::getIntensity>, std::greater<double> >,
		// sort by less<float> on signal_over_noise.
		ordered_non_unique<member<RawPeak, float, &RawPeak::getSignalOverNoise> >
		>
	> PeakListPtr;



	typedef PeakListPtr::nth_index<0>::type PeakListByID;
	typedef PeakListPtr::nth_index<1>::type PeakListByMZ;
	typedef PeakListPtr::nth_index<2>::type PeakListByIntensity;

	// The spectrum is the raw type of a vector of peaks.
	typedef std::vector<RawPeak> Spectrum;

	//typedef std::vector<RawPeak> Spectrum;
	//typedef std::vector<Spectrum> SpectrumList;

}



#endif   /* ----- #ifndef GAG_SPECTRUM_H_INC  ----- */
