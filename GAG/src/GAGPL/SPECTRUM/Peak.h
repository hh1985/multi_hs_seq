/********************************************************************
	created:	2012/11/12
	created:	12:11:2012   15:14
	filename: 	Peak.h
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	Peak
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/
#ifndef GAG_PEAK_H
#define GAG_PEAK_H

#include <boost/shared_ptr.hpp>
#include <iostream>


namespace gag
{
	// Forward declaration.
	struct Peak;
	typedef boost::shared_ptr<Peak> PeakPtr;
	
	struct Peak
	{
		double mz;
		double intensity;

		Peak(const double& pk_mz, const double& pk_intensity)
			: mz(pk_mz), intensity(pk_intensity)
		{}

		Peak()
			: mz(0.0), intensity(0.0)
		{}

		static inline bool mzSmaller(const Peak& left, const Peak& right)
		{
			return (left.mz <= right.mz);
		}
		static inline bool mzLarger(const Peak& left, const Peak& right)
		{
			return (left.mz > right.mz);
		}
		static inline bool intensitySmaller(const Peak& left, const Peak& right)
		{
			return (left.intensity <= right.intensity);
		}
		static inline bool intensityLarger(const Peak& left, const Peak& right)
		{
			return (left.intensity > right.intensity);
		}
	
		virtual inline bool empty()
		{
			return mz == 0.0;
		}

		virtual void printout(std::ostream& os) const;

		friend std::ostream& operator<<(std::ostream& os, const PeakPtr& pk)
		{
			pk->printout(os);
			return os;
		}
	};

}

#endif /* GAG_PEAK_H */