/********************************************************************
	created:	2012/11/12
	created:	12:11:2012   15:11
	filename: 	RichPeak.h
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	RichPeak
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef GAG_RICHPEAK_H
#define GAG_RICHPEAK_H

#include "GAGPL/SPECTRUM/Peak.h"
#include <boost/make_shared.hpp>

namespace gag
{
	// Forward declaration.
	struct RichPeak;
	typedef boost::shared_ptr<RichPeak> RichPeakPtr;

	struct RichPeak: public Peak
	{
		unsigned int id;
		double resolution;
		double signal_noise;
		int charge_state;

		bool pk_status;

		// ENV, TBD, NOISE
		std::string pk_type;
		
		// To facilitate the management of id.
		friend RichPeakPtr createRichPeak(const double& pk_mz, const double& pk_intensity, const double& pk_resolution, const double& pk_signal_noise);

		RichPeak() 
			: id(0), resolution(0.0), signal_noise(0.0) {}

		inline bool operator<(RichPeak& pk)
		{
			return id < pk.id;
		}

		RichPeak(const double& pk_mz, const double& pk_intensity, unsigned int pk_id, const double& pk_resolution, const double& pk_signal_noise)
			: Peak(pk_mz, pk_intensity), id(pk_id), resolution(pk_resolution), signal_noise(pk_signal_noise), pk_status(true), pk_type("TBD")
		{
		}

	private:		
		virtual void printout(std::ostream& os) const;

		friend std::ostream& operator<<(std::ostream& os, const RichPeakPtr& pk)
		{
			pk->printout(os);
			return os;
		}

	};

	
}


#endif /* GAG_RICHPEAK_H */