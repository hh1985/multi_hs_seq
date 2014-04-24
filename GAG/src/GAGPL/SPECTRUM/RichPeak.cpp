/********************************************************************
	created:	2012/11/13
	created:	13:11:2012   9:58
	filename: 	RichPeak.cpp
	file path:	GAGPL\SPECTRUM
	file base:	RichPeak
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/SPECTRUM/RichPeak.h"
#include <boost/make_shared.hpp>
#include <iostream>

namespace gag
{
	RichPeakPtr createRichPeak(const double& pk_mz, const double& pk_intensity, const double& pk_resolution, const double& pk_signal_noise)
	{
		static unsigned int pk_id = 1;
		return boost::make_shared<RichPeak>(pk_mz, pk_intensity, pk_id++, pk_resolution, pk_signal_noise);
	}

	void RichPeak::printout( std::ostream& os ) const
	{
		os << "ID: " << id << "\t"
			<< "MZ: " << mz << "\t"
			<< "Intensity: " << intensity << "\t"
			<< "Resolution: " << resolution << "\t"
			<< "S/N: " << signal_noise << std::endl;
	}

}


