/********************************************************************
	created:	2012/11/16
	created:	16:11:2012   16:41
	filename: 	Envelop.cpp
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	Envelop
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/SPECTRUM/Envelop.h"

namespace gag
{
	//EnvelopPtr createEnvelop(int charge, const EnvelopBoundary& bound)
	//{
	//	static unsigned int env_id = 0;
	//	return boost::make_shared<Envelop>(env_id++, charge, bound);
	//}
	EnvelopPtr createEnvelop(int charge)
	{
		static unsigned int env_id = 0;
		return boost::make_shared<Envelop>(env_id++, charge);
	}

	std::map<int, PeakPtr> Envelop::getTheoreticalPeaks()
	{
		std::map<int, PeakPtr> pk_map;

		// Currently, the program only process envelops with positive shift.
		int shift = 0;
		while(1)
		{
			PeakPtr pk = this->getTheoreticalPeak(shift);
			if(pk->mz != 0.0)
			{
				pk_map.insert(std::make_pair(shift, pk));
			} else {
				break;
			}
			shift++;
		}
		return pk_map;

	}

	PeakPtr Envelop::getTheoreticalPeak( int shift )
	{
		PeakPtr pk = theo_dist.getPeakByShift<peak_intensity>(shift);
		return pk;
	}


}