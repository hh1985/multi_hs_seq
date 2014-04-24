#ifndef GAG_MONOPEAK_H
#define GAG_MONOPEAK_H

#include "GAGPL/SPECTRUM/Peak.h"
#include <boost/make_shared.hpp>

namespace gag
{
	struct MonoPeak;
	typedef boost::shared_ptr<MonoPeak> MonoPeakPtr;

	struct MonoPeak: public Peak
	{
		int z;

		MonoPeak(const double& mz, const double& intensity, const int& charge_state)
			: Peak(mz, intensity), z(charge_state) {}
	};
}

#endif /* GAG_MONOPEAK_H */