/********************************************************************
	created:	2012/11/13
	created:	13:11:2012   14:54
	filename: 	Peak.cpp
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	Peak
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/SPECTRUM/Peak.h"
#include <iostream>

namespace gag
{
	void Peak::printout(std::ostream& os) const
	{
		std::cout.precision(10);
		os << "MZ: " << mz << "\t"
			<< "Intensity: " << intensity << "\t" << std::endl;
	}
}