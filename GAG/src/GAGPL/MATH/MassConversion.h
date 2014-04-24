/********************************************************************
	created:	2013/04/15
	created:	15:4:2013   0:31
	filename: 	MassConversion.h
	file path:	GAGPL\MATH
	file base:	MassConversion
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef GAG_MASSCONVERSION_H
#define GAG_MASSCONVERSION_H

#include <GAGPL/MISC/Param.h>
#include <GAGPL/CHEMISTRY/Composition.h>

using namespace param;
using namespace gag;


namespace msmath
{
	double calculateMass( double mz, int charge);

	double calculateMZ(double mass, int charge, int pre_charge /* = 0 */);
}


#endif /* GAG_MASSCONVERSION_H */