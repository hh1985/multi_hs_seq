/********************************************************************
	created:	2013/04/19
	created:	19:4:2013   22:56
	filename: 	MassConversion.cpp
	file path:	GAGPL\MATH
	file base:	MassConversion
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include <GAGPL/MATH/MassConversion.h>

namespace msmath
{
	double calculateMass( double mz, int charge)
	{
		// int mode = param.getParameter<int>("mode").first;

		// Be careful of the ion mode.
		Param& param = Param::Instance();
		double electron_mass = param.getParameter<double>("electron_mass").first;

		return abs(charge) * mz - charge * (Composition("H").getMass() -electron_mass);
	}

	double calculateMZ(double mass, int charge, int pre_charge /* = 0 */)
	{
		Param& param = Param::Instance();
		double electron_mass = param.getParameter<double>("electron_mass").first;

		int coef_h = (pre_charge == 0 ? charge : pre_charge);
		// TBD: this should be controlled by parameter.
		// int coef_e = (pre_charge == 0 ? 1 : 0);
		return (mass + coef_h * (Composition("H").getMass() - electron_mass))/abs(charge);
	}
}