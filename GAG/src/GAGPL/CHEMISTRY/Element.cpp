/********************************************************************
	created:	2013/11/07
	created:	7:11:2013   23:56
	filename: 	Element.cpp
	file path:	GAGPL\CHEMISTRY
	file base:	Element
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/CHEMISTRY/Element.h"

namespace gag
{

	double Element::getAverageMass() const
	{
		double avg_mass = 0.0;
		const IsotopeSetSequential& iso_index = isotopes.get<0>();

		for(IsotopeSetSequential::const_iterator iter = iso_index.begin(); iter != iso_index.end(); iter++)
		{
			avg_mass += iter->mass * iter->abundance;
		}

		return avg_mass;
	}

}