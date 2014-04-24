/********************************************************************
	created:	2013/07/04
	created:	4:7:2013   14:24
	filename: 	Polynomials.h
	file path:	GAGPL\MATH
	file base:	Polynomials
	file ext:	h
	author:		Han Hu
	
	purpose:	The header file stores the basic terms in polynomials.
*********************************************************************/

#ifndef GAG_POLYNOMIALS_H
#define GAG_POLYNOMIALS_H

#include <vector>

namespace msmath
{
	// The coefficients is ordered, from a(0), a(1), ... to a(n). The coefficients uniquely defines a polynomial equation.
	typedef std::vector<double> PolyCoef;

	// The vector of elementary symmetric polynomial.
	typedef std::vector<double> EleSymPolyVec;

	// The vector of power sum.
	typedef std::vector<double> PowerSumVec;

	// The parameter pair of esp and ps.
	typedef std::pair<EleSymPolyVec, PowerSumVec> PolyParam;
}

#endif /* GAG_POLYNOMIALS_H */
