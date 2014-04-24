/********************************************************************
	created:	2012/08/27
	created:	27:8:2012   14:18
	filename: 	VietesFormulae.h
	file path:	GAGPL\MATH
	file base:	VietesFormulae
	file ext:	h
	author:		Han Hu, hh1985@bu.edu

	purpose:	The class for applying Viete's formulae.
				For a general polynomial of degree n
				P(x) = a(0) + a(1)x + a(2)x^2 + ... a(n-1)x^(n-1) + a(n)x^n
				The relationship between the roots and coefficients is:
				x(1) + x(2) + ... + x(n-1) + x(n) = -a(n-1)/a(n)
				...
			    x(1)x(2) ... x(n) = (-1)^n * a(0)/a(n)
*********************************************************************/
#ifndef GAG_VIETESFORMULAE_H
#define GAG_VIETESFORMULAE_H

#include "GAGPL/MATH/Polynomials.h"

namespace msmath
{
	class VietesFormulae
	{
	private:
		PolyCoef& _coef;

	public:
		VietesFormulae(PolyCoef& coef)
			: _coef(coef) {}

		inline void setCoef(PolyCoef& coef)
		{
			_coef = coef;
		}

        // Generating the vector of elementary symmetric polynomial from 
        // coefficients of polynomial.
		EleSymPolyVec getElementarySymmetricFunctionFromCoef();
	};

}

#endif /* GAG_VIETESFORMULAE_H */