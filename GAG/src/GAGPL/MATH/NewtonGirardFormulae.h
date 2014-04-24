/*
 * =====================================================================================
 *
 *       Filename:  NewtonGirardFormulae.h
 *
 *    Description:  Implementation of Newton-Girard formulae, which allows
 *					the conversion between power sum and elementary
 *					symmetric polynomial.
 *
 *        Version:  1.0
 *        Created:  8/26/2012 2:22:42 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu, hh1985@bu.edu
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#ifndef GAG_NEWTONGIRARDFORMULAE_H
#define GAG_NEWTONGIRARDFORMULAE_H

#include <vector>
#include "GAGPL/MATH/Polynomials.h"

namespace msmath
{	
    // Basic principles:
    // 1. The implementation of Newton-Girard formula should be independent of Vietes formula;
    // 2. The class is not supposed to store the values of power sums and ESPs.
    class NewtonGirardFormulae
    {
    public:
        // Constructor. 
	    // NewtonGirardFormulae() {}

        // Preferred constructor
	    NewtonGirardFormulae(size_t var_num)
		    : _var_num(var_num) {}

	    inline void setOrder(size_t n)
	    {
		    _var_num = n;
	    }
	    inline double getOrder() const
	    {
		    return _var_num;
	    }

        // New: This function will automatically fill the missing values from 
        // ps_vec or esp_vec to the other.
        void updateParameters(PowerSumVec& ps_vec, EleSymPolyVec& esp_vec);
        inline void updateParameters(PolyParam& poly_param)
        {
            this->updateParameters(poly_param.second, poly_param.first);
        }

    private:
        // The highest order of the polynomial is determined by _coef.size()-1.
        //const PolyCoef& _coef; 
        size_t _var_num;

    private:
        // New: Updating power sum values given elementary symmetric polynomial.
        void updatePowerSum(PowerSumVec& ps_vec, const EleSymPolyVec& esp_vec);

        // New: Updating elementary symmetric polynomial given power sum values.
        void updateElementarySymmetricPoly(const PowerSumVec& ps_vec, EleSymPolyVec& esp_vec);
    };
}

#endif /* GAG_NEWTONGIRARDFORMULAE_H */





