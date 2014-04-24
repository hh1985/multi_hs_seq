/********************************************************************
	created:	2012/08/27
	created:	27:8:2012   15:04
	filename: 	VietesFormulaeTest.cpp
	file path:	test
	file base:	VietesFormulaeTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/MATH/VietesFormulae.h"
#include <boost/foreach.hpp>
#include <iostream>

int main()
{
	using namespace msmath;
	using namespace std;

    double mycoef[] = {3, 2, 0, 5, 1};
    PolyCoef poly_eq(mycoef, mycoef + sizeof(mycoef) / sizeof(double));

    VietesFormulae formulae(poly_eq);
    EleSymPolyVec ele_values = formulae.getElementarySymmetricFunctionFromCoef();

    // The correct value should be: 1, -5, 0, -2, 3.
    for(size_t i=0; i<ele_values.size(); i++)
    {
        cout << i << ": " << ele_values.at(i) << endl;
    }

    cout << "The End!" << endl;

    cin.get();

}