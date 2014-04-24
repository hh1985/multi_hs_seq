/*
 * =====================================================================================
 *
 *       Filename:  NewtonGirardFormulaeTest.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  8/26/2012 5:12:38 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu, hh1985@bu.edu 
 *   Organization:  Boston University.
 *
 * =====================================================================================
 */
#include "GAGPL/MATH/NewtonGirardFormulae.h"
#include "GAGPL/MATH/VietesFormulae.h"
#include <iostream>
#include <iterator>

using namespace msmath;
using namespace std;

int main() {
  
    // Create a vector of coefficient for n-order polynomial equation.
    double mycoef[] = {3, 2, 0, 5, 1};
    PolyCoef poly_eq(mycoef, mycoef + sizeof(mycoef) / sizeof(double));

    VietesFormulae viete(poly_eq);
    EleSymPolyVec esp_vec = viete.getElementarySymmetricFunctionFromCoef();

    NewtonGirardFormulae newton(poly_eq.size()-1);
    PowerSumVec ps_vec;
    newton.updateParameters(ps_vec, esp_vec);

    cout << "Power sum size: " << ps_vec.size() << endl;
    copy(ps_vec.begin(), ps_vec.end(), ostream_iterator<double>(cout, ", "));
    // 0 -5, 25, -131, 653

	cout << "\n";
    cout << "ESP size: " << esp_vec.size() << endl;
    copy(esp_vec.begin(), esp_vec.end(), ostream_iterator<double>(cout, ", "));
    // 1, -5, 0, -2, 3

	cout << "\n";
    cout << "************************************************" << endl;
    cout << "*       Test 2: increasing order               *" << endl;
    cout << "************************************************" << endl;

    for(size_t i = 0; i < 3; i++)
        esp_vec.push_back(0.0);

    newton.updateParameters(ps_vec, esp_vec);
    cout << "Power sum size: " << ps_vec.size() << endl;
    copy(ps_vec.begin(), ps_vec.end(), ostream_iterator<double>(cout, ", "));
    // 0, -5, 25, -131, 653, -3300, 16687, -84348

	cout << "\n";
    cout << "ESP size: " << esp_vec.size() << endl;
    copy(esp_vec.begin(), esp_vec.end(), ostream_iterator<double>(cout, ", "));
    // 1, -5, 0, -2, 3, 0, 0, 0
	
	cout << "\n";
    cout << "The End!" << endl;

    cin.get();
}
