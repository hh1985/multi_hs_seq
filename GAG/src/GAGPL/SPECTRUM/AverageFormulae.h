/*
 * =====================================================================================
 *
 *       Filename:  AverageFormulae.h
 *
 *    Description:  Class for calculating average formulae. Similar to averagine 
 *    algorithm, this class provide a way of estimating theoretical element composition 
 *    based on unit Da.
 *
 *        Version:  1.0
 *        Created:  9/12/2012 1:44:45 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (hh1985), 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef  GAG_AverageFormulae_H
#define  GAG_AveageFormulae_H

#include <GAGPL/CHEMISTRY/Composition.h>

namespace gag
{
	class AverageFormulae
	{
	private:
		std::map<std::string, float> avg_compo;

	public:
		AverageFormulae() {}
		AverageFormulae(const Composition& compo)
		{ this->update(compo);}
		
		// Calculate the avg_compo using given standard composition.
		void update(const Composition& compo);

		// Get the theoretical composition for current mass.
		Composition getComposition(double mass);
	};
}

#endif   /* ----- #ifndef AverageFormulae_INC  ----- */

