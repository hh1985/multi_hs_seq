/*
 * =====================================================================================
 *
 *       Filename:  AverageFormulae.cpp
 *
 *    Description:  Implementation for calculating average formulae.
 *
 *        Version:  1.0
 *        Created:  9/12/2012 1:49:50 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (hh1985@bu.edu), 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#include <math.h>

namespace gag
{

	void AverageFormulae::update(const Composition& compo)
	{
		if(!avg_compo.empty())
			avg_compo.clear();

		// Calculate the mass of the composition.
		double total_mass = compo.getMASS();
		
		double factor = total_mass;

		const std::map<std::string, int>& internal_compo = compo.get();

		std::map<std::string, int>::const_iterator const_iter = internal_compo.begin();
	
		for(; const_iter != internal_compo.end(); const_iter++)
			avg_compo.insert(std::make_pair(const_iter->first, const_iter->second / factor));
		
	}

	Composition AverageFormulae::getComposition(double mass)
	{
		if(avg_compo.empty())
			throw std::runtime_error("The average composition has to be defined!");

		Composition temp_compo;

		std::map<std::string, int>::iterator iter = avg_compo.begin();

		for(; iter != avg_compo.end(); avg_compo++)
			temp_compo.insert(std::make_pair(iter->first, std::ceil(iter->second * mass)));

		return temp_compo;
	}
}
