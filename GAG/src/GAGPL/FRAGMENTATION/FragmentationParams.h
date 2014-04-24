/*
 * =====================================================================================
 *
 *       Filename:  FragmentationParams.h
 *
 *    Description:  Functional group table loaded from xml drive.
 *
 *        Version:  1.0
 *        Created:  04/29/2012  9:51:00 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef GAG_FRAGMENTATIONPARAMS_H
#define GAG_FRAGMENTATIONPARAMS_H

#include <string>
#include <map>
#include <vector>
#include <GAGPL/CHEMISTRY/Composition.h>

namespace gag
{
	//typedef std::pair<int, Composition> CompositionShift;
	struct FragmentationParams
	{
			std::string type;
			std::vector<std::string> cleavage_shift;
			// Dissociation name and corresponding shifts.
			std::map<std::string, std::vector<std::string> > dis_shift;	
			// Internal cleavage might be considered as neutral loss, which reduces the number of fragments.
			// std::map<std::string, std::pair<int, int> > mass_loss; 
	};
}

#endif /* GAG_FRAGMENTATIONPARAMS_H */