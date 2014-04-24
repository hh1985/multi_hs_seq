/*
 * =====================================================================================
 *
 *       Filename:  Elements.h
 *
 *    Description:  This class stores the basic information of elements including
 *    							Symbol, mass, Alias name and so on. This information can come from 
 *    							config files storing in a separate element.conf file.
 *
 *        Version:  1.0
 *        Created:  4/22/2012 3:06:04 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef  GAG_ELEMENT_H_INC
#define  GAG_ELEMENT_H_INC

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/member.hpp>
#include <string>

namespace gag
{
  struct Isotope
  {
    double mass;
    float abundance;
		size_t neutron_shift;
		size_t nominal_mass;
  };
	// Multiply indexed Isotope set.
	using namespace boost;
	using namespace boost::multi_index;
	typedef multi_index_container<
		Isotope,
		indexed_by<
			sequenced<>, // index #0
			// Sort by less<size_t> on neuron_shift
			ordered_unique<member<Isotope, size_t, &Isotope::neutron_shift> >,
			// Sort by less<size_t> on rounded_mass
			ordered_unique<member<Isotope, size_t, &Isotope::nominal_mass> >
		>
	> IsotopeSet;

	typedef IsotopeSet::nth_index<0>::type IsotopeSetSequential;
	typedef IsotopeSet::nth_index<1>::type IsotopeSetByShift;
	typedef IsotopeSet::nth_index<2>::type IsotopeSetByMass;
  
	struct Element 
	{
		std::string symbol;
		std::string name;
		int atomicity;

		IsotopeSet isotopes;
		double getAverageMass() const;
	};
}

#endif   /* ----- #ifndef GAG_ELEMENT_H_INC  ----- */
