/*
 * =====================================================================================
 *
 *       Filename:  Unit.h
 *
 *    Description:  The base class for a unit in chemistry.
 *
 *        Version:  1.0
 *        Created:  04/30/2012  3:55:04 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */


#ifndef  GAG_UNIT_H
#define  GAG_UNIT_H

#include <GAGPL/CHEMISTRY/Composition.h>
//#include <GAGPL/CHEMISTRY/FunctionalGroup.h>

namespace gag
{
	// A wrapper of Composition and more.
	class Unit
	{
		public:
			
			// This vector is empty until explicitly added.
			//std::vector<Composition> subgroups;
			Composition compo;
			std::string type;

			Unit()
			{
			}
			
			Unit(const std::string& str)
				: compo(str),type()
			{
			}

			Unit(const Composition& cp)
				: compo(cp),type()
			{
			}

			inline Composition& getComposition()
			{
				return compo;
			}

			inline std::string getCompositionString() const
			{
				return compo.getCompositionString();
			}

			inline double getMass() const
			{
				return compo.getMass();
			}

			inline std::string getType() const
			{
				return type;
			}
			inline void add(const Composition& cp)
			{
				compo.add(cp);
			}
			inline void add(const std::string& str)
			{
				compo.add(str);
			}
			inline void deduct(const Composition& cp)
			{
				compo.deduct(cp);
			}
			inline void deduct(const std::string& str)
			{
				compo.deduct(str);
			}
	};
}

#endif   /* ----- #ifndef GAG_UNIT_H ----- */
