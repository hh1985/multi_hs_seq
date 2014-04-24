/*
 * =====================================================================================
 *
 *       Filename:  PeriodicTable.h
 *
 *    Description:  Periodic Table loaded from xml file.
 *
 *        Version:  1.0
 *        Created:  4/24/2012 2:53:40 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef  GAG_PERIODICTABLE_H_INC
#define  GAG_PERIODICTABLE_H_INC

#include <GAGPL/CHEMISTRY/Element.h>
#include <GAGPL/MISC/ConfigLoader.h>
#include <boost/noncopyable.hpp>
#include <map>

namespace gag
{
  class PeriodicTable: public ConfigLoader, private boost::noncopyable
  {
    private:
    
			//boost::property_tree::ptree params;
			std::map<std::string, Element> elements;

		protected:

			PeriodicTable() {
				load();
			}
			
		public:

			static PeriodicTable& Instance();
			
			void load(const std::string& filename = "./config/SampleParameterFile.xml");
			
			// Currently we only consider about the situation of retrieving element by symbol.
			const Element getElementBySymbol(const std::string& symbol) const;

      inline const Isotope getIsotope(const std::string& symbol) const
      {
        return this->getIsotopeByRelativeShift(symbol);
      }
      
			size_t getNumIsotope(const std::string symbol) const
      {
        return (*this).getElementBySymbol(symbol).isotopes.size();
      }

			// Get isotope by relative shift from the monoisotopic isotope (default = 0).
			const Isotope getIsotopeByRelativeShift(const std::string& symbol, const size_t shift=0) const;

			// Get isotope by nominal mass.
			const Isotope getIsotopeByNominalMass(const std::string& symbol, const size_t mass) const;

			size_t getMaxRelativeShift(const std::string& symbol) const;
  };
}

#endif   /* ----- #ifndef GAG_PERIODICTABLE_H_INC  ----- */
