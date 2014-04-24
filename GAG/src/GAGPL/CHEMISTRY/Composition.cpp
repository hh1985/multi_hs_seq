/*
 * =====================================================================================
 *
 *       Filename:  Composition.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  4/23/2012 10:12:49 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  HAN HU, hh1985@bu.edu
 *   Organization:  Bioinformatics Program, Boston University
 *
 * =====================================================================================
 */
#include "GAGPL/CHEMISTRY/Composition.h"
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <locale>

namespace  gag
{
	double Composition::getMass() const
	{
		std::map<std::string, int>::const_iterator iter1 = _composition.begin();
		double mass = 0.0;
		for(; iter1 != _composition.end(); iter1++)
		{
			const Isotope& iso = _ptable.getIsotopeByRelativeShift(iter1->first);
			// TBD: support for multiple isotopes.
			// std::vector<Isotope>::iterator iter2 = e.isotopes.begin();
			mass += (*iter1).second * iso.mass;
		}
		return mass;
	}

	std::string Composition::getCompositionString() const
	{
		std::string compo_string;

		// Generate the string.
		for(std::map<std::string, int>::const_iterator iter = _composition.begin(); iter != _composition.end(); iter++)
		{
			if(iter->second == 1)
				compo_string.append(iter->first);
			else if(iter->second > 1)
			{
				compo_string.append(iter->first);
				compo_string.append(boost::lexical_cast<std::string>(iter->second));
			} else {
				throw 32751;
			}
		}
		return compo_string;
	}

	void Composition::addElement(const std::string& elem, const int& num)
	{
		// If the string has been contained in _composition, simply add it up.
		// otherwise create a new one.
		std::map<std::string, int>::iterator iter = _composition.find(elem);
		
		if(iter != _composition.end())
		{
			iter->second += num;
		} else 
		{
			_composition.insert(std::make_pair(elem, num));
		}

		//updateMass();
		//updateString();
	}

	void Composition::update(const std::string& str)
  {
    // Split the composition into a linear combination of Elements.
    // and do format examination.
		if(str.length() == 0)
			throw 8888;
		
		std::locale loc;
		std::string temp_str;
		std::string temp_num;

		for(size_t i=0; i<str.length(); i++)
		{
			// Only three types of character are allowed:
			// Uppercase letter, lowercase letter and digits (0-9)

			// Better to add exception control of the data type here.
			
			if(std::isupper(str[i],loc))
			{
				temp_str=str[i];

				// Dump the string and number into composition.
				if(i == str.length()-1 || std::isupper(str[i+1],loc))
				{
					(*this).addElement(temp_str, 1);

					temp_str.clear();
					//temp_num.clear();

				} 

			} else if(std::islower(str[i],loc)) // There will always be only one lowercase letter.
			{

				temp_str.push_back(str[i]);
				
				// If the lowercase letter is the last one or the next is not a digit (should be uppercase letter)
				if(i == str.length()-1 || isupper(str[i+1],loc)) 
				{
					(*this).addElement(temp_str, 1);
					temp_str.clear();
				}
				
			} else if(std::isdigit(str[i],loc))
			{
				temp_num.push_back(str[i]);

				// It is unlikely that temp_num still has length 0. But be careful of that.
				if(i == str.length()-1 || isupper(str[i+1],loc)) 
				{
					(*this).addElement(temp_str, boost::lexical_cast<int>(temp_num));
					
					temp_str.clear();
					temp_num.clear();
				}
			}
		}

		//updateMass();
		//updateString();

		//return _composition;

  }
  
	void Composition::add(const Composition& cp)
  {
		
		const std::map<std::string, int>& cp_map = cp.get();
    
		for(std::map<std::string, int>::const_iterator iter1 = cp_map.begin(); iter1 != cp_map.end(); iter1++)
    {
      
			std::map<std::string, int>::iterator iter2 = _composition.find((*iter1).first);

      if(iter2 != _composition.end())
				(*iter2).second += (*iter1).second;
      
			else
				_composition.insert(std::make_pair(iter1->first,iter1->second));
    }

		//updateMass();
		//updateString();

  }
	
	void Composition::deduct(const Composition& cp)
  {

    const std::map<std::string, int>& cp_map = cp.get();

		for(std::map<std::string, int>::const_iterator iter1 = cp_map.begin(); iter1 != cp_map.end(); iter1++)
    {

			std::map<std::string, int>::iterator iter2 = _composition.find((*iter1).first);
      
			if(iter2 != _composition.end() && (*iter1).second < (*iter2).second)
        (*iter2).second -= (*iter1).second;
      else if((*iter1).second == (*iter2).second)
        _composition.erase(iter2);
      else
				throw std::runtime_error("Error: the deducted composition is invalid!");
    }

		//updateMass();
		//updateString();
  }

  // Update the value of _compostion.
  void Composition::add(const std::string& compo)
  {
		Composition cp(compo);
		this->add(cp);    

  }
 
	void Composition::deduct(const std::string& compo)
  {
		Composition cp(compo);
		this->deduct(cp);
  }

	void Composition::clear()
	{
		_composition.clear();
		//_mass = 0.0;
		//_compo_string.clear();
	}

	Composition::Composition(const std::string& str)
		: _ptable(PeriodicTable::Instance()), _composition()/*, _mass(0.0), _compo_string() */
	{
		// Update _composition and _mass
		//_ptable.load();
		(*this).update(str);
	}

	Composition::Composition()
		: _ptable(PeriodicTable::Instance()), _composition()/*, _mass(0.0), _compo_string()*/
	{
		//_ptable.load();
	}

	Composition& Composition::operator=(const Composition& rhs)
	{
		if(this != &rhs)
		{
			_composition = rhs._composition; // Copy constructor for map.
			//_mass = rhs._mass;
		}

		return *this;
	}

	size_t Composition::getMaxNumVariants() const
	{
		std::map<std::string, int>::const_iterator it = _composition.begin();

		size_t max_num = 0;

		for(; it!=_composition.end(); it++)
		{
			max_num += it->second * _ptable.getMaxRelativeShift(it->first);
		}

		return max_num;
	}

	double Composition::getAverageMass() const
	{
		double avg_mass(0.0);
		std::map<std::string, int>::const_iterator it = _composition.begin();

		for(; it!=_composition.end(); it++)
		{
			avg_mass += _ptable.getElementBySymbol(it->first).getAverageMass() * it->second;
		}

		return avg_mass;
	}
}

