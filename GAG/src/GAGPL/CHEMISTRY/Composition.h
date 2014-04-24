/*
 * =====================================================================================
 *
 *       Filename:  Composition.h
 *
 *    Description:  The concept of composition is a linear combination of element mass
 *
 *        Version:  1.0
 *        Created:  4/23/2012 2:51:58 PM
 *       Revision:  none
 *       Compiler:  msvc
 *         Author:  Han Hu, hh.earlydays@gmail.com 
 *   Organization:	Bioinformatics Program, Boston University  
 *
 * =====================================================================================
 */

#ifndef  GAG_COMPOSITION_H_INC
#define  GAG_COMPOSITION_H_INC

#include "GAGPL/CHEMISTRY/Element.h"
#include "GAGPL/CHEMISTRY/PeriodicTable.h"

namespace gag
{
	
	class Composition
	{
		
		private:
			// Periodic table storing the elements information.
			PeriodicTable& _ptable;
  
			//Storing the data.
			std::map<std::string, int> _composition;
      
			// Mass calculated directly from the composition. This is 
			// a general rule except the treatment of Fragmentation.
			// double _mass;

			// The string representative of the composition.
			// std::string _compo_string;
		
		protected:
			

		public:
			// Constructor. 
			Composition(const std::string&);
			Composition();
			
			// Update _composition and _mass
			void update(const std::string&);
			// update _mass
			// void updateMass();

			//void updateString();
			// get _mass
     		double getMass() const;

			double getAverageMass() const;

			// From map structure to string.
			std::string getCompositionString() const;

			void addElement(const std::string&, const int&);

			size_t getMaxNumVariants() const;

			// get _composition.
			inline const std::map<std::string, int>& get() const
			{
				return _composition;
			}

			void add(const Composition&);
			void add(const std::string&);
			
			void deduct(const Composition&);
			void deduct(const std::string&);

			void clear();
			inline bool empty()
			{
				return _composition.empty();
			}

			Composition& operator=(const Composition&);
			friend bool operator==(const Composition& left, const Composition& right)
			{
				return left.getCompositionString()==right.getCompositionString();
			}

			bool operator<(const Composition& compo)
			{
				if(_composition.empty())
					return true;
				std::map<std::string, int>::const_iterator const_iter = _composition.begin();
				for(; const_iter != _composition.end(); const_iter++)
				{
					std::map<std::string, int>::const_iterator it = compo.get().find(const_iter->first);
					if(it != compo.get().end()){
						if(const_iter->second > it->second){
							return false;
						}
					} else {
						return false;
					}
				}
				return true;
			}
		
	};
}




#endif   /* ----- #ifndef GAG_COMPOSITION_H_INC  ----- */
