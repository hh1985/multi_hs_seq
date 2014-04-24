/*
 * =====================================================================================
 *
 *       Filename:  GlycanSequenceFactory.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  4/23/2012 4:26:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef  GAG_GLYCANSEQUENCEFACTORY_H_INC
#define  GAG_GLYCANSEQUENCEFACTORY_H_INC

#include <GAG/FRAGMENTATION/Fragments.h>
/* Construct a genearal template so that the specific types of glycan can follow the
 * organization */
namespace gag
{
  template <typename T>
  class GlycanSequenceFactory
  {
    private:
      std::set<T> positions;
      std::vector<Monosaccharide> _sequence;
    protected:
      virtual ~GlycanSequenceFactory();
      virtual void addModification(T position, FunctionGroup& functiongroup ) = 0;
  };
}

#endif   /* ----- #ifndef GAG_GLYCANSEQUENCEFACTORY_H_INC  ----- */
