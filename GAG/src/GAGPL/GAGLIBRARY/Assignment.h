/********************************************************************
	created:	2014/04/25
	created:	25:4:2014   16:04
	filename: 	Assignment.h
	file path:	GAG\src\GAGPL\GAGLIBRARY
	file base:	Assignment
	file ext:	h
	author:		Han Hu
	
	purpose:	Assignment describes the structural information regarding 
            sulfate distribution: candidate modification sites and 
            modification number.
*********************************************************************/

#ifndef GAG_ASSIGNMENT_H
#define GAG_ASSIGNMENT_H

#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <GAGPL/SPECTRUM/MonoPeak.h>

namespace gag
{
  // Only concern about the neutral loss, 
  // It has nothing to do with the sulfate and acetate number.

  typedef std::map<std::string, int> NeutralLoss;

  class Assignment: public Unit {

  };
}


#endif /* ASSIGNMENT_H */