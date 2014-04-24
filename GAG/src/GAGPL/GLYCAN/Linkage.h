/*
 * =====================================================================================
 *
 *       Filename:  Linkage.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/ 4/2012  7:41:08 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */


#ifndef  GAG_LINKAGE_INC
#define  GAG_LINKAGE_INC

#include <string>
#include <GAGPL/CHEMISTRY/FunctionalGroup.h>

namespace gag
{
	struct Linkage
	{
		size_t nre_id;
		size_t start;
		size_t end;
		std::string type;

		Linkage() {}
		Linkage(const size_t& id, const size_t& s1, const size_t& s2, const std::string& tp)
			: nre_id(id),start(s1),end(s2),type(tp) 
		{}
	};

}

#endif   /* ----- #ifndef GAG_LINKAGE_INC ----- */
