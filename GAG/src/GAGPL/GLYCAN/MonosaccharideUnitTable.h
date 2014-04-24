/*
 * =====================================================================================
 *
 *       Filename:  MonosaccharideUnitTable.h
 *
 *    Description:  Load monosaccharide structure information from xml file.
 *
 *        Version:  1.0
 *        Created:  04/30/2012  9:57:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */



#ifndef  GAG_MONOSACCHARIDEUNITTABLE_H
#define  GAG_MONOSACCHARIDEUNITTABLE_H

#include <GAGPL/GLYCAN/Monosaccharide.h>
#include <GAGPL/MISC/ConfigLoader.h>
#include <boost/noncopyable.hpp>

#include <map>

namespace gag
{
	class MonosaccharideUnitTable : public ConfigLoader, private boost::noncopyable
	{
		private:
			std::map<std::string, Monosaccharide> monos;

		protected:

			MonosaccharideUnitTable() {
				load();
			}

		public:

			static MonosaccharideUnitTable& Instance();

			void load(const std::string& filename = "../config/monosaccharide.xml");

			Monosaccharide getMonosaccharideBySymbol(const std::string& symbol) ;

	};
}

#endif   /* ----- #ifndef GAG_MONOSACCHARIDEUNITTABLE_H ----- */
