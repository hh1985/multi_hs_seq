/*
 * =====================================================================================
 *
 *       Filename:  ConfigLoader.h
 *
 *    Description:  This file is for loading configuration file. This should only 
 *    							include the interface for accessing xml file. Originally implemented 
 *    							using boost_property_tree, but now consider using tinyXML2.
 *
 *        Version:  1.0
 *        Created:  04/26/2012 10:06:16 AM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef GAG_CONFIGLOADER_H
#define GAG_CONFIGLOADER_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

namespace gag
{
	class ConfigLoader 
	{
		public:
			// Abstract interface: loading XML parameter file.
			virtual void load(const std::string& filename) = 0;
			
			virtual ~ConfigLoader() {}

	};

}
#endif /* GAG_CONFIGLOADER_H */
