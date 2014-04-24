/********************************************************************
	created:	2012/11/12
	created:	12:11:2012   13:04
	filename: 	Param.cpp
	file path:	GAG\src\GAGPL\MISC
	file base:	Param
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/MISC/Param.h"
#include <iostream>


namespace param
{
	

	void Param::load( const std::string& filename )
	{
		std::cout << "Load it once!" << std::endl;

		using boost::property_tree::ptree;
		ptree pt;

		read_xml(filename, pt);

		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.SystemParameters"))
		{
			if(v.first == "Parameter")
			{
				this->setParameter(v.second.get<std::string>("Name"), 
					v.second.get<std::string>("Value"));
			}
		}

		pt.erase("parameters.SystemParameters");
	}

}