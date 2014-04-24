/********************************************************************
	created:	2012/11/10
	created:	10:11:2012   15:10
	filename: 	Param.h
	file path:	GAGPL\MISC
	file base:	Param
	file ext:	h
	author:		Han Hu
	
	purpose:	Basic module for defining parameters used in this project.
*********************************************************************/

#ifndef GAG_PARAM_H
#define GAG_PARAM_H

#include <boost/noncopyable.hpp>
#include <boost/lexical_cast.hpp>
#include "GAGPL/MISC/ConfigLoader.h"
#include <string>
#include <map>

namespace param
{
	// Singleton design pattern.
	class Param: boost::noncopyable, gag::ConfigLoader {
	public:

		static Param& Instance(const std::string& filename = "./config/systemparameters.xml") {
			static Param p(filename);
			return p;
		}

		template<typename T>
		void setParameter(const std::string& key, const T& value) {
			// Test if the parameter exist.
			if(params.find(key) != params.end())
				params[key] = boost::lexical_cast<std::string>(value);
			else
				params.insert(
					std::make_pair(key, boost::lexical_cast<std::string>(value)));
		}

		template<typename T>
		std::pair<T, bool> getParameter(const std::string& key, const T& defaultValue = T()) {
			std::map<std::string, std::string>::const_iterator i = params.find(key);
			if (i != params.end()) {
				try {
					return std::make_pair(boost::lexical_cast<T>(i->second), true);
				} catch (...) {
					//
					return std::make_pair(defaultValue, false);
				}
			} else {
				return std::make_pair(defaultValue, false);
			}
		}

	private:
		void load(const std::string& filename);

		Param(const std::string& filename) {
			load(filename);
		}

	private:
		// Storing the parameters information.
		std::map<std::string, std::string> params;
	};
}


#endif /* GAG_PARAM_H */