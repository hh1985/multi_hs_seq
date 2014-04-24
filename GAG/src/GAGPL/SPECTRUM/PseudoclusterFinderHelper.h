/*
 * =====================================================================================
 *
 *       Filename:  PseudoclusterFinderHelper.h
 *
 *    Description:  This file includes the basic data types for our analysis.
 *
 *        Version:  1.0
 *        Created:  4/2/2012 1:20:20 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef GAG_PSEUDOCLUSTERFINDERHELPER_H
#define GAG_PSEUDOCLUSTERFINDERHELPER_H

#include <map>
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/noncopyable.hpp>

namespace gag
{
				
	// Flags that indicate if a peak is already used in a feature.
	enum Flag {UNUSED, USED};
	enum Order {ACE, DEC};
	
	// This is my customized type of peaks, which includes the status information of a given peak.
	class Peak{
		public:
			Peak() 
				:_mz(), _intensity()
			{
			}
			Peak(const double& mz, const float& intensity)
				: _mz(mz), _intensity(intensity)
			{
			}
			Peak(const Peak& peak)
				: _mz(peak.getMZ()), _intensity(peak.getIntensity())
			{
			}
			// Why virtual destructor.
			virtual ~Peak() {}
			
			inline double getMZ() const {
				return _mz;
			}
			inline void setMZ(double mz) {
				_mz = mz;
			}
			inline float getIntensity() const {
				return _intensity;
			}
			inline void setIntensity(float intensity) {
				_intensity = intensity;
			}
			static inline bool mzSmaller(const Peak& left, const Peak& right)
			{
				return (left.getMZ() <= right.getMZ());
			}
			static bool inline mzLarger(const Peak& left, const Peak& right)
			{
				return (left.getMZ() > right.getMZ());
			}
			static bool intensitySmaller(const Peak& left, const Peak& right)
			{
				return (left.getIntensity() <= right.getIntensity());
			}
			static bool intensityLarger(const Peak& left, const Peak& right)
			{
				return (left.getIntensity() > right.getIntensity());
			}
			
		private:
			double _mz;
			float _intensity;
	};
	
	typedef std::vector<Peak> PeakList;
	
	// Singleton design pattern. Needs to understand it later.
	class Param: boost::noncopyable {
	public:

		static Param& instance() {
			static Param p;
			return p;
		}

		template<typename T>
		void setParameter(const std::string& key, const T& value) {
			params.insert(
					std::make_pair(key, boost::lexical_cast<std::string>(value)));
		}

		template<typename T>
		std::pair<T, bool> getParameter(const std::string& key, const T& defaultValue =
				T()) {
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
		Param() {

		}

	private:
		std::map<std::string, std::string> params;
	};

}

#endif // GAG_PSEUDOCLUSTERFINDERHELPER_H
