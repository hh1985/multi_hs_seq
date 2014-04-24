/*
 * =====================================================================================
 *
 *       Filename:  PseudoclusterFinder.h
 *
 *    Description:  This class is used to find pseudoclusters: peak series with given
 *    		    			even space. This finder is designed for glycan, which receives less 
 *    		    			limit from averagine but more benefits from theoretical library. The
 *  		    				theoretical library can be considered to be a quickfinder, like what
 *    		    			Yu has done.
 *
 *        Version:  1.0
 *        Created:  3/30/2012 4:21:24 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */


#ifndef GAG_PSEUDOCLUSTERFINDER_H
#define GAG_PSEUDOCLUSTERFINDER_H

#include <GAGPL/SPECTRUM/PseudoclusterFinderHelper.h>

namespace gag
{	
	/*  Define the basic structure of a peak */
	/*  The status member is used to indicate the including of pseudocluster
	 *  */
	/*  The spectrum should not include scan information since it is based on the average. */
	// The format of the input.
	
	// Enhanced version of the peaktype for isotope cluster.
	class Pseudocluster;
	typedef std::vector<Pseudocluster> PseudoclusterList;
	
	class IsotopePeak: public Peak
	{
		// The location information includes cluster id and shift from highest peak.  The shift can be either left or right, and the value can therefore be negative or positive.
		typedef std::map<size_t, int> PeakLocation;
		private:
			// The ID of the isotope cluster.
			int _flag;
			size_t _peakid;
			//size_t _clusterid;
			//int _shift;
			PeakLocation peakinfo;
		public:
			// Construct an isotopepeak.
			IsotopePeak(const Peak& pk, const size_t& peakid)
				:  Peak(pk), _flag(0), _peakid(peakid), peakinfo()
			{
			}
			virtual ~IsotopePeak() {}
			// Add cluster information.
			inline void addClusterInformation(const size_t& cid, const int& shift)
			{
				// See if the cluster id information has been recorded. If yes, just update the inforamtion
				std::map<size_t, int>::iterator iter = peakinfo.find(cid);
				if(iter != peakinfo.end())
					iter->second = shift;
				else
					peakinfo.insert(std::make_pair(cid, shift));
			}

			// If a cluster has been identified not as a candidate isotope cluster, just clean it.
			inline void cleanClusterInformation(const size_t& cid)
			{
				// Erase by key.
				peakinfo.erase(cid);
			}
			// Decide if a peak is included in a cluster.	
			bool includeCluster(const PseudoclusterList& clusterlist);
			
			inline void changeFlag()
			{
				_flag = _flag * (-1);
			}
			inline void setFlag(const size_t& flag_state)
			{
				_flag = flag_state;
			}
			// Search the shift in a given cluster.
			inline int getShift(size_t clusterid) const
			{
				return peakinfo.find(clusterid)->second;
			}
			inline size_t getPeakID() const
			{
				return _peakid;
			}

	};
	
	
	typedef std::vector<IsotopePeak> IsopeakList;
	
	//class Pseudocluster : public SpectrumType<IsotopePeak>
	class Pseudocluster
	{
		private:
			IsopeakList _peaks;
			size_t _clusterid;
			size_t _headpeakid;
			int max_peaks;
			int charge_state;
			int charge_type;
		public:
			Pseudocluster(const size_t& cid, const size_t& hpid, const int& mp, const int& chs, const int& cht = -1)
				: _peaks(), _clusterid(cid), _headpeakid(hpid), max_peaks(mp), charge_state(chs), charge_type(cht)
			{
			}
			Pseudocluster(IsopeakList& peaks)
				: _peaks(peaks), _clusterid(0), _headpeakid(0), max_peaks(0), charge_state(0), charge_type(-1) 
			{
			}
			Pseudocluster(IsopeakList& peaks, const size_t& cid, const size_t& hpid, const int& mp, const int& chs, const int& cht=-1)
				: _peaks(peaks), _clusterid(cid), _headpeakid(hpid), max_peaks(mp), charge_state(chs), charge_type(cht) 
			{
			}
			
			inline int getChargeState() const
			{
				return charge_state;
			}
			
			inline void setChargeState(const int chg)
			{
				charge_state = chg;
			}
			
			inline int getChargeType() const 
			{
				return charge_type;
			}
			
			inline void setChargeType(const int pn)
			{
				charge_type = pn;
			}
			
			inline size_t getHeadID() const
			{
				return _headpeakid;
			}

			inline int getMaxPeaks() const
			{
				return max_peaks;
			}

			bool isIncluded(const IsotopePeak ip) const;
			
			inline size_t getClusterID() const
			{
				return _clusterid;
			}
			
			inline void setClusterID(size_t cid)
			{
				_clusterid = cid;
			}
			
			inline void addPeak(const IsotopePeak& peak)
			{
				_peaks.push_back(peak);
			}
			
			inline void addPeakList(const IsopeakList& peaks)
			{
				_peaks.insert(_peaks.end(), peaks.begin(), peaks.end());
			}
			
			inline IsopeakList getPeakList() const
			{
				return _peaks;
			}

	};
	
	// Personally there will be a PseudoclusterList for each charge.
	class PseudoclusterFinder {
		private:
			Param& _param;
			// Peak list
			//Pseudoclusterlist& _clustergroup;
		public:
			// The real algorithm part.
			PseudoclusterFinder()
				: _param(Param::instance())
			{
			}
			void run(const std::vector<Peak>&, std::map<std::string, PseudoclusterList>&);
			//void setParameters(Param& );
			// inline ClusterList& getCluster() const {
				// return _clustergroup;
			//}
	};
}

#endif // GAG_PSEUDOCLUSTERFINDER_H

