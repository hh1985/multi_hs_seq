#include <GAGPL/SPECTRUM/PseudoclusterFinder.h>
#include <algorithm>
#include <math.h>
#include <iostream>

namespace gag 
{
	bool IsotopePeak::includeCluster(const PseudoclusterList& clusterlist)
	{
		// The situation of empty vector should be included.
		for(size_t idx=0; idx<clusterlist.size(); idx++)
		{
			PeakLocation::iterator it = peakinfo.find(clusterlist[idx].getClusterID());
			if(it!=peakinfo.end())
			{
				return true;
			}
		}
		return false;
	}

	/* The critical function which realize the algorithm */
	void PseudoclusterFinder::run(const std::vector<Peak>& plist, std::map<std::string, PseudoclusterList>& clusterlist_by_charge) 
	{
		/* Step 1: Divide the spectrum into small islands */
		
		// isolist stores the candidate isopeak list.
		IsopeakList isolist;
		// islands stores the list of candidate isopeak list.
		std::vector<IsopeakList> islands;
		

		for(size_t i=0; i<plist.size(); i++)
		{
			// construct an isopeak with given id.
			IsotopePeak ipeak(plist[i], i);
			// Add the isotopepeak into spectrum.
			isolist.push_back(ipeak);

			// End or gap found! The largest mass difference is limited. The distance is manually decided.
			double gap_distance = _param.getParameter<double>("gap_distance").first;
			
			if((i==plist.size()-1) || (plist[i+1].getMZ() - plist[i].getMZ()) > gap_distance)
			{
				// The number of isopeaks are not enough to form a isotope cluster.
				// Simply discard the cluster.
				if(isolist.size() <= 3)
					continue;
				islands.push_back(isolist);
				// Empty the vector and prepare for next cluster.
				isolist.clear();
			}
		}

		//int upper_charge = boost::any_cast<int>(_param.getValue("max_charge"));
		//int lower_charge = boost::any_cast<int>(_param.getValue("min_charge"));
		
		// :) TBD: consider the reasonable range of the upper_charge and lower_charge.
		int upper_charge = _param.getParameter<int>("max_charge").first;
		int lower_charge = _param.getParameter<int>("min_charge").first;
		
		//std::map<std::string, std::vector<IsopeakList> > islands_chg;
		/* Step 2: Iterate over each island, and generate candidate isotope clusters. */
		size_t cluster_index = 0;
		for(size_t k=0; k<islands.size(); k++)
		{	
			// Sort peaks based on their intensity values in a descending order.
			std::sort(islands[k].begin(), islands[k].end(), Peak::intensityLarger);
				
			// Go throught the most intensed peaks to deduce the charge, and search against the rest of the peaks.
			// pk_idx indicates the head peak.
			
			for(size_t pk_idx=0; pk_idx<islands[k].size()-3; pk_idx++)
			{
				
				// Is it possible that the two charge states have different charge types?
				for(int chg=upper_charge; chg>=lower_charge && chg!=0; chg--)
				{
					cluster_index++;	
					int shift = 0;
					// If the head peak has been included in any previous cluster with the same charge, 
					// skip this round.
					
					if(islands[k][pk_idx].includeCluster(clusterlist_by_charge[boost::lexical_cast<std::string>(chg)]))
						continue;
					
					// Update the isopeak information.
					// :)TBD
					
					// Add cluster id, head id and charge into cluster_chg.
					IsopeakList peaks_tmp;
					
					islands[k][pk_idx].addClusterInformation(cluster_index, shift);
					islands[k][pk_idx].setFlag(USED);
					
					/* Add candidate head peak. */
					peaks_tmp.push_back(islands[k][pk_idx]);
					
					//double ppm = boost::any_cast<double>(_param.getValue("tolerance"));
					double ppm0 = 1e-5;
					double ppm1 = _param.getParameter<double>("tolerance").first;
					double ppm2 = _param.getParameter<double>("A2tolerance").first;
		
					for(size_t fl_idx=pk_idx+1; fl_idx<islands[k].size(); fl_idx++)
					{
							double diff = islands[k][fl_idx].getMZ() - islands[k][pk_idx].getMZ();
							//shift= (int)(fabs(diff) + 0.5);
							// The shift value can be either positive or negative.
							shift = (int)floor(diff * (double)chg + 0.5);

							if(shift > 5)
								continue;
							double err = 0.0;
							if(shift == 1) // A+1 peak.
								err = (diff - 1.0035/(double)chg)/ppm1;
							else if(shift == 2) // A+2 peak.
								err = (diff - 1.999807169/(double)chg)/ppm2;
							else 
								err = (diff - shift / (double)chg)/ppm0;


							// Candidate isotope peak with fixed charge.
							if(fabs(err) > 2 * islands[k][pk_idx].getMZ())
								continue;
							
							//std::cout << "ID: " << islands[k][fl_idx].getPeakID() << " Peak: " << islands[k][fl_idx].getMZ() << " Shift: " << shift << std::endl;
							islands[k][fl_idx].addClusterInformation(cluster_index, shift);
							islands[k][fl_idx].setFlag(USED);
							peaks_tmp.push_back(islands[k][fl_idx]);
					}
					
					//std::cout << "Head peak: " << islands[k][pk_idx].getMZ() << " Charge: " << chg 
						//<< " Number of peaks: " << peaks_tmp.size() << std::endl; 
					/* If peaks_tmp meets the basic requirement of candidate isocluster. */
					//int min_peaks = boost::any_cast<int>(_param.getValue("min_peaks"));
					
					
					/* resort the peak list based on their m/z values in an ascending order. */
					std::sort(peaks_tmp.begin(), peaks_tmp.end(), Peak::mzSmaller);
					
					//for(size_t i =0; i< peaks_tmp.size(); i++)
					//{
					//	std::cout << peaks_tmp[i].getMZ() << " " << std::endl;
					//}
					
					/* Calculate maximum number of continuous peak */
					int current_shift = peaks_tmp[0].getShift(cluster_index);
					int current_head = peaks_tmp[0].getPeakID();
					double current_mz = peaks_tmp[0].getMZ();

					size_t current_count = 1;
					size_t max_pks = 1;

					size_t min_num_peaks = _param.getParameter<int>("min_num_peaks").first;
					if(peaks_tmp.size() < min_num_peaks)
					{
						peaks_tmp.clear();

						//for(size_t fl_idx=pk_idx+1; fl_idx<islands[k].size(); fl_idx++)
						//{
						//	islands[k][fl_idx].cleanClusterInformation(cluster_index);
						//}
						continue;
					}

					for(size_t mx=1; mx<peaks_tmp.size(); mx++)
					{
						//std::cout << "Cluster: " << cluster_index << " Shift: " << peaks_tmp[mx].getShift(cluster_index) << std::endl;
						if(current_shift >= 0 && (peaks_tmp[mx].getShift(cluster_index) - current_shift) == 1)
						{
							// Extend the number.	
							current_count += 1;

							// Update max_num;
							if(current_count > max_pks)
								max_pks = current_count;
						} else if((peaks_tmp[mx].getShift(cluster_index) - current_shift) > 1)
						{	
							// Reset the current_head and current_count
							//std::cout << current_mz << " " << current_count << std::endl;
							current_head = peaks_tmp[mx].getPeakID();
							current_count = 1;
							current_mz = peaks_tmp[mx].getMZ();
						}
						// Update current_shift;
						current_shift = peaks_tmp[mx].getShift(cluster_index);
					}
					
					//std::cout << "Max peaks: " << max_pks << std::endl;

					
					size_t max_continuous_peaks = _param.getParameter<int>("max_continuous_peaks").first;

					if(max_pks >=max_continuous_peaks )
					{

						/* Add cluster information */
						Pseudocluster pcluster(peaks_tmp,cluster_index, pk_idx, max_pks, chg);

						clusterlist_by_charge[boost::lexical_cast<std::string>(chg)].push_back(pcluster);
						
						std::cout.precision(10);
						//std::cout << "Found one cluster: " << std::endl;
						//std::cout << cluster_index << "\t" << pk_idx << "\t" << chg << "\t" << max_pks << "\t" << islands[k][pk_idx].getMZ() << std::endl;
						/*for(size_t i =0; i< peaks_tmp.size(); i++)
						{
							std::cout << peaks_tmp[i].getMZ() << "\t" << peaks_tmp[i].getIntensity() 
								<< "\t" << peaks_tmp[i].getShift(cluster_index) << std::endl;
						}*/
						
						
					}
					peaks_tmp.clear();
				}
				
			}			
		
		}
		
	}

}
