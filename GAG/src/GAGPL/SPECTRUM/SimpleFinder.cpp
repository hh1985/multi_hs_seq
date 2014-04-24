/********************************************************************
	created:	2013/01/10
	created:	10:1:2013   9:15
	filename: 	SimpleFinder.cpp
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	SimpleFinder
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include <algorithm>
#include "GAGPL/SPECTRUM/SimpleFinder.h"
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>

namespace gag
{
	SimpleFinder::SimpleFinder(RichList& spec)
		: spectrum(spec), param(Param::Instance()), pt(PeriodicTable::Instance())
	{
		int precursor_charge = param.getParameter<int>("precursor_charge").first;

		double pre_mz = param.getParameter<double>("precursor_mz").first;
		precursor_mass = calculateMass(pre_mz, precursor_charge);

		/* 
		Guess the range of average formula. Usually we assume the difference between composition estimated from the base peak and the one estimated from monoisotopic peak can be safely ignored.
			No Sulfur (100D)		High Sulfur (100D) 
			C: 3.7238523        1.9961864
			H: 5.4425534        2.9942796
			O: 2.8645018        3.1938983
			N: 0.2864502        0.1996186
			S: 0                0.5988559
		*/
		no_sulfur_ave.insert(std::make_pair("C", 3.7238523));
		no_sulfur_ave.insert(std::make_pair("H", 5.4425534));
		no_sulfur_ave.insert(std::make_pair("O", 2.8645018));
		no_sulfur_ave.insert(std::make_pair("N", 0.2864502));
		no_sulfur_ave.insert(std::make_pair("S", 0.0));

		max_sulfur_ave.insert(std::make_pair("C", 1.9961864));
		max_sulfur_ave.insert(std::make_pair("H", 2.9942796));
		max_sulfur_ave.insert(std::make_pair("O", 3.1938983));
		max_sulfur_ave.insert(std::make_pair("N", 0.1996186));
		max_sulfur_ave.insert(std::make_pair("S", 0.5988559));
		run();
	}
	
	void SimpleFinder::run()
	{
		// There is a need to divide the spectrum into separate islands at this moment.
		RichPeakListBySignalOverNoise& pks_sn = spectrum.getPeakListByType<peak_signal_noise>();

		// The S/N threshold for monoisotopic peak.
		double mono_sn_threshold = param.getParameter<double>("signal_noise").first;

		RichPeakListBySignalOverNoise::iterator end_sn = pks_sn.upper_bound(mono_sn_threshold);		
		RichPeakListBySignalOverNoise::iterator iter_sn = pks_sn.begin();
		//end_iter--;

		// Start from the peak with highest s/n
		for(; iter_sn != end_sn; iter_sn++) {
			this->findNextEnvelop(spectrum, *iter_sn);
		}

		std::cout << "Before optimizing, output all the envelop information..." << std::endl;
		BOOST_FOREACH(EnvelopPtr env, env_pool)
			this->printEnvelop(env);

		std::cout << std::endl;

		//std::cout << "After optimizing..." << std::endl;
		//Estimate parameters for each envelop and update the peak information.
		this->optimizeEnvelopSet(env_pool);
		this->optimizeSuspiciousPeaks();

	}

	void SimpleFinder::findNextEnvelop(RichList& pk_list, RichPeakPtr base_pk, bool raw_spec /* = true */)
	{
		// Remove the neighbor peak.
		std::cout << "Check if peak " << base_pk->mz << " is in the black list" << std::endl;
	
		if(!base_pk->pk_status) {
			std::cout << "Peak " << base_pk->mz << " is a noise peak." << std::endl;
			return;
		}

		double mz = base_pk->mz;

		// Process the apodization problem.
		this->cleanNeighborNoise(pk_list, base_pk);

		int precursor_charge = param.getParameter<int>("precursor_charge").first;
		
		double signal_noise = param.getParameter<double>("signal_noise").first;

		int sign = (precursor_charge > 0 ? 1 : -1);
		bool suspicious_pk = true;

		for(int z = abs(precursor_charge); z>=1; z--) {
	
			// The actual mass is allowed to be larger than precursor mass.
			//double mass = calculateMass((*iter_intensity)->mz, sign * z);
			//if(mass > precursor_mass) continue;

			// Not a qualified peak due to its poor intensity.
			if(base_pk->signal_noise < signal_noise) continue;
			
			std::cout << "Identify base peak: " << base_pk->mz << " Charge: " << z;
			std::cout << (precursor_charge > 0 ? "+" : "-") << std::endl;

			bool status = this->extendEnvelop(pk_list, base_pk, z, raw_spec);
			
			if(status) { 
				suspicious_pk = false;
				// If the S/N of the base peak is higher than 200.
				if(base_pk->signal_noise > 200)
					this->detectHarmonicCluster(base_pk, sign * z);
			}

		}

		if(suspicious_pk)
			setPeakType(base_pk, "ISO");
	}

	void SimpleFinder::optimizeEnvelopSet(std::set<EnvelopPtr>& env_set)
	{
		// Start from low m/z to high m/z.
		RichList pk_list = env_ref.getRichPeakList(env_set);
		RichPeakListByMZ& pk_mz = pk_list.getPeakListByType<peak_mz>();

		// Guess the number of sulfate.
		BOOST_FOREACH(EnvelopPtr env, env_set)
			this->updateEnvelopParameter(env);

		for(RichPeakListByMZ::iterator mz_iter = pk_mz.begin(); 
			mz_iter != pk_mz.end(); mz_iter++)
		{
			RichPeakPtr pk = *mz_iter;
			double mz = pk->mz;
			std::cout << "Check peak " << pk->mz << std::endl;

			// 1. For each peak, check if there is any new envelop occurred at this peak.
			std::vector<EnvEntry> new_env_entries = env_ref.getOccurredEnvelopEntries(pk, NEW);
			if(new_env_entries.size() == 0) 
				continue;

			// Update the fitting score for all the new envelop set.
			BOOST_FOREACH(EnvEntry& env_entry, new_env_entries) {
				this->calculateFittingScore(env_entry.env);
			}

			// Check if this is a peak from previously identified envelop.
			std::vector<EnvEntry> old_env_entries = env_ref.getOccurredEnvelopEntries(pk, OLD);
				
			// 1. Sort the envelop vector by their fitting scores.
			// Since the value of the fitting score will be frequently updated,				// it might not be a good idea to keep it in a multi-index container
			std::sort(new_env_entries.begin(), new_env_entries.end(), EnvEntry::scoreLarger);

			// Check if the new envelops are actually just subset of the old 
			// envelops.
			if(old_env_entries.size() != 0) {
				// a. Declared as abnormal by old_env_entry.
				// b. Identified new peak.
				// c. Meet the threshold of the envelop.
				BOOST_FOREACH(EnvEntry& new_entry, new_env_entries)
				{
						bool fake_envelop = false;
						BOOST_FOREACH(EnvEntry old_entry, old_env_entries)
						{
							bool false_positive = this->isRedundantEnvelop(new_entry.env, old_entry.env, old_entry.getShift());
							
							// If this peak is considered by any of the old envelops as								// a true one, move to the next step. 
							if(false_positive) {
								fake_envelop = true;
								break;
							}
						}

						if(fake_envelop) new_entry.info->entry_status = FALSE;
				}
			}

			/* Rule of thumb */
			// 1. If there is only one new envelop, just accept the new envelop.
			if(new_env_entries.size()==1) {
				// No longer new envelop.
				if(new_env_entries[0].info->entry_status != FALSE) {
					// Identify if this is a noise envelop.
					if(!isNoiseEnvelop(new_env_entries[0].env))
						this->acceptEnvelop(new_env_entries[0].env);
					
				}
				continue;
			}

			// 2. If there is more than one new envelop, which means multiple charge states, select the one with higher fitting score as the basis, sequentially merge the other into it, and see if the overall fitting score will be improved.

			// Upgrade the score of the system.
			// Calculate the previous system score.
			double system_score = 0.0;
			
			// 2. Sequentially merge one into the basis and decide if the performance is improved.
			std::vector<EnvEntry> added_entries;
			BOOST_FOREACH(EnvEntry& new_entry, new_env_entries)
			{
				if(new_entry.info->entry_status == FALSE)
					continue;

				// Assign the fitting score to the system score.
				if(system_score == 0.0) {
					system_score = new_entry.env->fitting_score;
					this->acceptEnvelop(new_entry.env);

					added_entries.push_back(new_entry);
					continue;
				}

				// 2.a Find the corresponding fitting score. Only considering the situation where there are multiple charge states for the same mono-isotopic peak. 
				double new_score = updateEnvelopTree(old_env_entries, added_entries, new_entry, pk);	

				// Decide if the envelop is proper to be added into the 
				// envelop set.
				if(new_score <= system_score) {
					new_entry.info->entry_status = FALSE;
				}else {
					/* Accept the candidate envelop.*/
					// No longer new envelop.
					this->acceptEnvelop(new_entry.env);
					added_entries.push_back(new_entry);
					system_score = new_score;
				}
			}

		}
	}

	bool SimpleFinder::extendEnvelop( RichList& pk_list, RichPeakPtr base_pk, int charge, bool raw_spec)
	{
		// a. Using A+1 and A+2 peak to identify the most likely composition.
		// b. Try to extend the peak using the estimated composition.
		// c. Notice that sometimes, the A+2 peak might have missing information due to the defective peak picking algorithm.
		// d. This algorithm requires that A peak has to present.
		
		/* A + 1 peak */
		// Assume it follows the no sulfate model, get the boundary.
		// Guess the number of sulfate. Notice that there might be overlapping.
		// Calculate the no sulfate model. Notice that the information has been corrected by charge.

		AggregatedIsotopicVariants higher_theo = this->estimateDistribution(base_pk->mz, charge);
		//AggregatedIsotopicVariants lower_theo = this->estimateDistribution(base_pk->mz, charge, ave2);

		// Get the isotope information.
		const Isotope s32 = pt.getIsotopeByNominalMass("S", 32);
		const Isotope s34 = pt.getIsotopeByNominalMass("S", 34);
		const Isotope c12 = pt.getIsotopeByNominalMass("C", 12);
		const Isotope c13 = pt.getIsotopeByNominalMass("C", 13);
		const Isotope h = pt.getIsotopeByNominalMass("H",1);
		double electron_mass = param.getParameter<double>("electron_mass").first;
		double low_sn_threshold = param.getParameter<double>("lower_bound_sn").first;

		// The difference between mass of H and peak(A -> A+1). This variable is for re-use convinience.
		double mz_diff = (h.mass + electron_mass)/(double)charge - higher_theo.getMassDifferenceByShift<peak_intensity>(0, 1);

		// Temporary variable for storing the peak information.
		// 1. Shift;
		// 2. Peak;
		// 3. Status.
		std::multimap<int, std::pair<RichPeakPtr, std::string> > pk_map;
		pk_map.insert(std::make_pair(0, std::make_pair(base_pk, "Normal")));

		int shift = 1;

		while(1) {
			// Theoretically, the (A+n)' should always be present.
			double mz0 = base_pk->mz + higher_theo.getMassDifferenceByShift<peak_intensity>(0,shift);

			std::string status = "Normal";
			if(shift == 1) {			
				std::set<RichPeakPtr> pks_set = this->getClosestRichPeaks(pk_list, mz0, raw_spec, low_sn_threshold);
				
				if(pks_set.size() == 0) { // C+H peak for different shift.
									
					double mz1 = mz0 + mz_diff;
					
					pks_set = this->getClosestRichPeaks(pk_list, mz1, raw_spec, low_sn_threshold);
					status = "C+H";

					if(pks_set.size() == 0 && shift > 1) { // C+2H
						double mz2 = mz0 + 2.0 * mz_diff;
						pks_set = this->getClosestRichPeaks(pk_list, mz2, raw_spec, low_sn_threshold);
						status = "C+2H";
						
						if(pks_set.size() == 0)
							break;
					} else if(pks_set.size() == 0 && shift == 1) {
						break;
					}
					
				}
				
				BOOST_FOREACH(RichPeakPtr pk, pks_set) {
					pk_map.insert(std::make_pair(shift, std::make_pair(pk, status)));
				}

			} else if(shift == 2) {
				// If 34S peak is present.
				double mz1 = base_pk->mz + (s34.mass - s32.mass)/(double)charge;
				std::set<RichPeakPtr> pks_set1 = this->getClosestRichPeaks(pk_list, mz1, raw_spec, low_sn_threshold);

				BOOST_FOREACH(RichPeakPtr pk, pks_set1) {
					pk_map.insert(std::make_pair(shift, std::make_pair(pk, "Sulfur")));
				}

				// (A+n)' peak.
				std::set<RichPeakPtr> pks_set2 = this->getClosestRichPeaks(pk_list, mz0, raw_spec, low_sn_threshold);

				BOOST_FOREACH(RichPeakPtr pk, pks_set2) {
					pk_map.insert(std::make_pair(shift, std::make_pair(pk, "Normal")));
				}		

				if(pks_set1.size() == 0 && pks_set2.size() == 0) {
					// C+H peak.
					double mz2 = mz0 + mz_diff;

					std::set<RichPeakPtr> pks_set3 = this->getClosestRichPeaks(pk_list, mz2, raw_spec, low_sn_threshold);
					BOOST_FOREACH(RichPeakPtr pk, pks_set3) {
						pk_map.insert(std::make_pair(shift, std::make_pair(pk, "C+H")));
					}

					if(pks_set3.size() == 0) {
						// C + 2H
						double mz3 = mz0 + 2.0 * mz_diff;
						std::set<RichPeakPtr> pks_set4 = this->getClosestRichPeaks(pk_list, mz3, raw_spec, low_sn_threshold);
						BOOST_FOREACH(RichPeakPtr pk, pks_set4) {
							pk_map.insert(std::make_pair(shift, std::make_pair(pk, "C+2H")));
						}
						if(pks_set4.size() == 0) 
							break;
						
					}
				}		

			} else { // Shift >= 2.
				// Left boundary.
				int sulfur_multiplier = shift / 2;
				size_t count = 0;
				for(int i = sulfur_multiplier; i>=0; i--) {
					int carbon_multiplier = shift - i * 2;
					double mzx = base_pk->mz + i * (s34.mass - s32.mass)/(double)charge + carbon_multiplier * (c13.mass - c12.mass)/(double)charge;
					std::set<RichPeakPtr> pks_set = this->getClosestRichPeaks(pk_list, mzx, raw_spec, low_sn_threshold);
					count += pks_set.size();
					BOOST_FOREACH(RichPeakPtr pk, pks_set) {
						pk_map.insert(std::make_pair(shift, std::make_pair(pk, "C+S")));
					}
          
          // Also consider the H shift.
          double mz1 = mzx + mz_diff;
          pks_set = this->getClosestRichPeaks(pk_list, mz1, raw_spec, low_sn_threshold);
          if(pks_set.size() != 0) {
            status = "C+H";
            BOOST_FOREACH(RichPeakPtr pk, pks_set) {
              pk_map.insert(std::make_pair(shift, std::make_pair(pk, status)));
            }
            //pk_map.insert(std::make_pair(shift, std::make_pair(pk, "C+S")));
          } 
				}
				// No peak has been found in this range.
				if(count == 0) {
          break;
        }
					
			}

			shift++;
		}

		// No additional peak has been found.
		if(pk_map.size() == 1)
			return false;
		
		EnvelopPtr env_ptr = createEnvelop(charge);

		double param_confidence = 0.0;

		std::multimap<int, std::pair<RichPeakPtr, std::string> >::iterator map_iter = pk_map.begin();
		for(; map_iter != pk_map.end(); map_iter++)
		{
			int pk_pos = map_iter->first; 
			RichPeakPtr pk_ptr = map_iter->second.first;
			this->cleanNeighborNoise(pk_list, pk_ptr);
			std::string pk_status = map_iter->second.second;

			EntryStatus entry_status = (pk_pos == 0 ? NEW : OLD);
			
			InfoPeakPtr pk_infor = boost::make_shared<InfoPeak>(pk_ptr->intensity, pk_status, entry_status);

			env_ref.addDictionaryReference(pk_ptr, env_ptr, pk_pos, pk_infor);
		}

		int max_shift = shift - 1; 
		env_ptr->max_shift = max_shift;

		std::cout << "Found candidate envelop -- mono peak: " << base_pk->mz << " charge: " << env_ptr->charge_state << std::endl;
		env_pool.insert(env_ptr);

		//calculateFittingScore(env_ptr);

		return true;

	}

	AggregatedIsotopicVariants SimpleFinder::estimateDistribution( double mz, int charge, bool sulfur)
	{
		double mass = calculateMass(mz, charge);
		double coef = mass / 100.0;
		Composition compo;

		if(sulfur) {
			// Correct sulfur number. 
			AveragineFormulae::iterator sulfur_iter = max_sulfur_ave.find("S");
			int num = (int)ceil(sulfur_iter->second * coef);

			double new_mass = mass - (double)num * Composition("SO3").getMass();
			
			double new_coef = new_mass / 100.0;
			
			AveragineFormulae::iterator iter = max_sulfur_ave.begin();
			for(; iter != max_sulfur_ave.end(); iter++) {
				if(iter->first != "S") {
					compo.addElement(iter->first, (int)floor(iter->second * new_coef + 0.5));
				} else {
					compo.addElement(iter->first, num);
				}
			}
			IsotopicDistribution iso(compo);
			return iso.getAggregatedIsotopicVariants(charge);
		} else {
			AveragineFormulae::iterator iter = no_sulfur_ave.begin();
			for(; iter != no_sulfur_ave.end(); iter++)
				compo.addElement(iter->first, (int)floor(iter->second * coef + 0.5));

			IsotopicDistribution iso(compo);
			return iso.getAggregatedIsotopicVariants(charge);
		}


	}

	AggregatedIsotopicVariants SimpleFinder::estimateDistribution( double mz, int charge, int sulfur_num )
	{
		double mass = calculateMass(mz, charge);
		double new_mass = mass - (double)sulfur_num * Composition("SO3").getMass();

		double new_coef = new_mass / 100.0;

		AveragineFormulae::iterator iter = no_sulfur_ave.begin();
		Composition compo;
		for(; iter != no_sulfur_ave.end(); iter++) {
			if(iter->first == "S") {
				compo.addElement(iter->first, sulfur_num);
			} else if(iter->first == "O") {
				compo.addElement(iter->first, sulfur_num * 3);
			} else {
				compo.addElement(iter->first, (int)floor(iter->second * new_coef + 0.5));
			}
		}
		IsotopicDistribution iso(compo);
    //AggregatedIsotopicVariants agr_var = iso.getAggregatedIsotopicVariants(charge);
		return iso.getAggregatedIsotopicVariants(charge);
	}


	std::set<RichPeakPtr> SimpleFinder::getClosestRichPeaks(RichList& pk_list, double expected_mz, bool raw_spec, double sn_threshold, double error)
	{
		RichPeakListByMZ& pks_mz = pk_list.getPeakListByType<peak_mz>();

		RichPeakListByMZ::iterator mz_increase_iter, mz_decrease_iter;
		mz_increase_iter = pks_mz.lower_bound(expected_mz);
		mz_decrease_iter = mz_increase_iter;

		std::set<RichPeakPtr> pk_set;
		std::map<double, RichPeakPtr> pk_map;

		// Examine the peak towards two directions.
		while(1) {
			if(mz_increase_iter == pks_mz.end())
				break;

			RichPeakPtr pk = *mz_increase_iter;

			bool meet_threshold = raw_spec ? true : (pk->signal_noise > sn_threshold); 

			double win_error = this->getWindowError(pk, error);
			
			if(abs(pk->mz - expected_mz) < win_error) {
				if(meet_threshold)
					pk_map.insert(std::make_pair(pk->intensity, pk));
			} else {
				break;
			}

			mz_increase_iter++;
			if(mz_increase_iter == pks_mz.end())
				break;
		}

		mz_decrease_iter--;
		while(1) {
			RichPeakPtr pk = *mz_decrease_iter;

			double win_error = this->getWindowError(pk, error);
			
			bool meet_threshold = raw_spec ? true : (pk->signal_noise > sn_threshold);

			if(abs(pk->mz - expected_mz) < win_error)	{		
				if(meet_threshold)
					pk_map.insert(std::make_pair(pk->intensity, pk));
			} else {
				break;
			}

			if(mz_decrease_iter == pks_mz.begin())
				break;
			else
				mz_decrease_iter--;
		}

		//pk_set.insert(selected_pk);
		if(pk_map.size()>0) {
			RichPeakPtr selected_pk = pk_map.rbegin()->second;
			pk_set.insert(selected_pk);
		}
		
		return pk_set;
	}


	void SimpleFinder::updateEnvelopParameter( EnvelopPtr env )

	{
		// 0. Check env id.
		unsigned int env_id = env->id;
		std::cout << "Check env #" << env_id << std::endl;

		// 1. Get all the peaks.
		EnvEntry base_entry = env_ref.getEntryByShift(env, 0);
		
		std::vector<EnvEntry> entry_vec = env_ref.getEntryByEnvelop(env);

		AggregatedIsotopicVariants higher_theo = this->estimateDistribution(base_entry.pk->mz, env->charge_state);
		AggregatedIsotopicVariants lower_theo = this->estimateDistribution(base_entry.pk->mz, env->charge_state, true);

		int max_sulfur = this->guessMaxSulfurNumber(base_entry.pk->mz, env->charge_state);
		int max_num = param.getParameter<int>("max_sulfur_num").first;

		const Isotope s34 = pt.getIsotopeByNominalMass("S", 34);
	 
		std::pair<int, int> sulfur_range = std::make_pair(0,max_sulfur);

		// Gradually estimate the sulfur number. The default value is 0.
		BOOST_FOREACH(EnvEntry entry, entry_vec)
		{
			// For given shift, estimate the range of the intensities.
			int pk_pos = entry.getShift();
			
			if(pk_pos == 0) continue;

			// Expected abundance.
			double non_sulfur_model_abundance = base_entry.pk->intensity * higher_theo.getPeakByShift<peak_intensity>(pk_pos)->intensity;
			double sulfur_model_abundance = base_entry.pk->intensity * lower_theo.getPeakByShift<peak_intensity>(pk_pos)->intensity;

			double max_abundance = non_sulfur_model_abundance;
			double min_abundance = sulfur_model_abundance;
			if(min_abundance > max_abundance)
				swap(max_abundance, min_abundance);

			// Narrow down the range of sulfur_range.
			// For tandem ms analysis, only use A+1 and A+2 peak.
			// The information from A+1 peak might help increase the lower limit 
			if(pk_pos == 1) { // odd peak, increase the lower limit of sulfur number.
				
				if(entry.info->status != "Normal") continue;

				if((entry.info->adjusted_abundance-max_abundance)/max_abundance >0.2 ) {
					// Significantly higher than upper boundary.
					// Nothing to do with the range.
				} else if((entry.info->adjusted_abundance-min_abundance)/min_abundance < -0.2) {
					// Significantly lower than expected.
					// env->suspicious = true;
				} else {
					// Normal situation. Estimate sulfur_num from the abundance.
					double lambda = (entry.info->adjusted_abundance-min_abundance)/(max_abundance - min_abundance);

					// Convert lambda to sulfur_num.
					if(sulfur_model_abundance == min_abundance)
						lambda = 1.0 - lambda;

					int sulfur_num = (int)floor(lambda * max_sulfur + 0.5);

					//entry.env->sulfur_num = sulfur_num > max_num ? max_num : sulfur_num;
					if(sulfur_num < sulfur_range.second)
						sulfur_range.second = sulfur_num;
					
				}
			} else if(pk_pos == 2) { // even peak, decrease the upper limit of sulfur number. Notice that one of the peak might be missing.
				if(entry.info->status == "Sulfur") {
					// TBD: decide if the sulfur peak and non-sulfur peak will separate.
					// If present, estimate sulfur using this information.
					int sulfur_num = (int)floor(entry.info->adjusted_abundance/base_entry.info->adjusted_abundance / pow(s34.abundance, pk_pos/2) + 0.5);
				
					if(sulfur_num > sulfur_range.first && sulfur_num < max_num)
						sulfur_range.first = sulfur_num;
					//entry.env->sulfur_num = sulfur_num > max_num ? max_num : sulfur_num;
					// The information from sulfur peak have enough confidence.
					break;

				} else if(entry.info->status == "Normal") {
					// Assume the peak and 34S peak are able to be resolved.
					// Define the intensity range for current peak.


					// Expected abundance of the non-sulfur-peak.
					double non_sulfur_peak_abundance = sulfur_model_abundance - base_entry.info->adjusted_abundance * pow(s34.abundance, pk_pos);
					// Do the same thing for A+1 peak.
					max_abundance = non_sulfur_model_abundance;
					min_abundance = non_sulfur_peak_abundance;

					if(max_abundance < min_abundance) 
						swap(max_abundance, min_abundance);
					
					if((entry.info->adjusted_abundance - max_abundance)/max_abundance > 0.2) {
						// Significantly higher than expected.
						// Do nothing.
					} else if((entry.info->adjusted_abundance - min_abundance)/min_abundance < -0.2) {
						// Significantly lower than expected.
						// Do nothing.
					} else {

						double lambda = (entry.info->adjusted_abundance-min_abundance)/(max_abundance - min_abundance);

						// Sometimes, the distance between max_abundance and min_abundance is too small, which will cause a lot of troubles.
						if(lambda < 0.0) continue;

						if(non_sulfur_peak_abundance == min_abundance)
							lambda = 1.0 - lambda;

						// Convert lambda to sulfur_num.
						int sulfur_num = (int)floor(lambda * max_sulfur + 0.5);
						if(sulfur_num < sulfur_range.second)
							sulfur_range.second = sulfur_num;
						//entry.env->sulfur_num = sulfur_num > max_num ? max_num : sulfur_num;
					}

				} else {
					// C+H peak. Nothing to do with the sulfur number.
				}
			} 
		}
		
		// Get the sulfur number estimation by averaing all possible 
		env->sulfur_num = (sulfur_range.first + sulfur_range.second)/2;

		// Unlikely situation.
		if(sulfur_range.first > sulfur_range.second) 
			env->env_status = FP;

		// Update the theoretical distribution for current envelop.
		env->theo_dist = this->estimateDistribution(base_entry.pk->mz, env->charge_state, env->sulfur_num);

		std::map<int, PeakPtr> pk_map = env->getTheoreticalPeaks();

		for(std::map<int, PeakPtr>::iterator iter = pk_map.begin(); 
			iter != pk_map.end(); iter++)
		{
			// Theo peak information.
			int shift = iter->first; PeakPtr theo_pk = iter->second;
			double theo_abundance = base_entry.info->adjusted_abundance * theo_pk->intensity;
			
			// No need to estimate base peak.
			if(shift == 0) continue;

			// Be careful of the situation that the shift for the envelop cannot be observed.
			std::vector<EnvEntry> exp_pks = env_ref.getPeaksByShift(env, shift);

			// Depending on the status of the peak, the theoretical information will be recalculated. It might be possible that some of the peaks are missing, it is better to calculate the mean of the shift for each position.
			// If the sulfur number is 0, we assume there is no peak split.
			if(env->sulfur_num == 0) {
				BOOST_FOREACH(EnvEntry& pk_entry, exp_pks) {
					if(pk_entry.info->status == "Normal") {
						pk_entry.info->relative_shift = (pk_entry.info->adjusted_abundance - theo_abundance) / theo_abundance;
					}
				}
			} else {
				if(shift == 1 || shift == 2) {
					BOOST_FOREACH(EnvEntry& pk_entry, exp_pks) {
						if(shift == 1) {
							if(pk_entry.info->status == "Normal") {
								pk_entry.info->relative_shift = (pk_entry.info->adjusted_abundance - theo_abundance)/ theo_abundance;
							}
						} else if(shift == 2){
							// Here we assume the peak will split.
							if(pk_entry.info->status == "Sulfur") {
								// Calculate the theoretical abundance and update the relative shift.
								theo_abundance = base_entry.info->adjusted_abundance * pow(s34.abundance, env->sulfur_num);
								pk_entry.info->relative_shift = (pk_entry.info->adjusted_abundance - theo_abundance) / theo_abundance;
							} else if(pk_entry.info->status == "Normal"){
								// Notice that theo_abundance should be larger than 0.
								theo_abundance = base_entry.info->adjusted_abundance * (theo_pk->intensity - pow(s34.abundance, env->sulfur_num));
								pk_entry.info->relative_shift = (pk_entry.info->adjusted_abundance - theo_abundance) / theo_abundance;
							} else {
								// C+H peak. Nothing to do with the relative shift.
							}
						} 
					}
				} else {
					// A + n peak
					double total_shift_abundance = 0.0;
					BOOST_FOREACH(EnvEntry& pk_entry, exp_pks) {
						if(pk_entry.info->status == "Normal" 
							|| pk_entry.info->status == "Sulfur") {
								total_shift_abundance += pk_entry.info->adjusted_abundance; 
						}
					}
					BOOST_FOREACH(EnvEntry& pk_entry, exp_pks) {
						pk_entry.info->relative_shift = (pk_entry.info->adjusted_abundance - theo_abundance)/ theo_abundance;
					}
				}
			}

		}


	}

	bool SimpleFinder::isRedundantEnvelop( EnvelopPtr new_env, EnvelopPtr old_env, int old_shift)
	{
		bool state = false;

    std::set<RichPeakPtr> old_pks = env_ref.getPeaksByEnvelop(old_env);
    std::set<RichPeakPtr> new_pks = env_ref.getPeaksByEnvelop(new_env);

    if(std::includes(old_pks.begin(), old_pks.end(), new_pks.begin(), new_pks.end())) {
      // Calculate the sub fitting score of the subset from old_pks.
      if(new_env->fitting_score < old_env->fitting_score)
        state = true;

    } else {
      if(new_env->fitting_score < old_env->fitting_score && (old_env->charge_state % new_env->charge_state == 0))
        state = true;

      EnvEntry base_entry = env_ref.getEntryByShift(old_env, 0);
      //  //EnvEntry new_base = env_ref.getEntryByShift(new_env, 0);
      EnvEntry current_entry = env_ref.getEntryByShift(old_env, old_shift);
      //  // Check if it is different from theoretical abundance.
      AggregatedIsotopicVariants higher_theo = this->estimateDistribution(base_entry.pk->mz, old_env->charge_state);
      AggregatedIsotopicVariants lower_theo = this->estimateDistribution(base_entry.pk->mz, old_env->charge_state, true);

      double non_sulfur_abundance = higher_theo.getPeakByShift<peak_intensity>(old_shift)->intensity * base_entry.info->adjusted_abundance;
      double sulfur_abundance = lower_theo.getPeakByShift<peak_intensity>(old_shift)->intensity * base_entry.info->adjusted_abundance;

      double expected_abundance(0.0);
      if(current_entry.info->status == "Normal") {
        expected_abundance = std::max(non_sulfur_abundance,sulfur_abundance);
      } else {
        if(current_entry.info->status == "Sulfur")
          state = true;

        return state;
      }

      double upper_shift = (current_entry.pk->intensity - expected_abundance)/(expected_abundance);

      double max_intensity_shift = param.getParameter<double>("max_intensity_shift").first;
      
      if(upper_shift < max_intensity_shift)
      	state = true;
    }

		return state;
	}

	void SimpleFinder::printEnvelop( EnvelopPtr env )
	{
		std::cout << "Envelop ID: " << env->id << std::endl;
		std::cout << "Charge: " << env->charge_state << std::endl;
    std::cout << "Fitting score: " << env->fitting_score << std::endl;
		std::cout << std::endl;
		// Peak information.
		std::cout << "Covered peaks:" << std::endl;
		std::cout << "Shift\tMZ\tEXP_ABD\tSTATUS\tSTATUS2" << std::endl;
    
		// Get all peaks.
		std::vector<EnvEntry> pk_entries = env_ref.getEntryByEnvelop(env);
		std::cout.precision(5);
		BOOST_FOREACH(EnvEntry& pk_entry, pk_entries)
		{
			std::cout << pk_entry.getShift() << "\t" ;
			std::cout << std::fixed << pk_entry.pk->mz << "\t";
			std::cout	<< pk_entry.pk->intensity << "\t" << pk_entry.info->status << std::endl;
		}
		std::cout << std::endl;
	}

	std::multimap<RichPeakPtr, EnvelopPtr> SimpleFinder::getEnvelops(EntryStatus status)
	{
		std::vector<EnvEntry> entry_vec = env_ref.getBaseEntriesByStatus();
		
		//EnvDictByEnvID& env_by_eid = env_
		std::multimap<RichPeakPtr, EnvelopPtr> env_map;
		BOOST_FOREACH(EnvEntry& env_entry, entry_vec)
		{
			// Get the base peak for the envelop.
			RichPeakPtr base_pk = env_entry.pk;
			env_map.insert(std::make_pair(base_pk, env_entry.env));
		}
		return env_map;
	}

	double SimpleFinder::getWindowError( RichPeakPtr pk, double error )
	{
		// If error is specified.
		if(error != -1.0) return error * pk->mz;

		double internal_accuracy = param.getParameter<double>("internal_accuracy").first;
		
		//double resolving_power = pk->resolution;
		//double coef = pk->mz / resolving_power;
		double coef = (pk->signal_noise > 15 && pk->mz > 400.0) ? 1.2 : 1.5;
    //double coef = 1.5;
		return 2e-6 * internal_accuracy * pk->mz * coef;
	}

	bool SimpleFinder::isNoiseEnvelop( EnvelopPtr env )
	{	
		unsigned int env_id = env->id;
		// Identification of noise envelop.
		// Situation 1. Weired cluster shape. envelop size has to be larger than 3. This is not noise envelop.
		std::set<int> shift_set = env_ref.getShiftSet(env);
		RichPeakPtr base_pk = env_ref.getBasePeakForEnvelop(env);

		if(shift_set.size() > 2) {
			BOOST_FOREACH(int shift, shift_set)
			{
				std::vector<EnvEntry> pk_entry_vec = env_ref.getPeaksByShift(env, shift);
				BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec)
				{
					if(pk_entry.info->relative_shift < -0.3) {
						env->env_status = MISC;
						// Misc envelop is not noise envelop.
						return false;
					}

				}
			}
		}

		// Situation 2. Matching to noise peak.
		// 2.a size == 2
		if(shift_set.size() == 2) {
			// Get A+1 peak.
			std::vector<EnvEntry> pk_entry_vec = env_ref.getPeaksByShift(env, 1);
			if(pk_entry_vec.size() == 1 && pk_entry_vec[0].info->status == "C+H") {
				// Get theoretical abundance.
				double theo_abundance = env_ref.getTheoreticalPeak(env, 1)->intensity;
				if(theo_abundance - pk_entry_vec[0].info->adjusted_abundance > 0)
					return true;
			}

			RichPeakPtr base_pk = env_ref.getBasePeakForEnvelop(env);
			std::vector<EnvEntry> old_env_entries = env_ref.getOccurredEnvelopEntries(base_pk, OLD);
			// 2.b A + 1 peak is significantly lower than expected.
			if(pk_entry_vec[0].info->relative_shift < -0.2 || old_env_entries.size() > 0)		{
				return true;
			}
		}

		return false;

	}

	void SimpleFinder::cleanNeighborNoise( RichList& pk_list, RichPeakPtr pk/*, std::set<RichPeakPtr>& blacklist*/ )
	{
		RichPeakListByMZ& pks_mz = pk_list.getPeakListByType<peak_mz>();

		RichPeakListByMZ::iterator mz_increase_iter, mz_decrease_iter;
		mz_increase_iter = pks_mz.lower_bound(pk->mz);
		mz_decrease_iter = mz_increase_iter;

		mz_increase_iter++;

		double win_error = 0.012;

		while(1) {
			if(mz_increase_iter == pks_mz.end())
				break;

			RichPeakPtr neighbor_pk = *mz_increase_iter;

			if(abs(neighbor_pk->mz - pk->mz) < win_error && (neighbor_pk->intensity / pk->intensity < 0.1)) {
				//blacklist.insert(neighbor_pk);
				//neighbor_pk->pk_status = false;
				this->setNoisePeak(neighbor_pk);
			} else if(abs(neighbor_pk->mz - pk->mz) >= win_error)
				break;

			mz_increase_iter++;
			if(mz_increase_iter == pks_mz.end())
				break;
		}

		// This is the begin of the peak list.
		if(mz_decrease_iter == pks_mz.begin())
			return;

			mz_decrease_iter--;
		while(1) {
			RichPeakPtr neighbor_pk = *mz_decrease_iter;

			if(abs(neighbor_pk->mz - pk->mz) < win_error && (neighbor_pk->intensity / pk->intensity < 0.1)) {
				
				//neighbor_pk->pk_status = false;
				this->setNoisePeak(neighbor_pk);
			} else if(abs(neighbor_pk->mz - pk->mz) >= win_error)
				break;

			if(mz_decrease_iter == pks_mz.begin())
				break;
			else
				mz_decrease_iter--;
		}
	}

	double SimpleFinder::updateEnvelopTree( std::vector<EnvEntry>& old_entries, std::vector<EnvEntry>& added_entries, EnvEntry& test_entry, RichPeakPtr pk )
	{
		double rest_intensity = 0.0;
		BOOST_FOREACH(EnvEntry& entry, old_entries)
			rest_intensity += entry.info->adjusted_abundance;

		rest_intensity = test_entry.pk->intensity - rest_intensity;

		// Construct the vector using theoretical isotopic distributions. Scale the second one and get the best fit.
		double scaling_factor = 0.0; 
		double scaling_range[2] = {0.0, 1.0};

		double total_score = this->calculateLinearAssociationScore(added_entries);
		
		std::vector<EnvEntry> temp_entries = *(&added_entries);
		temp_entries.push_back(test_entry);

		while(1)
		{
			// Based on scaling factor, update the peak abundance.
			std::vector<EnvEntry> pk_entry_vec = env_ref.getEntryByEnvelop(test_entry.env);
			BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec)
			{
				pk_entry.info->adjusted_abundance *= scaling_factor;
			}

			BOOST_FOREACH(EnvEntry& env_entry, added_entries)
			{
				std::vector<EnvEntry> pk_entry_vec2 = env_ref.getEntryByEnvelop(env_entry.env);
				BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec2)
				{
					pk_entry.info->adjusted_abundance *= (1-scaling_factor);
				}
			}

			// Update scaling factor. It depends on the peaks associated with the second test_entry.
			double new_score = this->calculateLinearAssociationScore(temp_entries);

			double ratio0 = this->calculateScalingRatio(temp_entries);
			double ratio1 = this->calculateScalingRatio(test_entry);
			
			int scale_index = ratio0 > ratio1 ? 1 : 0;

			// Keep a copy of the current scaling factor.
			double old_scaling_factor = scaling_factor;

			// Update the scaling factor of the new envelop.
			scaling_factor = (scaling_range[scale_index] + scaling_factor)/2.0;
			
			// Update the boundary of the scale range.
			scaling_range[1-scale_index] = old_scaling_factor;

			if(abs(new_score - total_score) < 1e-6) 
				break;
			else 
				total_score = new_score;
		}

		double total_fitting_score = 0.0;
		BOOST_FOREACH(EnvEntry& env_entry, temp_entries) {
			// Reset the scaling factor to 1.
			env_entry.env->scaling_factor = 1.0;
			total_fitting_score += this->calculateFittingScore(env_entry.env);
		}
		return total_fitting_score;
		
	}

	double SimpleFinder::calculateFittingScore( EnvelopPtr env )
	{	
		// Get base entry.
		double fitting_score = 0.0;

		std::set<int> shift_set = env_ref.getShiftSet(env);

		// If sulfur peak is present, re-estimate the abundance.

		BOOST_FOREACH(int shift, shift_set)
		{
			double expected_abundance = env_ref.getTheoreticalPeak(env, shift)->intensity;
			double experimental_abundance = env_ref.getExperimentalAbundanceByShift(env, shift);

			// The fitting score comes from 3 parts: 
			//1. the shift of intensity: abs(E-T)/T.
			double intensity_shift = (expected_abundance - experimental_abundance)/ expected_abundance;
			if(intensity_shift > 0.0 && intensity_shift <= 1.0)
				intensity_shift = 1.0 - intensity_shift;
			else if(intensity_shift <= 0.0 && intensity_shift >= -1.0)
				intensity_shift = sqrt(1.0 + intensity_shift);
			else if(intensity_shift > 1)
				intensity_shift = 1.0;
			else
				intensity_shift = 0.0;

			// 2. the shift of m/z. Since this is for high resolution data, the shift of m/z should be ignored.

			// 3. the square root of expected_abundance.
			double weight = sqrt(expected_abundance);

			fitting_score += weight * intensity_shift;
		}

		env->fitting_score = fitting_score;
		return fitting_score;

	}

	double SimpleFinder::calculateFittingScore(std::set<EnvelopPtr>& env_set)
	{
		double total_fitting_score = 0.0;
		BOOST_FOREACH(EnvelopPtr env, env_set)
			total_fitting_score += this->calculateFittingScore(env);
		return total_fitting_score;
	}

	double SimpleFinder::calculateLinearAssociationScore( std::vector<EnvEntry>& entries )
	{
		// Calculate the value of cos theta.
		double inner_prod = 0.0; double norm_term[2] = {0.0,0.0};
		BOOST_FOREACH(EnvEntry& env_entry, entries)
		{
			std::set<int> shift_set = env_ref.getShiftSet(env_entry.env);

			BOOST_FOREACH(int shift, shift_set)
			{
				// Calculate the accumulated abundance for multiple peaks mapped into the same shift position.
				std::vector<EnvEntry> pk_entry_vec = env_ref.getPeaksByShift(env_entry.env, shift);
				double accumulated_abundance = 0.0;
				BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec)
					accumulated_abundance += pk_entry.info->adjusted_abundance;

				double theo_abundance = env_entry.env->getTheoreticalPeak(shift)->intensity;
				inner_prod += theo_abundance * accumulated_abundance;

				norm_term[0] += pow(theo_abundance, 2);
				norm_term[1] += pow(accumulated_abundance, 2);

			}
		}

		return inner_prod / (sqrt(norm_term[0]) * sqrt(norm_term[1]));

	}

	double SimpleFinder::calculateScalingRatio( std::vector<EnvEntry>& entries )
	{
		double adjusted_total = 0.0; double theo_total = 0.0;

		BOOST_FOREACH(EnvEntry& env_entry, entries)
		{
			std::set<int> shift_set = env_ref.getShiftSet(env_entry.env);
			BOOST_FOREACH(int shift, shift_set) {
				// Get the acutal peaks for given position. There might be more than 1 which can be mapped into one peak.
				std::vector<EnvEntry> pk_entry_vec = env_ref.getPeaksByShift(env_entry.env, shift);

				BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec)
					adjusted_total += pk_entry.info->adjusted_abundance;

				// There will be only one theoretical peak for each shift position.
				theo_total += env_entry.env->getTheoreticalPeak(shift)->intensity;
			}
		}

		return adjusted_total / theo_total;
	}

	double SimpleFinder::calculateScalingRatio( EnvEntry& env_entry )
	{
		// Compare the total intensity of the suggested envleop over the theoretical intensity.
		double adjusted_total = 0.0; double theo_total = 0.0;
		
		std::set<int> shift_set = env_ref.getShiftSet(env_entry.env);
		BOOST_FOREACH(int shift, shift_set) {
			// Get the acutal peaks for given position. There might be more 
			// than 1 which can be mapped into one peak.
			std::vector<EnvEntry> pk_entry_vec = env_ref.getPeaksByShift(env_entry.env, shift);

			BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec)
				adjusted_total += pk_entry.info->adjusted_abundance;
			
			// There will be only one theoretical peak for each shift position.
			theo_total += env_entry.env->getTheoreticalPeak(shift)->intensity;
		}

		return adjusted_total/theo_total;

	}

	int SimpleFinder::guessMaxSulfurNumber( double mz, int charge_state )
	{
		double mass = calculateMass(mz, charge_state);
		double coef = mass / 100.0;

		AveragineFormulae::iterator sulfur_iter = max_sulfur_ave.find("S");
		return (int)ceil(sulfur_iter->second * coef);
	}

	void SimpleFinder::acceptEnvelop( EnvelopPtr env )
	{
		// 1. Set the status of the envelop as true positive.
		env->env_status = TP;
		// 2. Set all the peaks associated to the envelop as OLD and the peak type as ENV.
		std::vector<EnvEntry> pk_entry_vec = env_ref.getEntryByEnvelop(env);
		
		
		BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec) {
			pk_entry.info->entry_status = OLD;
			setPeakType(pk_entry.pk, "ENV");
			//pk_entry.pk->pk_type = "ENV";
		}
		
	}

	void SimpleFinder::setNoisePeak( RichPeakPtr pk )
	{
		pk->pk_status = false;
		setPeakType(pk, "NOISE");
	}

	void SimpleFinder::optimizeSuspiciousPeaks()
	{
		// 1. For each suspicious peak, iterate over all charge state, try to decide if there is any peak that can follow into the range.
		RichPeakListBySignalOverNoise& pks_sn = spectrum.getPeakListByType<peak_signal_noise>();

		// The S/N threshold for monoisotopic peawk.
		double mono_sn_threshold = param.getParameter<double>("signal_noise").first;

		RichPeakListBySignalOverNoise::iterator end_sn = pks_sn.upper_bound(mono_sn_threshold);		
		RichPeakListBySignalOverNoise::iterator iter = pks_sn.begin();
		//RichPeakListBySignalOverNoise::iterator iter;


		/*RichPeakListByType& pks_type = spectrum.getPeakListByType<peak_type>();

		std::pair<RichPeakListByType::iterator, RichPeakListByType::iterator> p = pks_type.equal_range("ISO");*/

		int precursor_z = param.getParameter<int>("precursor_charge").first;
		int charge_sign = precursor_z > 0 ? 1 : -1;

		//for(RichPeakListByType::iterator iter = p.first; iter != p.second; iter++)
		for(; iter != end_sn; iter++ )
		{
			RichPeakPtr base_pk = *iter;
			if(base_pk->pk_type != "ISO")
			{
				continue;
			}

			int count = 0;
			for(int z = 1; z < abs(precursor_z); z++)
			{
				// Consider only A+1 peak.
				std::pair<PeakPtr, PeakPtr> pk_pair = getPeakRangeByShift(base_pk, charge_sign*z, 1, -1);
				double mz0 = pk_pair.second->mz;
				std::set<RichPeakPtr> pk_set = getClosestRichPeaks(spectrum, mz0, true, 0.0);

				if(pk_set.size() == 0)
					continue;

				// Estimate the intensity shift.
				RichPeakPtr pk0 = *(pk_set.begin());

				double intensity_shift = (pk0->intensity - pk_pair.first->intensity) / pk_pair.first->intensity;

				if(intensity_shift > -0.2) {
					count++;
					// A candidate envelop.
					und_pk.insert(std::make_pair(base_pk, z));
				}

			}

			if(count == 0)
				setPeakType(base_pk, "NOISE");

		}
	}

	std::pair<PeakPtr, PeakPtr> SimpleFinder::getPeakRangeByShift( RichPeakPtr base_pk, int charge, int shift, int sulfur /*= -1*/ )
	{
		double pre_mass = calculateMass(base_pk->mz, charge);

		AggregatedIsotopicVariants iso_no_sulfur = this->estimateDistribution(base_pk->mz, charge, 0);
		int max_sulfur = (sulfur == -1 ? guessMaxSulfurNumber(base_pk->mz, charge) : sulfur);
		AggregatedIsotopicVariants iso_max_sulfur = this->estimateDistribution(base_pk->mz, charge, max_sulfur);

		double mz_max_sulfur = base_pk->mz + iso_max_sulfur.getMassDifferenceByShift<peak_intensity>(0, shift)/(float)abs(charge);
		double mz_no_sulfur = base_pk->mz + iso_no_sulfur.getMassDifferenceByShift<peak_intensity>(0, shift)/(float)abs(charge);

		double int_max_sulfur = base_pk->intensity * iso_max_sulfur.getPeakByShift<peak_intensity>(shift)->intensity;
		double int_no_sulfur = base_pk->intensity * iso_no_sulfur.getPeakByShift<peak_intensity>(shift)->intensity;

		PeakPtr pk_max_sulfur = boost::make_shared<Peak>(mz_max_sulfur, int_max_sulfur);
		PeakPtr pk_no_sulfur = boost::make_shared<Peak>(mz_no_sulfur, int_no_sulfur);

		return std::make_pair(pk_max_sulfur, pk_no_sulfur);

	}

	void SimpleFinder::setPeakType( RichPeakPtr pk, const std::string& new_type )
	{
		// 1. Find the iterator of the peak.
		RichPeakListByID& pks_id = spectrum.getPeakListByType<peak_id>();
		RichPeakListByID::iterator id_iter = pks_id.find(pk->id);

		// 2. Project the iterator to type.
		RichPeakListByType::iterator type_iter = spectrum.getPeakContainer().project<peak_type>(id_iter);

		// 2. Modify the type of the element referred by the iterator.
		spectrum.modifyType(type_iter, new_type);
	}

	void SimpleFinder::detectHarmonicCluster(RichPeakPtr base_pk, int charge)
	{
		//std::set<RichPeakPtr> harmonic_cluster;

		// The S/N threshold for monoisotopic peawk.
		double mono_sn_threshold = param.getParameter<double>("signal_noise").first;

		/*double cur_mz = param.getParameter<double>("precursor_mz").first;
		int cur_charge = param.getParameter<int>("precursor_charge").first;*/

		double cur_mass = calculateMass(base_pk->mz, charge);

		int coef = 2;
		// Get the minimum m/z.
		RichPeakListByMZ& pks_mz = spectrum.getPeakListByType<peak_mz>();
		double min_mz = (*(pks_mz.begin()))->mz;

		while(1)
		{
			// Get new m/z based on adjusted charge state.
			
			int new_charge = coef * charge;
			double new_mz = calculateMZ(cur_mass, new_charge, charge);

			coef++;

			if(new_mz < min_mz)
				break;

			// Find the base peak.
			std::set<RichPeakPtr> pk_set = getClosestRichPeaks(spectrum, new_mz, false, mono_sn_threshold, 3e-5);

			if(pk_set.size() == 0)
				continue;
			
			RichPeakPtr base_pk = *(pk_set.begin());
			this->cleanNeighborNoise(spectrum, base_pk);

			if(base_pk->signal_noise < mono_sn_threshold)
				continue;

			bool status = extendEnvelop(spectrum, base_pk, abs(new_charge), false);

			if(status) {
				//harmonic_cluster.insert(base_pk);
				std::cout << "Find harmonic cluster: " << base_pk->mz << " " << new_charge << std::endl;
			}
		}
		
		//return harmonic_cluster;
	}

}