/********************************************************************
	created:	2013/01/10
	created:	10:1:2013   9:07
	filename: 	SimpleFinder.h
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	SimpleFinder
	file ext:	h
	author:		Han Hu
	
	purpose:	The program requires that the monoisotopic peak has to be observed. It will focus on the identification of A+1 and A+2 peak.  The program assumes composition can be explained as x(CH2O) + y(SO3). 
*********************************************************************/

#ifndef GAG_SIMPLEFINDER_H
#define GAG_SIMPLEFINDER_H

#include <boost/shared_ptr.hpp>
#include "GAGPL/SPECTRUM/Envelop.h"
#include "GAGPL/SPECTRUM/EnvelopReference.h"
#include "GAGPL/MISC/Param.h"

namespace gag
{	
	typedef std::map<std::string, double> AveragineFormulae;
	using namespace param;

	class SimpleFinder
	{
	//typedef std::vector<RichList> IslandList;
	private:
		Param& param;
		PeriodicTable& pt;

		// No sulfate.
		AveragineFormulae no_sulfur_ave;
		// High sulfate
		AveragineFormulae max_sulfur_ave;

		EnvelopReference env_ref;

		std::set<EnvelopPtr> env_pool;

		// Undetermined peaks and corresponding maximum charge states.
		std::multimap<RichPeakPtr, int> und_pk;

		RichList& spectrum;

		double precursor_mass;

	private:
		void run();

		//void findNextEnvelop(RichList& pk_list, RichPeakListByIntensity::iterator iter_intensity);
		void findNextEnvelop(RichList& pk_list, RichPeakPtr base_pk, bool raw_spec = false);
		
		// This function will only collect as more candidate envelops as possible. Notice that the sign of the charge is 1.
		bool extendEnvelop(RichList& pk_list, RichPeakPtr base_pk, int charge , bool raw_spec);

		// The calculation of fitting score refers to MS-Deconv. The fitting score of the envelop will be updated and returned.
		double calculateFittingScore(EnvelopPtr env);
		double calculateFittingScore(std::set<EnvelopPtr>& env_set);
		
		AggregatedIsotopicVariants estimateDistribution(double mz, int charge, bool sulfur = false);
		// Given the envelop information, estimate the theoretical distribution
		AggregatedIsotopicVariants estimateDistribution(double mz, int charge, int sulfur_num);

		// Find all possible peaks within a range configured by parameter file. Currently just find the most abundant one to prevent excessive calculation.
		std::set<RichPeakPtr> getClosestRichPeaks(RichList& pk_list, double expected_mz, bool raw_spec, double sn_threshold, double error = -1.0);

		// Optimize the parameter using peaks only matched to the current shift.
		bool isRedundantEnvelop(EnvelopPtr new_env, EnvelopPtr old_env, int old_shift);

		AveragineFormulae generateFormulae(double sulfur_num);

		// Estimate sulfur number and calculate the shift difference for peaks in the same envelop.
		void updateEnvelopParameter(EnvelopPtr env);

		//void removeEnvelop(EnvelopPtr env);

		// TBD: Convert the signal/noise of current peak into corresponding window error.
		double getWindowError(RichPeakPtr pk, double error = -1.0);

		// This function will clean all the apodization peaks surrounding significant peak.
		void cleanNeighborNoise(RichList& pk_list, RichPeakPtr pk/*, std::set<RichPeakPtr>& blacklist*/);
		
		//void setNoisePeak(RichPeakPtr pk);
		void setNoisePeak(RichPeakPtr pk);
		void setPeakType(RichPeakPtr pk, const std::string& new_type);

		//void splitPeak(std::vector<EnvEntry>& old_entries, EnvelopPtr new_env, );

		// Update the overall fitting score of the whole tree. Be careful that the updating will be strictly limited to the envelops involved, no others will be involved. The parameter of the new envelop might be altered. A hyper scaling-factor will be append to the old entries.
		double updateEnvelopTree(std::vector<EnvEntry>& old_entries, std::vector<EnvEntry>& added_entries, EnvEntry& test_entry, RichPeakPtr pk);

		// Calculate the linear relationship between the exp cluster and theo cluster. The method will take the vectors and merge them together based on the actual overlapping situation.
		double calculateLinearAssociationScore(std::vector<EnvEntry>& entries);

		// Calculate the scaling factor of the actual peaks over the theo_peaks.
		double calculateScalingRatio(std::vector<EnvEntry>& entries);
		double calculateScalingRatio(EnvEntry& entry);

		std::pair<double, double> calculateAccumutativeAbundance();

		void updatePeakAbuandance(std::vector<EnvEntry>& entries, RichPeakPtr pk);
		
		int guessMaxSulfurNumber(double mz, int charge_state);

		void acceptEnvelop(EnvelopPtr env);

	public:
		SimpleFinder(RichList& spec);
	
		// Identify if this is a true envelop, including misc envelop.
		bool isNoiseEnvelop(EnvelopPtr env);
		// mz value and corresponding envelop.
		std::multimap<RichPeakPtr, EnvelopPtr> getEnvelops(EntryStatus status = OLD);
		// To be implemented. The function will return the results that store the information for library matching.
		RichList getMonoisotopicPeakList();

		// Get the peak/charge the pattern of which can not be observed. The charge is the maximum value it is possible to achieve.
		inline std::multimap<RichPeakPtr, int>& getUndeterminedPeaks()
		{
			return und_pk;
		}

		// Optimization of current envelop set.
		void optimizeEnvelopSet(std::set<EnvelopPtr>& env_set);

		void optimizeSuspiciousPeaks();

		void printEnvelop(EnvelopPtr env);

		// If the sulfur number is configured, use the number instead.
		// Otherwise, use the default value.
		std::pair<PeakPtr, PeakPtr> getPeakRangeByShift(RichPeakPtr base_pk, int charge, int shift, int sulfur = -1);

		// Once an envelop has been identified, check the corresponding harmonic clusters.
		void detectHarmonicCluster(RichPeakPtr base_pk, int charge);

	};

}




#endif /* GAG_SIMPLEFINDER_H */