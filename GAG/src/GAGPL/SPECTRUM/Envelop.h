/********************************************************************
	created:	2012/11/16
	created:	16:11:2012   15:31
	filename: 	Envelop.h
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	Envelop
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef GAG_ENVELOP_H
#define GAG_ENVELOP_H

#include "GAGPL/SPECTRUM/PeakList.h"
#include "GAGPL/SPECTRUM/IsotopicDistribution.h"
#include <map>

namespace gag
{
	// In order to uniquely specify an envelop, the base pk id,
	// charge states and averagine model has been defined.
	struct Envelop;
	typedef boost::shared_ptr<Envelop> EnvelopPtr;

	// New -- the start of new envelop.
	// OLD -- the extend of old envelop.
	// FALSE -- not an envelop.
	// NOISE -- not an envelop.
	// The member entry_status is set to be FALSE in the beginning. Once it is recognized as start of candidate envelop, the status is set as new, and as old for non-start candidate envelop.
	enum EntryStatus {NEW, OLD, FALSE, NOISE};
	
	// Any temporary information about one peak can be appended here!!!
	struct InfoPeak
	{
		//RichPeakPtr pk_ptr;
		// Can only be in [0,1]
		// double lambda_mz;
		// Can be < 0 (High sulfate) or > 1 (Low sulfate)
		// double lambda_abd;
		// Confidence score associated with the peak position.
		double prob;
		
		EntryStatus entry_status;

		// Estimated abundance after splitting (absolute value). It is initialized as experimental abundance.
		double adjusted_abundance;
		// Once the adjusted_abundance is decided, the relative_shift should be calculate.
		double relative_shift;

		// Estimated abundance from lambda. This is a absolute value.
		//double expected_abundance;

		// Relative share in the envelop.
		double relative_share;

		// Like "Sulfur", "Normal", "C+H"
		std::string status;

		// Should be abandoned in the future.
		InfoPeak(double a_abd, /*double param_mz, double param_abd, */ double p)
			: adjusted_abundance(a_abd), relative_shift(0.0), prob(p), relative_share(1.0), status("Normal")
		{}
		InfoPeak(double a_abd, const std::string& stat, EntryStatus status)
			: adjusted_abundance(a_abd), relative_shift(0.0), status(stat), entry_status(status)
		{}

	};

	typedef boost::shared_ptr<InfoPeak> InfoPeakPtr;

	class EnvelopBoundary
	{
	public:
		EnvelopBoundary(const AggregatedIsotopicVariants& ave1, const AggregatedIsotopicVariants ave2)
			: bound(std::make_pair(ave1, ave2))
		{}
		EnvelopBoundary() {}

		// No sulfate model.
		inline AggregatedIsotopicVariants& getUpperBound()
		{
			return bound.first;
		}
		// Sulfate model.
		inline AggregatedIsotopicVariants& getLowerBound()
		{
			return bound.second;
		}
		PeakPtr getTheoreticalPeak(double sulfur_num, int shift);
		double getMZDistance(double sulfur_num, int ori_shift, int new_shift);
	private:
		std::pair<AggregatedIsotopicVariants, AggregatedIsotopicVariants> bound;
	};

	enum EnvelopStatus {TP, FP, MISC};

	struct Envelop
	{

		// Envelop ID: has to be non-negative.
		// ID 0 is set as initial value.
		unsigned int id;

		// The base peak information can be given.
		// This information should be maintained by the reference table.
		// std::map<int, InfoPeak> pk_cluster;

		// The value can only be positive. The sign is set in param.
		int charge_state;

		// lambda_mz and confidence.
		//double lambda_mz;
		//double lambda_conf;

		// New parameters: scaling factor and sulfur number.
		double scaling_factor;
		int sulfur_num;

		//std::pair<int, int> sulfur_range;

		// Fitting score.
		double fitting_score;

		// If this is an suspicious envelop.
		// bool suspicious;

		EnvelopStatus env_status;

		int max_shift;

		// Notice: The cost of calculating theoretical distrubtion is high!
		AggregatedIsotopicVariants theo_dist;
		
		// EntryStatus env_status;

		// The boundary is defined by the upper and lower theoretical peak clusters.
		// EnvelopBoundary env_bound;
		
		//inline EnvelopBoundary& getBoundary()
		//{
		//	return env_bound;
		//}

		//std::pair<int, int> getShiftRange();
		
		// This is a relative value. Use the sulfur_num as parameter to estimate the averagine and calculate the peak information.
		PeakPtr getTheoreticalPeak(int shift);

		std::map<int, PeakPtr> getTheoreticalPeaks();

		int getMaxSulfurNum();

		inline double getTheoreticalMZDistance(int ori_shift, int new_shift)
		{
			return this->getTheoreticalPeak(new_shift)->mz - this->getTheoreticalPeak(ori_shift)->mz;
		}

		friend EnvelopPtr createEnvelop(int charge);

		// A simpler version.
		Envelop(unsigned int eid, unsigned int charge)
			: id(eid), charge_state(charge), scaling_factor(1.0), sulfur_num(0), env_status(FP), fitting_score(0.0)
		{
		}

		//void printEnvelop();
	};

}


#endif /* GAG_ENVELOP_H */