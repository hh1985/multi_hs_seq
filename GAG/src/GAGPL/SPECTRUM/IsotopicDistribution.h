/********************************************************************
	created:	2012/08/23
	created:	23:8:2012   16:19
	filename: 	IsotopicDistribution.h
	file path:	GAGPL\SPECTRUM
	file base:	IsotopicDistribution
	file ext:	h
	author:		Han Hu
	
	purpose:	Calculating theoretical isotopic distribution.
*********************************************************************/
#ifndef GAG_ISOTOPICDISTRIBUTION_H
#define GAG_ISOTOPICDISTRIBUTION_H

#include "GAGPL/CHEMISTRY/Composition.h"
#include "GAGPL/SPECTRUM/IsotopicConstants.h"
#include "GAGPL/SPECTRUM/PeakList.h"

namespace gag
{
	// The aggregated masses and corresponding probabilities.
	// Note that the double type for probability should be converted 
	// to float type in order for the abundance matching procedure.

	// Based on the mode, charge state and mz, calculate the original mass.
	double calculateMass(double mz, int charge);

	// Calculate mz from mode, mass and charge. By default the number of proton equals to the charge state value.
	double calculateMZ(double mass, int charge, int pre_charge = 0);

	// Implementation of BRAIN algorithm.
	class IsotopicDistribution
	{
	public:
		IsotopicDistribution()
			: _iso_const(IsotopicConstants::Instance()) {}
		//IsotopicDistribution(size_t order)
		//: _iso_const(IsotopicConstants::Instance()), _var_num(order) {}
		IsotopicDistribution(Composition& compo, int order = 0)
			: _iso_const(IsotopicConstants::Instance()), _compo(compo)
		{
			// Implicitly build the distribution.
			createMonoPeak();
			setOrder(order);
			updatePhiValueVec();
		}

		void update(Composition& compo, int order=-1)
		{
			_compo = compo;
			createMonoPeak();
			setOrder(order);
			updatePhiValueVec();
		}

		EleSymPolyVec calculateProbabilityVec();

		std::vector<double> calculateCenterMassVec(const EleSymPolyVec& prob_vec);

		// The peak cluster can be adjusted by charge. Notice that the charge here can be either positive or negative.
		AggregatedIsotopicVariants getAggregatedIsotopicVariants(int signed_charge = 0);

		inline double getAverageMass() const
		{
			return avg_mass;
		}

		inline std::string getCompositionString() const
		{
			return _compo.getCompositionString();
		}

		inline void clear()
		{
			_iso_const.clear();
		}

	private:
		// phi_vec[0] should be assigned 0.
		IsotopicConstants& _iso_const;

		// Composition of the molecular formula.
		Composition _compo;

		Peak monopeak;

		// The average mass calculated based on returned isotopic peaks. 
		// It can be used to control the number of calculated peaks.
		double avg_mass;

		// The number of isotopic variants.
		size_t _var_num;

		// Updated.
		PowerSumVec phi_value_vec;
		//PowerSumVec modified_phi_value_vec;

		// The phi_vec[index] has been assigned, just return the value,
		// otherwise, calculate it from the scratch.
		double calculatePhiValue(size_t order);
		inline PowerSumVec calculatePhiValueVec()
		{
			return phi_value_vec;
		}
		// Updated: Calculating phi value vec and modified phi value vec at the
		// same time.
		void updatePhiValueVec();

		double calculateModifiedPhiValue(const std::string& symbol, size_t order);
		PowerSumVec calculateModifiedPhiValueVec(const std::string& symbol);


		// Calculating the information of monoisotopic peak. 
		void createMonoPeak();

		void setOrder(int order);

		void updateIsotopicConstants();


	};

}

#endif /* GAG_ISOTOPICDISTRIBUTION_H */