/********************************************************************
	created:	2012/08/27
	created:	27:8:2012   8:40
	filename: 	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\src\GAGPL\SPECTRUM\IsotopicDistribution.cpp
	file path:	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\src\GAGPL\SPECTRUM
	file base:	IsotopicDistribution
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/
#include "GAGPL/SPECTRUM/IsotopicDistribution.h"
#include "GAGPL/MISC/Param.h"

#include <time.h>
#include <algorithm>
#include <sstream>
#include <fstream>

#include <boost/make_shared.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>

namespace gag
{
	void IsotopicDistribution::createMonoPeak()
	{
		// Get the monoisotopic isotope for each element of the composition.
		monopeak.mz = _compo.getMass();

		const std::map<std::string, int> ele_count = _compo.get();

		std::map<std::string, int>::const_iterator it = ele_count.begin();

		for(; it != ele_count.end(); it++)
		{
			// For each element, get the value of prob^coef.
			PeriodicTable& ptable = PeriodicTable::Instance();
			Isotope iso = ptable.getIsotopeByRelativeShift(it->first);
			monopeak.intensity += it->second * log(iso.abundance);
		}
		monopeak.intensity = exp(monopeak.intensity);

	}

	double IsotopicDistribution::calculatePhiValue(size_t order)
	{
		double phi_value = 0.0;

		// For each element, get the phi factor.
		const std::map<std::string, int> ele_count = _compo.get();

		std::map<std::string, int>::const_iterator it = ele_count.begin();

		// Iterate over all elements.
		for(; it != ele_count.end(); it++)
		{
			phi_value += _iso_const.getNthElementPowerSum(it->first, order) * it->second;
		}

		return phi_value;
	}

	//PowerSumVec IsotopicDistribution::calculatePhiValueVec()
	//{
	//	PowerSumVec phi_value_vec;
	//	// phi(0) = 0.
	//	phi_value_vec.push_back(0);

	//	for(size_t i=1; i<=_var_num; i++)
	//	{
	//		phi_value_vec.push_back(this->calculatePhiValue(i));
	//	}
	//		
	//	return phi_value_vec;
	//}

	EleSymPolyVec IsotopicDistribution::calculateProbabilityVec()
	{
		//clock_t t = clock(); clock_t last_t = t;
		PowerSumVec phi_vec = this->calculatePhiValueVec();
		//t = clock()-last_t; last_t = clock();
		//std::cout << "Step 1: calculating phi value vec: " << (double)t/CLOCKS_PER_SEC << "s\n";

		NewtonGirardFormulae newton(_compo.getMaxNumVariants());
		EleSymPolyVec prob_vec;
		newton.updateParameters(phi_vec, prob_vec);
		//t = clock()-last_t; last_t = clock();
		//std::cout << "Step 2: newton: " << (double)t/CLOCKS_PER_SEC << "s\n";

		for(size_t i=0; i<prob_vec.size(); i++)
		{
			int sign = (i%2 == 0 ? 1 : -1);
			// q(j) = q(0) * e(j) * (-1)^j.
			prob_vec[i] = prob_vec[i] * monopeak.intensity * sign;
		}

		return prob_vec;
	}

	double IsotopicDistribution::calculateModifiedPhiValue(const std::string& symbol, size_t order)
	{
		
		double phi_value = phi_value_vec.at(order) - _iso_const.getNthElementPowerSum(symbol, order) + _iso_const.getNthModifiedElementPowerSum(symbol, order);

		return phi_value;
	}

	PowerSumVec IsotopicDistribution::calculateModifiedPhiValueVec(const std::string& symbol)
	{
		PowerSumVec modified_phi_value_vec;
		// phi(0) = 1.
		modified_phi_value_vec.push_back(0);

		for(size_t i=1; i<=_var_num; i++)
		{
			modified_phi_value_vec.push_back(this->calculateModifiedPhiValue(symbol, i));
		}

		return modified_phi_value_vec;
	}

	std::vector<double> IsotopicDistribution::calculateCenterMassVec(const EleSymPolyVec& prob_vec)
	{
		std::vector<double> mass_vec;

		std::map<std::string, int> ele_count = _compo.get();

		PeriodicTable& ptable = PeriodicTable::Instance();

		NewtonGirardFormulae newton(_compo.getMaxNumVariants());

		std::map<std::string, EleSymPolyVec> esp_map;
		//clock_t t = clock(); 
		for(std::map<std::string, int>::const_iterator it = ele_count.begin();
			it != ele_count.end(); it++)
		{
			const Isotope& mono = ptable.getIsotopeByRelativeShift(it->first);

			PowerSumVec& ps_vec = this->calculateModifiedPhiValueVec(it->first);
			EleSymPolyVec esp_vec;
			newton.updateParameters(ps_vec, esp_vec);

			esp_map.insert(std::make_pair(it->first, esp_vec));
		}
		//t = clock()-t;
		//std::cout << "Modified Phi Vec: " << (double)t/CLOCKS_PER_SEC << "s\n";

		for(size_t i = 0; i <= _var_num; i++)
		{
			std::map<std::string, EleSymPolyVec>::iterator it = esp_map.begin();
			int sign = (i%2 == 0 ? 1 : -1);
			double center_mass = 0.0;
			for(; it != esp_map.end(); it++)
			{
				const Isotope& mono = ptable.getIsotopeByRelativeShift(it->first);

				center_mass += ele_count[it->first] * sign * it->second.at(i) * monopeak.intensity * mono.mass;
			}
			// m(j) = sum(m(jk) * p(jk))/sum(p(jk))
			mass_vec.push_back(center_mass/prob_vec.at(i));
		}

		return mass_vec;
	}

	void IsotopicDistribution::setOrder(int order)
	{
		// Controlling the number of variants, and prevent unnecessary calculation.
		size_t max_num = _compo.getMaxNumVariants();
		if(order == 0) {
			// Default value calculated by formula 2.			
			int heu = std::max((int)ceil(abs(2 * (_compo.getMass() - _compo.getAverageMass()))), 50);
			_var_num = (size_t)heu > max_num ? max_num : (size_t)heu;
		} else if(order > 0)
			_var_num = (size_t)order > max_num ? max_num : (size_t)order;
		else
			throw std::runtime_error("Error: invalid order number.");

		_var_num--;

		// Adding elements from the composition into _iso_const.
		updateIsotopicConstants();
	}

	AggregatedIsotopicVariants IsotopicDistribution::getAggregatedIsotopicVariants(int charge /* = 0 */)
	{
		//clock_t t = clock(); clock_t last_t = t;
		EleSymPolyVec prob_vec = this->calculateProbabilityVec();

		//t = clock()-last_t; last_t = clock();
		//std::cout << "Prob Vec: " << (double)t/CLOCKS_PER_SEC << "s\n";

		std::vector<double> center_mass_vec = this->calculateCenterMassVec(prob_vec);
		//t = clock()-last_t; last_t = clock();
		//std::cout << "Mass Vec: " << (double)t/CLOCKS_PER_SEC << "s\n";

		AggregatedIsotopicVariants peakset;

		avg_mass = 0.0; 
		double sum(0.0);


		for(size_t i=0; i<= _var_num; i++)
		{
			double adjusted_mz = (charge == 0 ? center_mass_vec.at(i): calculateMZ(center_mass_vec.at(i), charge));

			PeakPtr peak = boost::make_shared<Peak>(adjusted_mz, prob_vec.at(i));
			peakset.addPeak(peak);	

			// Code for calculating avg_mass
			avg_mass += adjusted_mz * prob_vec.at(i);
			sum += prob_vec.at(i);
		}

		avg_mass /= sum;

		return peakset;

	}

	void IsotopicDistribution::updateIsotopicConstants()
	{
		const std::map<std::string, int> ele_count = _compo.get();

		std::map<std::string, int>::const_iterator it = ele_count.begin();

		// Build the constants for all elements.
		for(; it != ele_count.end(); it++)
			_iso_const.addElement(it->first);

		_iso_const.updateOrder(_var_num);
	}

	void IsotopicDistribution::updatePhiValueVec()
	{
		phi_value_vec.clear();
		phi_value_vec.push_back(0);
		for(size_t i=1; i<=_var_num; i++)
		{
			phi_value_vec.push_back(this->calculatePhiValue(i));
		}
	}

	double calculateMass( double mz, int charge)
	{	
		// Be careful of the ion mode.
		param::Param& param = param::Param::Instance();
		double electron_mass = param.getParameter<double>("electron_mass").first;

		return abs(charge) * mz - charge * (Composition("H").getMass() -electron_mass);
	}

	double calculateMZ(double mass, int charge, int pre_charge /* = 0 */)
	{
		param::Param& param = param::Param::Instance();
		double electron_mass = param.getParameter<double>("electron_mass").first;

		int coef_h = (pre_charge == 0 ? charge : pre_charge);

		// TBD: this should be controlled by parameter.
		// int coef_e = (pre_charge == 0 ? 1 : 0);
		return (mass + coef_h * (Composition("H").getMass() - electron_mass))/abs(charge);
	}
}