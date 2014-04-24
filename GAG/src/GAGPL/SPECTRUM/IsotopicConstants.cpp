#include "GAGPL/SPECTRUM/IsotopicConstants.h"

#include <numeric>
#include <stdexcept>

namespace gag 
{
	using namespace std;
	PolyParam IsotopicConstants::updateCoefsOfElement( const Element& e, bool has_mass /*= false*/ )
	{
		// Get the coefficients of given element
		IsotopeSetSequential::reverse_iterator it = e.isotopes.rbegin();

		// Get the highest order of isotope neutron shift.
		size_t max_n = it->neutron_shift;

		// temp_vec stores the polynomial coef.
		std::vector<double> temp_vec;

		for(; it != e.isotopes.rend(); it++)
		{
			size_t current_order = max_n - it->neutron_shift;

			double coef = (has_mass == true ? it->mass : 1.0);

			// Empty shift.
			if(current_order > temp_vec.size()) {

				for(size_t i=temp_vec.size(); i<current_order; i++)
				{
					temp_vec.push_back(0.0);					
				}
				temp_vec.push_back(it->abundance * coef);

			} else if(current_order == temp_vec.size()){

				temp_vec.push_back(it->abundance * coef);		
			} else {

				throw std::runtime_error("The list of neutron shift is not ordered.");
			}

		}

		VietesFormulae viete(temp_vec);
		EleSymPolyVec esp_vec = viete.getElementarySymmetricFunctionFromCoef();
		NewtonGirardFormulae newton(temp_vec.size()-1);
		PowerSumVec ps_vec;

		newton.updateParameters(ps_vec, esp_vec);	

		return make_pair(esp_vec, ps_vec);
	}

	void IsotopicConstants::updateOrder( size_t order )
	{
		if(order <= _var_num) return;

		_var_num = order;
		update();
	}

	void IsotopicConstants::update()
	{
		ElementConstants::iterator it = _ele_const.begin();
		for(; it != _ele_const.end(); it++)
		{
			// No need to update.
			if(_var_num <= it->second.nth_order)
				continue;

			for(size_t i = it->second.nth_order+1; i <= _var_num; i++)
			{
				it->second.ele_param.first.push_back(0.0);
				it->second.ele_mass_param.first.push_back(0.0);
			}

			// Updating nth_order.
			it->second.nth_order = it->second.ele_param.first.size()-1;
			NewtonGirardFormulae newton(it->second.nth_order);
			newton.updateParameters(it->second.ele_param);
			newton.updateParameters(it->second.ele_mass_param);
		}
	}

	void IsotopicConstants::addElement(const std::string& symbol )
	{
		// If element has been added, ignore the rest part.
		ElementConstants::iterator it = _ele_const.find(symbol);
		if(it!=_ele_const.end())
			return;

		Element ele = ptable.getElementBySymbol(symbol);

		PhiConstants phi_const;
		phi_const.nth_order = ele.isotopes.rbegin()->neutron_shift;;
		phi_const.ele_param = this->updateCoefsOfElement(ele);
		phi_const.ele_mass_param = this->updateCoefsOfElement(ele, true);

		_ele_const.insert(std::make_pair(symbol, phi_const));
	}

	double IsotopicConstants::getNthElementPowerSum( const std::string& symbol, size_t order )
	{
		ElementConstants::iterator it = _ele_const.find(symbol);

		if(it == _ele_const.end())
			throw runtime_error("The specified element is not found!");

		PowerSumVec& ps_vec = it->second.ele_param.second;

		if(ps_vec.size() <= order )
			throw runtime_error("The order number is invalid! Update it first");

		return ps_vec.at(order);
	}

	double IsotopicConstants::getNthModifiedElementPowerSum( const std::string& symbol, size_t order )
	{
		ElementConstants::iterator it = _ele_const.find(symbol);

		if(it == _ele_const.end())
			throw runtime_error("The specified element is not found!");

		PowerSumVec& ps_vec = it->second.ele_mass_param.second;

		if(ps_vec.size() <= order )
			throw runtime_error("The order number is invalid! Update it first");

		return ps_vec.at(order);
	}
}


