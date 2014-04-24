/*
 * =====================================================================================
 *
 *       Filename:  FragmentationTable.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/29/2012 10:16:24 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/FRAGMENTATION/FragmentationTable.h>

namespace gag
{
	FragmentationTable& FragmentationTable::Instance()
	{
		static FragmentationTable ftt;
		return ftt;
	}

	void FragmentationTable::load(const std::string& filename)
	{

		using boost::property_tree::ptree;
		ptree pt;

		read_xml(filename, pt);
		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.MassShift"))
		{
			std::string value = v.second.get<std::string>("Value");
			int lowerlimit = v.second.get<int>("LowerLimit");
			int upperlimit = v.second.get<int>("UpperLimit");
			//std::pair<int, int>
			mass_loss.insert(std::make_pair(value, std::make_pair(lowerlimit, upperlimit)));
			
			
		}
		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.CleavageTypes"))
		{
			FragmentationParams fp;
			fp.type = v.second.get<std::string>("Name");
			BOOST_FOREACH(ptree::value_type &s, v.second.get_child("Shift"))
			{
				if(s.first == "Value")
					fp.cleavage_shift.push_back(s.second.data());
				else if(s.first == "Dissociation") {
					std::string name = s.second.get<std::string>("Name");
					std::vector<std::string> shift_vals;
					std::pair<ptree::assoc_iterator, ptree::assoc_iterator> ret2 = s.second.equal_range("Value");
					for(ptree::assoc_iterator iter = ret2.first; iter != ret2.second; iter++)
					{
						shift_vals.push_back(iter->second.data());
            fp.cleavage_shift.push_back(iter->second.data());
					}
					fp.dis_shift.insert(std::make_pair(name, shift_vals));
				}
			}

			fragmentation_params.insert(std::make_pair(fp.type, fp));	
			
		}


		pt.erase("parameters.CleavageTypes");

	}

	FragmentationParams FragmentationTable::getFragmentationParams(const std::string& cleavage_type) const
	{
		std::map<std::string, FragmentationParams>::const_iterator i = fragmentation_params.find(cleavage_type);
		return i != fragmentation_params.end()? i->second : FragmentationParams();
	}

	CompositionSigned FragmentationTable::getCompositionSigned(const std::string& str) const
	{
		if(str.empty())
			return CompositionSigned();
		CompositionSigned compo_signed;
		std::string::const_iterator it = str.begin();
		if((*it) == '-') { // Negative
			compo_signed.first = -1;
			compo_signed.second = str.substr(1);
		}	else {
			compo_signed.first = 1;
			compo_signed.second= str;
		}
		return compo_signed;
	}
	
	CompositionShift FragmentationTable::getCleavageShift(const std::string& cleavage_type) const
	{
		std::map<std::string, FragmentationParams>::const_iterator i = fragmentation_params.find(cleavage_type);

		CompositionShift compo_shift;
		if(i != fragmentation_params.end()) {
			std::vector<std::string>::const_iterator const_iter = i->second.cleavage_shift.begin();
			for(; const_iter != i->second.cleavage_shift.end(); const_iter++)
				compo_shift.push_back((*this).getCompositionSigned(*const_iter));
		}
		return compo_shift;
	}

	CompositionShift FragmentationTable::getCleavageShift(const std::string& cleavage_type, const std::string& dis_type) const
	{
		std::map<std::string, FragmentationParams>::const_iterator i = fragmentation_params.find(cleavage_type);
		
		CompositionShift compo_shift;

		if(i != fragmentation_params.end()) {
			std::map<std::string, std::vector<std::string> >::const_iterator j = i->second.dis_shift.find(dis_type);
			
			if(j != i->second.dis_shift.end())				
				for(std::vector<std::string>::const_iterator iter = j->second.begin(); iter != j->second.end(); iter++)
					compo_shift.push_back((*this).getCompositionSigned(*iter));
		} 

		return compo_shift;
	}

	MassLossWindow FragmentationTable::getMassLoss() const
	{
		std::map<std::string, std::pair<int, int> >::const_iterator const_iter = mass_loss.begin();
		MassLossWindow mlw;

		for(; const_iter != mass_loss.end(); const_iter++)
		{
			MassLoss ml;
			ml.loss_compo = Composition(const_iter->first);
			ml.lower = const_iter->second.first;
			ml.upper = const_iter->second.second;
			mlw.push_back(ml);
		}
		return mlw;
	}
}