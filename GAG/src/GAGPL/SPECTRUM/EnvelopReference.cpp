/********************************************************************
	created:	2012/11/24
	created:	24:11:2012   11:44
	filename: 	EnvelopReference.cpp
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	EnvelopReference
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/
#include "GAGPL/SPECTRUM/EnvelopReference.h"
#include <algorithm>

namespace gag
{

	void EnvelopReference::addDictionaryReference(RichPeakPtr pk, EnvelopPtr env, int shift, InfoPeakPtr info)
	{
		EnvEntry env_entry(pk, env, info, shift);
		env_dict.insert(env_entry);
	}

	std::vector<EnvEntry> EnvelopReference::getEntryByPeak( RichPeakPtr pk )
	{
		std::vector<EnvEntry> entry_vec;
		EnvDictByPeakID& env_by_pid = env_dict.get<1>();

		for(EnvDictByPeakID::iterator env_iter = env_by_pid.equal_range(pk->id).first; env_iter != env_by_pid.equal_range(pk->id).second; env_iter++)
			entry_vec.push_back(*env_iter);

		return entry_vec;
	}

	std::vector<EnvEntry> EnvelopReference::getEntryByEnvelop( EnvelopPtr env )
	{
		// Using the number of sulfate to calculate the distribution.

		std::vector<EnvEntry> entry_vec;
		EnvDictByEnvID& env_by_eid = env_dict.get<0>();
		std::pair<EnvDictByEnvID::iterator, EnvDictByEnvID::iterator> p = env_dict.equal_range(env->id);
		for(EnvDictByEnvID::iterator iter = p.first; iter != p.second; iter++)
			entry_vec.push_back(*iter);

		return entry_vec;
	}

	EnvEntry EnvelopReference::getEntryByShift( EnvelopPtr env, int shift )
	{
		EnvDictionary::iterator iter = env_dict.find(boost::make_tuple(env->id, shift));

		return iter != env_dict.end() ? *iter : EnvEntry();
	}

	RichPeakPtr EnvelopReference::getBasePeakForEnvelop( EnvelopPtr env )
	{
		EnvEntry entry = this->getEntryByShift(env, 0);
		return entry.pk;
	}

	PeakPtr EnvelopReference::getTheoreticalPeak( EnvelopPtr env, int shift )
	{
		// 1. Get the base peak of the envelop.
		EnvEntry base_entry = this->getEntryByShift(env,0);
    double base_intensity = env->getTheoreticalPeak(0)->intensity;

		// 2. adjust the mz using base peak information.
		double mz = base_entry.pk->mz + env->getTheoreticalMZDistance(0, shift);
		double intensity = env->getTheoreticalPeak(shift)->intensity / base_intensity;

		// 3. Update the intensity to the real scale.
		intensity = base_entry.info->adjusted_abundance * intensity;

		return boost::make_shared<Peak>(mz, intensity);

	}

	void EnvelopReference::printEnvelopInformation(EnvelopPtr env)
	{
		std::vector<EnvEntry> entry_vec = this->getEntryByEnvelop(env);

		std::cout << "Envelop ID: " << env->id << env->fitting_score << std::endl;

		BOOST_FOREACH(EnvEntry& entry, entry_vec)
		{
			std::cout << "Shift: " << entry.getShift() << " MZ: " << entry.pk->mz << " Adjusted ABD: " << entry.info->adjusted_abundance << std::endl;
		}
	}

	bool EnvelopReference::isEnvelopIncluded( RichPeakPtr pk, EnvelopPtr env )
	{
		std::vector<EnvEntry> entry_vec = this->getEntryByPeak(pk);
		BOOST_FOREACH(EnvEntry& entry, entry_vec)
		{
			if(entry.env->id == env->id)
				return true;
		}
		return false;
	}

	double EnvelopReference::getTheoreticalAbundance( EnvEntry& pk_entry )
	{
		//RichPeakPtr base_pk = this->getBasePeakForEnvelop(entry.env);
		EnvEntry base_entry = this->getEntryByShift(pk_entry.env, 0);
		return base_entry.info->adjusted_abundance * pk_entry.env->getTheoreticalPeak(pk_entry.getShift())->intensity;
	}

	RichList EnvelopReference::getRichPeakList( std::set<EnvelopPtr>& env_set )
	{
		// Get all the peaks with entry status showing as OLD or NEW
		//RichList pk_list;
		std::set<RichPeakPtr> pk_set;
		BOOST_FOREACH(EnvelopPtr env, env_set)
		{
			std::vector<EnvEntry> entry_vec = this->getEntryByEnvelop(env);
			BOOST_FOREACH(EnvEntry pk_entry, entry_vec)
				pk_set.insert(pk_entry.pk);
		}

		RichList pk_list;
		BOOST_FOREACH(RichPeakPtr pk, pk_set)
			pk_list.addPeak(pk);
		return pk_list;
	}

	RichList EnvelopReference::getRichPeakList(std::vector<EnvEntry>& env_entry_vec)
	{
		std::set<RichPeakPtr> pk_set;
		BOOST_FOREACH(EnvEntry& env_entry, env_entry_vec)
		{
			// Get all the peaks belongs to the same envelop.
			std::vector<RichPeakPtr> temp_set;
			std::vector<EnvEntry> entry_vec = this->getEntryByEnvelop(env_entry.env);
			BOOST_FOREACH(EnvEntry pk_entry, entry_vec)
				temp_set.push_back(pk_entry.pk);

			std::set_union(pk_set.begin(), pk_set.end(),temp_set.begin(), temp_set.end(), std::inserter(pk_set, pk_set.begin()));
		}

		RichList pk_list;
		BOOST_FOREACH(RichPeakPtr pk, pk_set)
			pk_list.addPeak(pk);
		return pk_list;
	}
	/*void EnvelopReference::removeEnvelop( EnvelopPtr env )
	{
		EnvDictByEnvID& env_by_eid = env_dict.get<0>();
		std::pair<EnvDictByEnvID::iterator, EnvDictByEnvID::iterator> p = env_by_eid.equal_range(env->id);
		for(EnvDictByEnvID::iterator iter = p.first; iter != p.second; iter++)
			env_by_eid.erase(iter);
	}*/

	std::vector<EnvelopPtr> EnvelopReference::getOccurredEnvelops( RichPeakPtr pk, EntryStatus new_occur )
	{
		std::vector<EnvelopPtr> env_set;
		EnvDictByPeakID& env_by_pid = env_dict.get<1>();
		std::pair<EnvDictByPeakID::iterator, EnvDictByPeakID::iterator> p = env_by_pid.equal_range(boost::make_tuple(pk->id, new_occur));

		for(EnvDictByPeakID::iterator iter = p.first; iter != p.second; iter++)
			env_set.push_back(iter->env);

		return env_set;
	}

	std::vector<EnvEntry> EnvelopReference::getOccurredEnvelopEntries( RichPeakPtr pk, EntryStatus new_occur )
	{
		std::vector<EnvEntry> env_set;
		EnvDictByPeakID& env_by_pid = env_dict.get<1>();
		std::pair<EnvDictByPeakID::iterator, EnvDictByPeakID::iterator> p = env_by_pid.equal_range(boost::make_tuple(pk->id, new_occur));

		for(EnvDictByPeakID::iterator iter = p.first; iter != p.second; iter++)
			env_set.push_back(*iter);

		return env_set;
	}

	std::set<int> EnvelopReference::getShiftSet( EnvelopPtr env )
	{
		std::set<int> shift_set;
		std::vector<EnvEntry> pk_entry_vec = this->getEntryByEnvelop(env);
		BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec)
			shift_set.insert(pk_entry.getShift());

		return shift_set;
	}

	double EnvelopReference::getExperimentalAbundanceByShift( EnvelopPtr env, int shift )
	{
		double intensity = 0.0;
		std::vector<EnvEntry> pk_set = this->getPeaksByShift(env, shift);
		BOOST_FOREACH(EnvEntry& pk_entry, pk_set)
			intensity += pk_entry.info->adjusted_abundance;

		return intensity;
	}

	std::vector<EnvEntry> EnvelopReference::getPeaksByShift( EnvelopPtr env, int shift )
	{
		std::vector<EnvEntry> pk_entry_vec;
		EnvDictByEnvID& env_by_eid = env_dict.get<0>();
		std::pair<EnvDictByEnvID::iterator, EnvDictByEnvID::iterator> p = env_by_eid.equal_range(boost::make_tuple(env->id, shift));

		for(EnvDictByEnvID::iterator iter = p.first; iter != p.second; iter++)
			pk_entry_vec.push_back(*iter);

		return pk_entry_vec;
	}

	std::set<RichPeakPtr> EnvelopReference::getPeaksByEnvelop( EnvelopPtr env )
	{
		std::set<RichPeakPtr> pk_set;
		EnvDictByEnvID& env_by_eid = env_dict.get<0>();
		std::pair<EnvDictByEnvID::iterator, EnvDictByEnvID::iterator> p = env_by_eid.equal_range(env->id);
		for(EnvDictByEnvID::iterator iter = p.first; iter != p.second; iter++)
			pk_set.insert(iter->pk);

		return pk_set;

	}

	std::vector<EnvEntry> EnvelopReference::getBaseEntriesByStatus( EnvelopStatus status )
	{
		std::vector<EnvEntry> env_entry_vec;
		EnvDictByStatus& env_by_status = env_dict.get<2>();

		std::pair<EnvDictByStatus::iterator, EnvDictByStatus::iterator> p = env_by_status.equal_range(boost::make_tuple(0));
		while(p.first != p.second) {
			if(p.first->env->env_status == TP)
				env_entry_vec.push_back(*(p.first));
			p.first++;
		}

		return env_entry_vec;
	}


	//std::set<EnvelopPtr> EnvelopReference::getConnectingEnvelops( EnvelopPtr env )
	//{
	//	std::set<EnvelopPtr> env_vec;

	//	std::vector<EnvEntry> entry_vec = this->getEntryByEnvelop(env);

	//	BOOST_FOREACH(EnvEntry& entry, entry_vec)
	//	{
	//		// Get all connecting envelops for current peak.
	//		std::vector<EnvEntry> pk_entry_vec = this->getEntryByPeak(entry->pk);
	//		BOOST_FOREACH(EnvEntry& pk_entry, pk_entry_vec)
	//		{
	//			std::set<EnvelopPtr>::iterator env_iter = env_vec.find(pk_entry);
	//			if(env_iter != env_vec.end())
	//				env_vec.insert(*env_iter);
	//		}
	//		
	//	}
	//	return env_vec;
	//}

	//std::set<RichPeakPtr> EnvelopReference::getSharedPeaks( EnvelopPtr env1, EnvelopPtr env2 )
	//{
	//	std::set<RichPeakPtr> pk_vec;

	//	std::set<RichPeakPtr> pk_vec1 = this->getPeaksbyEnvelop(env1);
	//	std::set<RichPeakPtr> pk_vec2 = this->getPeaksbyEnvelop(env2);

	//	// Get set_intersection.
	//	std::set_intersection(pk_vec1.begin(), pk_vec1.end(), pk_vec2.begin(), pk_vec2.end(), std::inserter(pk_vec, pk_vec.begin()));

	//	return pk_vec;

	//}



}