/********************************************************************
	created:	2012/11/16
	created:	16:11:2012   16:46
	filename: 	EnvelopReference.h
	file path:	GAG\src\GAGPL\SPECTRUM
	file base:	EnvelopReference
	file ext:	h
	author:		Han Hu
	
	purpose:	Storing the mapping relationship between envelop id and 
	peak id.
*********************************************************************/

#ifndef GAG_ENVELOPREFERENCE_H
#define GAG_ENVELOPREFERENCE_H

#include "GAGPL/SPECTRUM/Envelop.h"
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <set>

namespace gag
{
	// Any information containing the peak can be append here.
	struct EnvEntry
	{
		EnvelopPtr env;
		RichPeakPtr pk;

		int pk_shift;

		InfoPeakPtr info;

		// Indicate if the entry is empty.
		bool status;
	
		EnvEntry(RichPeakPtr p, EnvelopPtr e, InfoPeakPtr i, int shift)
			: pk(p), env(e), pk_shift(shift), info(i), status(true)
		{}
		EnvEntry() : status(false) {}

		inline unsigned int getPeakID() const
		{
			return pk->id;
		}
		inline unsigned int getEnvelopID() const
		{
			return env->id;
		}
		inline int getShift() const
		{
			return pk_shift;
		}
		inline bool isEmpty() const
		{
			return !status;
		}
		inline EntryStatus getEntryStatus() const
		{
			return info->entry_status;
		}
		inline double getFittingScore() const
		{
			return env->fitting_score;
		}
		inline EnvelopStatus getEnvelopStatus() const
		{
			return env->env_status;
		}

		static bool scoreLarger(const EnvEntry& env_entry1, const EnvEntry& env_entry2)
		{
			return env_entry1.getFittingScore() > env_entry2.getFittingScore();
		}
	};

	using namespace ::boost;
	using namespace ::boost::multi_index;
	
	struct env_id{};
	struct pk_shift{};
	struct envelop_status{};
	
	typedef boost::multi_index_container<
		EnvEntry,
		indexed_by<
			ordered_non_unique<
			composite_key<
				EnvEntry,
				// Partial search should also be useful.
				const_mem_fun<EnvEntry, unsigned int, &EnvEntry::getEnvelopID>,
				const_mem_fun<EnvEntry, int, &EnvEntry::getShift> >
			>,
			ordered_non_unique<
			composite_key<
				EnvEntry,
				const_mem_fun<EnvEntry, unsigned int, &EnvEntry::getPeakID>,
				const_mem_fun<EnvEntry, EntryStatus, &EnvEntry::getEntryStatus> >
			>,
			ordered_non_unique<
			composite_key<
				EnvEntry,
				const_mem_fun<EnvEntry, int, &EnvEntry::getShift>,
				const_mem_fun<EnvEntry, EnvelopStatus, &EnvEntry::getEnvelopStatus> >
			>
		>
	> EnvDictionary;
	
	typedef nth_index<EnvDictionary, 0>::type EnvDictByEnvID;
	typedef nth_index<EnvDictionary, 1>::type EnvDictByPeakID;
	typedef nth_index<EnvDictionary, 2>::type EnvDictByStatus;

	class EnvelopReference
	{
	public:
		void addDictionaryReference(RichPeakPtr pk, EnvelopPtr env, int shift, InfoPeakPtr info);

		inline EnvDictionary& getDictionary()
		{
			return env_dict;
		}

		std::vector<EnvEntry> getEntryByPeak(RichPeakPtr pk);
		std::vector<EnvEntry> getEntryByEnvelop(EnvelopPtr env);
		std::set<RichPeakPtr> getPeaksByEnvelop(EnvelopPtr env);

		//void removeEnvelop(EnvelopPtr env);
		
		// TBD: the peaks should be non-redundant.
		RichList getRichPeakList(std::set<EnvelopPtr>& env_set);
		RichList getRichPeakList(std::vector<EnvEntry>& env_entry_set);

		RichPeakPtr getBasePeakForEnvelop(EnvelopPtr env);
		
		EnvEntry getEntryByShift(EnvelopPtr env, int shift);
		// Modified version of getting peak list.
		std::vector<EnvEntry> getPeaksByShift(EnvelopPtr env, int shift);
		
		// This method will calculate the sum of the peak intensities for the 
		// same shift.
		double getExperimentalAbundanceByShift(EnvelopPtr env, int shift);

		// Get the calculated theoretical peak.
		PeakPtr getTheoreticalPeak(EnvelopPtr env, int shift);

		// Get theoretical abundance using information from base peak.
		double getTheoreticalAbundance(EnvEntry& pk_entry);
		void printEnvelopInformation(EnvelopPtr env);

		// If the peak is included the given envelop.
		bool isEnvelopIncluded(RichPeakPtr pk, EnvelopPtr env);

		// Get all envelops where the function isNewEnvelop() of the entry info item says true.
		std::vector<EnvelopPtr> getOccurredEnvelops(RichPeakPtr pk, EntryStatus new_occur);
		std::vector<EnvEntry> getOccurredEnvelopEntries(RichPeakPtr pk, EntryStatus new_occur);
		std::vector<EnvEntry> getBaseEntriesByStatus(EnvelopStatus status = TP);

		// Get the set of shift values.
		std::set<int> getShiftSet(EnvelopPtr env);

	private: 
		EnvDictionary env_dict;
	};



}

#endif /* GAG_ENVELOPREFERENCE_H */