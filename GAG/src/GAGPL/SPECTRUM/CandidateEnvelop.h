/********************************************************************
	created:	2012/09/06
	created:	6:9:2012   10:31
	filename: 	CandidateEnvelops.h
	file path:	GAGPL\SPECTRUM
	file base:	CandidateEnvelops
	file ext:	h
	author:		Han Hu
	
	purpose:	The class can generate a list of candidate envelops.
*********************************************************************/
#ifndef GAG_CANDIDATEENVELOP_H
#define GAG_CANDIDATEENVELOP_H

#include <GAGPL/SPECTRUM/Spectrum.h>
#include <boost/bimap/multiset_of.hpp>
#include <map>

namespace gag
{
	struct PeakNode
	{
		// The score of each node is independent from others.
		float score;
		// It is necessary to keep the copy, since there might be area split process.
		RawPeak pk_ref;
		//PeakNode* prev_node;
		PeakNode(RawPeak& pk)
			:pk_ref(pk), score(0.0) {}
	};

	class Envelop
	{
	private:
		size_t _id;
		int _charge;
		// When the key value is 0, it is the base peak.
		std::multimap<int, PeakNode> shift_map;
	public:
		Envelop() {}
		Envelop(int charge)
			: _charge(charge) {}

		size_t getID() const;
		int getChargeState();
		RawPeak& getPeakByShift(int shift);
		double getPeakMass(int shift);
		int getMinShift();
		int getMaxShift();
		float getMaxScore();
		// Given the shift value, the function will try to find candidate isotopic peaks.
		void exploreIsotopicPeak(int shift);
		void addIsotopicPeak(int shift, PeakNode& pk);
		bool containIsotopicPeak(int shift);
		bool empty();

		// From the minimum to the maximum, if any of the peak is missing, filled with 0.
		// The value might be modified, so use a copy here.
		std::multimap<int, PeakNode>& getPeakNodes();
	};

	// The left key is the envelop id and the right key is the peak id.
	typedef boost::bimap<
		tagged<multiset_of<size_t>, env_id>,
		tagged<multiset_of<size_t>, pk_id>,
	> EnvelopReference;

	typedef EnvelopReference::value_type ClassifiedPeak;

}


#endif	/* GAG_CANDIDATEENVELOPS_H */