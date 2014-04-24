/********************************************************************
	created:	2013/05/18
	created:	18:5:2013   17:13
	filename: 	SequenceSpace.h
	file path:	GAGPL\GAGLIBRARY
	file base:	SequenceSpace
	file ext:	h
	author:		Han Hu
	
	purpose:	The module is responsible for comparing different types
						of calculation methods and output them with given format.
*********************************************************************/

#ifndef GAG_SEQUENCESPACE_H
#define GAG_SEQUENCESPACE_H


#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/SPECTRUM/MonoPeak.h>
#include <GAGPL/GAGLIBRARY/SequencePrediction.h>

namespace gag
{
  using namespace std;

	typedef vector<ModificationSequence> SequenceList;

  // The code and corresponding scores.
	typedef map<std::string, vector<double>> SequenceScoreList;

	// This code doesn't specify the position of sulfation on the ring, but only the number of sulfate groups each ring.
	typedef SequenceCode MergedCode;

	class SequenceSpace
	{
	public:
		// Constructor
		SequenceSpace(GlycanSequencePtr gs_ptr, bool flag = false)
			:gs(gs_ptr), merged(flag)
		{
    }

		// This function should generate the list of modification patterns for given modification combination.
		SequenceList generateModificationCombination();

		// Conversion from sequence to sequence code.
		SequenceCode getSequenceCode(const ModificationSequence& mod_seq);
		// Conversion from sequence to merged sequence code.
		MergedCode getMergedCode(const ModificationSequence& mod_seq);
		// Conversion from sequence to string representation.
		std::string getSequenceString(const ModificationSequence& mod_seq);
		// Conversion from sequence to merged string representation.
		std::string getMergedCodeString(const ModificationSequence& mod_seq);
		// Conversion from sequence code to sequence.
		ModificationSequence getSequenceMap(const std::string& seq_code);

    // Check if the sequence fit the biosynthetic rule.
    // If 3-O is 1, then 2-N and 6-O have to be 1.
    // If 6-O is 1, then 2-N has to be 1.
    //bool isQualifiedSequence(const ModificationSequence& mod_seq);

    inline void setStatus(bool status)
    {
      merged = status;
    }

    // This function will continue generate new modification sequence
    // If there is no more sequence, return false.
    ModificationSequence nextModificationSequence();

    //std::multimap<MonoPeakPtr, NodeItem> filterData(std::multimap<MonoPeakPtr, NodeItem>& matched_results, const ModificationSequence& seq);

    // Generate the beginning mod_seq and initialize container_seq at the same time.
    // Generate mod_seq from prototype_seq.
    ModificationSequence getCurrentSequence() const;

	private:
		// Function for internal usage.
		void nextModificationSymbol(set<string>& mod_set, set<string>::iterator current_iter, SequenceList& seq_list, ModificationSequence mod_seq, Modifier mod_manager);

    // If true, new candidate sequence has been found.
    // If false, no new candidate sequence on the next level based on current configuration.
    bool nextModificationSymbol(set<string>& mod_set, std::set<string>::iterator mod_iter, ModificationSites occupied_sites);

    // This function will generate a pool of available candidate sites.
    void updateSites(set<string>& mod_set, std::set<string>::iterator mod_iter, ModificationSites prev_occupied_sites);

    void initializeModificationSequence();

		// Update the modification status table stored in gs.
		void applyeModificationStatus(const ModificationSequence& mod_seq);

	private:
		GlycanSequencePtr gs;

    // The container sequence storing the total available candidate sites.
    std::map<std::string, std::vector<ModificationPosition>> container_seq;

    // The prototype sequence storing the currently occupied candidate sites.
    std::map<std::string, std::vector<ModificationPosition>> prototype_seq;

    bool merged;
		
	};
}
#endif /* GAG_SEQUENCESPACE_H */