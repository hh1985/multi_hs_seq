/********************************************************************
	created:	2014/01/15
	created:	15:1:2014   22:54
	filename: 	NaiveMethods.h
	file path:	GAGPL\GAGLIBRARY
	file base:	NaiveMethods
	file ext:	h
	author:		Han Hu
	
	purpose:	Methods for comparison with HS-SEQ.
*********************************************************************/

#ifndef GAG_NAIVEMETHODS_H
#define GAG_NAIVEMETHODS_H

#include <GAGPL/GAGLIBRARY/GeneralCleavage.h>
#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/SPECTRUM/MonoPeak.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/GAGLIBRARY/CleavageTree.h>

namespace gag
{
  
  struct ModificationKey
  {
    ModificationKey() : key_size(0) {}
    
    size_t key_size;
    std::string mod_key;
    //std::vector<std::pair<std::string, int>> mod_key;

    // Convert the information into string.
    void addKey(const ModificationSites& mod_sites, int mod_num)
    {
      mod_key.append(modificationString(mod_sites));
      mod_key.append(":");
      mod_key.append(boost::lexical_cast<std::string>(mod_num));
      mod_key.append("~");
      //mod_key.push_back(std::make_pair(modificationString(mod_sites), mod_num));
      key_size++;
    }

    bool operator<(const ModificationKey& rkey) const
    {
      // This struct requires that the two keys should have the identical size.
      if(key_size != rkey.key_size)
        throw std::runtime_error("Invalid comparison between the modification keys!");

      return mod_key < rkey.mod_key;
    }

    void print()
    {
      std::cout << mod_key << "\n";
    }
  };

  typedef std::set<ModificationKey> ModificationKeySet;

  class NaiveMethods
  {
  public:
    // Constructor.
    NaiveMethods(std::multimap<MonoPeakPtr, NodeItem>& pk_list, GlycanSequencePtr seq, int flag = 1)
      : gs(seq), data(pk_list), param(Param::Instance()), status(flag) {}

    inline void setStatus(int flag)
    {
      status = flag;
    }

    // Naive method 1 for calculating the coverage of mass values for each candidate sequence.
    double calculateCoverage(const ModificationSequence& mod_seq, int cleavage_num = 2);

    // Naive method 2 for calculating the number of golden pair. Notice to remove redundant pairs.  Basically, if the modification sites of two assignments are complementary to each other, we define this as a pair.
    double calculateGoldenPairNum(const ModificationSequence& mod_seq);

    std::multimap<MonoPeakPtr, NodeItem> filterData(const ModificationSequence& mod_seq);


  private:
    GlycanSequencePtr gs;
    //const ModificationSequence& _mod_seq;
    std::multimap<MonoPeakPtr, NodeItem>& data;
    Param& param;
    int status;
  };
}

#endif /* GAG_NAIVEMETHODS_H */