/********************************************************************
	created:	2012/06/26
	created:	26:6:2012   15:53
	filename: 	Modification.cpp
	file path:	GAG\src\GAGPL\CHEMISTRY
	file base:	Modification
	file ext:	cpp
	author:		Han Hu
	
	purpose: Modification class should be simple enough to be independent
	         of any system libraries.
*********************************************************************/

#include <GAGPL/CHEMISTRY/Modification.h>
#include <iostream>

namespace gag
{
	//ModificationPosition::mod_id = 0;

	std::vector<size_t> Modification::getModificationSites( const std::string& fg )
	{
		std::map<std::string, std::vector<size_t> >::iterator iter = mod_rule.find(fg);
		if(iter != mod_rule.end())
			return iter->second;
		else
			return std::vector<size_t>();
	}

	std::string ModificationPosition::printString() const
	{
		std::string mod_string;
		mod_string.append(boost::lexical_cast<std::string>(macro_pos.branch_id));
		mod_string.append("-");
		mod_string.append(boost::lexical_cast<std::string>(macro_pos.mono_id));
		mod_string.append("-");
		mod_string.append(boost::lexical_cast<std::string>(site_id));
		return mod_string;
	}

	ModificationSites getSiteDifference( const ModificationSites& m1, const ModificationSites& m2 )
	{
		ModificationSites diff_sites;
		std::set_difference(m1.begin(), m1.end(), m2.begin(), m2.end(), std::inserter(diff_sites, diff_sites.end()));
		return diff_sites;
	}

	void printModificationSites( const ModificationSites& mod_sites )
	{
    if(mod_sites.size() == 0) {
      std::cout << "Empty!\n";
      return;
    }

		for(ModificationSites::const_iterator iter = mod_sites.begin(); 
			iter != mod_sites.end(); iter++)
		{
			std::cout << iter->printString() << " ";
		}
    std::cout << "\n";
	}

	ModificationSites getSiteIntersection( const ModificationSites& m1, const ModificationSites& m2 )
	{
		ModificationSites inter_sites;
		std::set_intersection(m1.begin(), m1.end(), m2.begin(), m2.end(), std::inserter(inter_sites, inter_sites.end()));
		return inter_sites;
	}

  std::string modificationString( const ModificationSites& mod_sites )
  {
    std::string mod_str;
    if(mod_sites.size() != 0) {
      //std::cout << "Empty!\n";
      for(ModificationSites::const_iterator iter = mod_sites.begin(); 
        iter != mod_sites.end(); iter++)
      {
        mod_str.append(iter->printString());
        mod_str.append("/");
      }
    }

    return mod_str;
  }
  // A includes B
  bool containSubset( const ModificationSites& ms1, const ModificationSites& ms2 )
  {
    return std::includes(ms1.begin(), ms1.end(), ms2.begin(), ms2.end());
  }

  std::ostream& operator<<(std::ostream& os, const ModificationSites& mod_sites)
  {
      if(mod_sites.size() == 0) {
          os << "Empty!";
          return os;
      }

      for(ModificationSites::const_iterator iter = mod_sites.begin(); iter != mod_sites.end(); iter++)
      {
          os << iter->printString() << " ";
      }
      //os << "\n";

      return os;
  }

}
