#include <GAGPL/GAGLIBRARY/Assignment.h>

namespace gag
{
    ostream& operator<<(ostream& os, const Assignment& assignment)
    {
        os << "Fragment:" << "\n";
        os << *(assignment.fg);

        os << "MZ:" << assignment.pk->mz << "\t" << "Z:" << assignment.pk->z << "\n";
        for(auto iter = assignment.neu_loss.begin(); iter != assignment.neu_loss.end(); iter++)
        {
            os << iter->first << ":" << iter->second << "\t";
        }

        for(auto iter = assignment.mod_count.begin(); iter != assignment.mod_count.end(); iter++)
        {
            os << iter->first << ":" << iter->second << "\t";
        }
        os << "\n";

        return os;
    }


    bool Assignment::isNRECleavage()
    {
      string frag_type = fg->getFragmentType();
      if(frag_type == "A" || frag_type == "B" || frag_type == "C")
        return true;
      else
        return false;
    }

    bool Assignment::isRECleavage()
    {
      string frag_type = fg->getFragmentType();
      if(frag_type == "X" || frag_type == "Y" || frag_type == "Z")
        return true;
      else
        return false;
    }

}