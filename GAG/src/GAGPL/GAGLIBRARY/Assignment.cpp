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


}