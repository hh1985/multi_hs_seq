#include <GAGPL/GAGLIBRARY/AssignmentPool.h>

namespace gag
{


    void AssignmentPool::addAssignment(AssignmentPtr ptr)
    {
        _pool.insert(ptr);
    }

    set<AssignmentPtr> AssignmentPool::getAssignmentsByMass(double mass)
    {
        AssignmentsByMass& mass_index = _pool.get<theo_mass>();
        auto p = mass_index.equal_range(mass);
        set<AssignmentPtr> assign_set;
        while(p.first != p.second)
        {
            assign_set.insert(*(p.first));
        }
        return assign_set;
    }

    set<AssignmentPtr> AssignmentPool::getAssignmentsByModificationSites(const std::string& mod, const ModificationSites& mod_sites)
    {
        set<AssignmentPtr> assign_set;
        if(mod == "Ac") {
            AssignmentsByACSites& ac_index = _pool.get<ac_sites>();
            auto p = ac_index.equal_range(mod_sites);
            while(p.first != p.second)
            {
                assign_set.insert(*(p.first));
            }
        } else if(mod == "SO3") {
            AssignmentsBySulfateSites& sulfate_index = _pool.get<sulfate_sites>();
            auto p = sulfate_index.equal_range(mod_sites);
            while(p.first != p.second)
            {
                assign_set.insert(*(p.first));
            }
        }
        return assign_set;
    }

    map<int, set<AssignmentPtr>> AssignmentPool::sortAssignmentsByUniqueness(ModificationSites&)
    {
        map<double, set<AssignmentPtr>> mass_map;
        AssignmentsByMass& mass_index = _pool.get<theo_mass>();
        for(auto iter = mass_index.begin(); iter != mass_index.end(); iter++)
        {
            mass_map[iter->getMass()].insert(*iter);
        }

        map<int, set<AssignmentPtr>> uni_map;

        for(auto mass_iter = mass_map.begin(); mass_iter != mass_map.end(); mass_iter++)
        {
            auto uni_pair = this->calculateUniquenessValue(mass_iter->second);
            uni_map.insert(uni_pair);
        }
        
        return uni_map;
    }

    pair<int, set<AssignmentPtr>> AssignmentPool::calculateUniquenessValue(set<AssignmentPtr> assign_set)
    {
        // Separate the assignments by their cleavage type.
        map<string, set<AssignmentPtr>> type_map;

        for(auto iter = assign_set.begin(); iter != assign_set.end(); iter++)
        {
            type_map[iter->getCleavageType].insert(*iter);
        }

        

    }



}
