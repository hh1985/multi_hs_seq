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

    set<Backbone> AssignmentPool::selectQualifiedAssignments( const string& mod_symbol )
    {
      AssignmentsByMass& mass_index = _pool.get<theo_mass>();
      map<ModificationSites, Backbone> grouped_assignments;
      for(auto mass_iter = mass_index.begin(); mass_iter != mass_index.end(); mass_iter++)
      {
        ModificationSites mod_sites = mass_iter->getBackboneModificationSites(mod_symbol);
        auto bone_iter = grouped_assignments.find(mod_sites);
        if(bone_iter == grouped_assignments.end()) {
          Backbone single_bone(*mass_iter, mod_symbol);
          grouped_assignments.insert(make_pair(mod_sites, single_bone));
        } else {
          bone_iter->second.addAssignment(*mass_iter, mod_symbol);
        }
      }

      set<Backbone> bone_set;
      for(auto iter = grouped_assignments.begin(); iter != grouped_assignments.end(); iter++)
      {
        bone_set.insert(iter->second);
      }

      return bone_set;

    }

}
