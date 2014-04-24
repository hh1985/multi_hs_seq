#include "GAGPL/MATH/NewtonGirardFormulae.h"

namespace msmath
{
    // New:
    void NewtonGirardFormulae::updateElementarySymmetricPoly(const PowerSumVec& ps_vec, EleSymPolyVec& esp_vec)
    {
        // Based on the definition of ESP
        // k * e(k) = sum(j: 1 to k) [(-1)^(j-1) * e(k-j) * p(j)]
        // Notice: the size of esp_vec increases during the loop!
        size_t begin_k = esp_vec.size();
        size_t end_k = ps_vec.size();
        for(size_t k = begin_k; k < end_k; k++)
        {
            // If k = 0, e(k) = 1.
            // If k > n, e(k) = 0.
            if(k == 0) {
                esp_vec.push_back(1.0);
                continue;
            } else if(k > _var_num) {
                esp_vec.push_back(0.0);
                continue;
            }

            double temp_ele = 0.0;

            for(size_t j=1; j<=k; j++)
            {
                int sign = (j%2 == 1 ? 1 : -1);

                temp_ele += sign * ps_vec.at(j) * esp_vec.at(k-j); 
            }

            temp_ele /= k;

            esp_vec.push_back(temp_ele);
        }

    }
    // New:
    void NewtonGirardFormulae::updatePowerSum(PowerSumVec& ps_vec, const EleSymPolyVec& esp_vec)
    {
        size_t end_k = esp_vec.size();
        size_t begin_k = ps_vec.size();
        // p(k) = sum(j: 1 to k-1) [(-1)^(j-1) * e(j) * p(k-j)] + (-1)^k * e(k) * k.
        for(size_t k= begin_k; k<end_k; k++)
        {
            if(k == 0) {
                ps_vec.push_back(0);
                continue;
            }

            double temp_ps = 0.0;

            // From 1 to k-1
            int sign = -1;
            for(size_t j=1; j<k; j++)
            {
                sign *= -1;
                temp_ps += sign * esp_vec.at(j) * ps_vec.at(k-j);
            }

            sign *= -1;
            temp_ps += sign * esp_vec.at(k) * k;

            ps_vec.push_back(temp_ps);
        }

    }
    // New:
    void NewtonGirardFormulae::updateParameters( PowerSumVec& ps_vec, EleSymPolyVec& esp_vec )
    {
        // The default size for esp_vec should be _var_num + 1.
        if(ps_vec.size() > esp_vec.size())
            this->updateElementarySymmetricPoly(ps_vec, esp_vec);
        else if(ps_vec.size() < esp_vec.size())
            this->updatePowerSum(ps_vec, esp_vec);
        else
            ; // Do nothing.
    }
}
