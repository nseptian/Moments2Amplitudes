#pragma once

#include "MomentHelper.h"
#include <map>
#include <stdexcept>

class MassDependentMoments {
    public:
        MassDependentMoments() = default;

        // Function to set the mass-dependent moments
        void SetMoments(const Double_t massValue, MomentHelper& moments) {
            _moments[massValue] = moments;
        }

        // Function to get the mass-dependent moments
        MomentHelper GetMoments(const Double_t massValue) const {
            auto it = _moments.find(massValue);
            if (it != _moments.end()) {
                return it->second;
            } else {
                throw std::runtime_error("Mass value not found in moments map.");
            }
        }

        void PrintMoments() {
            std::cout << "Mass-dependent moments size: " << _moments.size() << std::endl;
            for (auto it = _moments.begin(); it != _moments.end(); ++it) {
                std::cout << std::endl;
                std::cout << "Mass: " << it->first << std::endl;
                it->second.PrintVals();
            }
        }
        
        size_t GetNMassBins() {
            return _moments.size();
        }

    private:
        std::map<Double_t, MomentHelper> _moments;
};