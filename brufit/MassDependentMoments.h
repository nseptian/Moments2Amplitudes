#pragma once

#include <MomentHelper.h>

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
                throw runtime_error("Mass value not found in moments map.");
            }
        }

        void PrintMoments() {
            cout << "Mass-dependent moments size: " << _moments.size() << endl;
            for (auto it = _moments.begin(); it != _moments.end(); ++it) {
                cout << endl;
                cout << "Mass: " << it->first << endl;
                it->second.PrintVals();
            }
        }

    private:
        map<Double_t, MomentHelper> _moments;

};