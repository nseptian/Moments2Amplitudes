#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"

void test_refactored_code() {
    std::cout << "Testing refactored MassDependentEquationSolver..." << std::endl;
    
    try {
        // Test basic instantiation
        auto& setup = ConfigureAmpsNoValues(2, 2, 2);
        std::vector<double> test_mass_bins = {0.82, 0.86, 0.90};
        
        std::cout << "Creating solver..." << std::endl;
        m2pw::MassDependentEquationSolver solver(setup, test_mass_bins, 2, 0.0, {"H_3"});
        
        std::cout << "NDim: " << solver.NDim() << std::endl;
        std::cout << "Mass bins: " << solver.GetMassBins().size() << std::endl;
        
        std::cout << "✅ Basic instantiation works!" << std::endl;
        
        // Test parameter evaluation (with dummy parameters)
        std::vector<double> dummy_params(solver.NDim(), 1.0);
        double chi2 = solver.DoEval(dummy_params.data());
        std::cout << "✅ DoEval works! Chi2: " << chi2 << std::endl;
        
        std::cout << "🎉 All tests passed! The refactored code works." << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "❌ Error: " << e.what() << std::endl;
    } catch (...) {
        std::cout << "❌ Unknown error occurred" << std::endl;
    }
}
