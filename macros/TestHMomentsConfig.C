#include "MassDependentFitter.h"
#include <iostream>
#include <vector>
#include <iomanip>

void TestHMomentsConfig() {
    using namespace m2pw;
    
    std::cout << "=== Testing H Moments Configuration Interface ===" << std::endl;
    
    // Test the three specifically requested H moment configurations
    std::cout << "\n1. Testing L=4 Only Configuration:" << std::endl;
    auto l4OnlyConfig = MassDependentFitter::CreateL4OnlyConfig();
    
    for (int L = 0; L <= 6; ++L) {
        std::cout << "   - Include L=" << L << ": " << (l4OnlyConfig.ShouldIncludeL(L) ? "YES" : "NO") << std::endl;
    }
    
    std::cout << "\n2. Testing L=2,4 Only Configuration:" << std::endl;
    auto l2l4OnlyConfig = MassDependentFitter::CreateL2L4OnlyConfig();
    
    for (int L = 0; L <= 6; ++L) {
        std::cout << "   - Include L=" << L << ": " << (l2l4OnlyConfig.ShouldIncludeL(L) ? "YES" : "NO") << std::endl;
    }
    
    std::cout << "\n3. Testing Include All Configuration:" << std::endl;
    auto includeAllConfig = MassDependentFitter::CreateIncludeAllConfig();
    
    for (int L = 0; L <= 6; ++L) {
        std::cout << "   - Include L=" << L << ": " << (includeAllConfig.ShouldIncludeL(L) ? "YES" : "NO") << std::endl;
    }
    
    std::cout << "\n=== Configuration Summary Table ===" << std::endl;
    
    // Show how different configurations affect the number of included moments
    std::vector<std::pair<TString, MassDependentFitter::HMomentsConfig>> configs = {
        {"Include All", MassDependentFitter::CreateIncludeAllConfig()},
        {"L=4 Only", MassDependentFitter::CreateL4OnlyConfig()},
        {"L=2,4 Only", MassDependentFitter::CreateL2L4OnlyConfig()}
    };
    
    std::cout << "\n4. Configuration Summary:" << std::endl;
    std::cout << "   Config Name     | L=0 | L=1 | L=2 | L=3 | L=4 | L=5 | L=6 |" << std::endl;
    std::cout << "   --------------- | --- | --- | --- | --- | --- | --- | --- |" << std::endl;
    
    for (const auto& [name, config] : configs) {
        std::cout << "   " << std::left << std::setw(15) << name << " |";
        for (int L = 0; L <= 6; ++L) {
            std::cout << "  " << (config.ShouldIncludeL(L) ? "Y" : "N") << "  |";
        }
        std::cout << std::endl;
    }
    
    std::cout << "\n=== Testing H Moment Name Parsing ===" << std::endl;
    
    // Test equation name parsing (simplified version)
    std::vector<TString> testHMomentNames = {
        "H_0_0_0",      // L=0
        "H_0_2_0",      // L=2
        "H_1_2_1",      // L=2
        "H_0_4_0",      // L=4
        "H_1_4_2",      // L=4
        "H_0_6_4",      // L=6
        "normalise",    // Not an H moment
        "other_eqn"     // Not an H moment
    };
    
    std::cout << "\n5. Expected L values from H moment equation names:" << std::endl;
    for (const auto& name : testHMomentNames) {
        std::cout << "   " << std::setw(12) << name << " -> ";
        
        if (name.BeginsWith("H_")) {
            // Simple parsing logic (H_alpha_L_M format)
            std::unique_ptr<TObjArray> parts(name.Tokenize("_"));
            if (parts->GetEntries() >= 3) {
                TString L_str = ((TObjString*)parts->At(2))->GetString();
                if (L_str.IsDigit()) {
                    int L = L_str.Atoi();
                    std::cout << "L=" << L;
                    
                    // Test each configuration
                    std::cout << " (L4Only: " << (l4OnlyConfig.ShouldIncludeL(L) ? "Y" : "N");
                    std::cout << ", L2L4Only: " << (l2l4OnlyConfig.ShouldIncludeL(L) ? "Y" : "N");
                    std::cout << ", All: " << (includeAllConfig.ShouldIncludeL(L) ? "Y" : "N") << ")";
                } else {
                    std::cout << "Invalid L";
                }
            } else {
                std::cout << "Invalid format";
            }
        } else {
            std::cout << "Not H moment";
        }
        std::cout << std::endl;
    }
    
    std::cout << "\n=== Usage Examples ===" << std::endl;
    
    std::cout << "\n6. Practical Usage Examples:" << std::endl;
    
    std::cout << "\n   Example 1: Fit only L=4 H moments" << std::endl;
    std::cout << "   // This would exclude all H moments except those with L=4" << std::endl;
    std::cout << "   auto massConfig = MassDependentFitter::CreateDefaultConfig();" << std::endl;
    std::cout << "   auto hConfig = MassDependentFitter::CreateL4OnlyConfig();" << std::endl;
    std::cout << "   auto fitter = MassDependentFitter(setup, massBins, 2, 0.0, {}, massConfig, hConfig);" << std::endl;
    
    std::cout << "\n   Example 2: Fit only L=2,4 H moments" << std::endl;
    std::cout << "   // This would include only H moments with L=2 or L=4" << std::endl;
    std::cout << "   auto massConfig = MassDependentFitter::CreateMultipleA2Config();" << std::endl;
    std::cout << "   auto hConfig = MassDependentFitter::CreateL2L4OnlyConfig();" << std::endl;
    std::cout << "   auto fitter = MassDependentFitter(setup, massBins, 2, 0.0, {}, massConfig, hConfig);" << std::endl;
    
    std::cout << "\n   Example 3: Include all H moments (default)" << std::endl;
    std::cout << "   // This includes all H moments (same as not specifying hConfig)" << std::endl;
    std::cout << "   auto massConfig = MassDependentFitter::CreateL0L2Config();" << std::endl;
    std::cout << "   auto hConfig = MassDependentFitter::CreateIncludeAllConfig();" << std::endl;
    std::cout << "   auto fitter = MassDependentFitter(setup, massBins, 2, 0.0, {}, massConfig, hConfig);" << std::endl;
    
    std::cout << "\n=== Benefits of H Moment Selection ===" << std::endl;
    
    std::cout << "\n7. Why Select Specific H Moments:" << std::endl;
    std::cout << "   - Reduce computational complexity by excluding high-L moments" << std::endl;
    std::cout << "   - Focus on specific partial wave contributions" << std::endl;
    std::cout << "   - Improve fit stability by using only well-measured moments" << std::endl;
    std::cout << "   - Study systematic effects of different L ranges" << std::endl;
    std::cout << "   - Compare fits with different moment selections" << std::endl;
    
    std::cout << "\n=== Test Complete ===" << std::endl;
}