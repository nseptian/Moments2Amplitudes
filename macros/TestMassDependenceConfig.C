#include "MassDependentFitter.h"
#include <iostream>
#include <vector>
#include <iomanip>

void TestMassDependenceConfig() {
    using namespace m2pw;
    
    std::cout << "=== Testing MassDependenceConfig Interface ===" << std::endl;
    
    // Test 1: Default configuration
    std::cout << "\n1. Testing Default Configuration:" << std::endl;
    auto defaultConfig = MassDependentFitter::CreateDefaultConfig();
    
    std::cout << "   - L=0 mass dependent: " << (defaultConfig.IsMassDependent(0) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=1 mass dependent: " << (defaultConfig.IsMassDependent(1) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=2 mass dependent: " << (defaultConfig.IsMassDependent(2) ? "YES" : "NO") << std::endl;
    
    std::cout << "   - L=0 mass independent: " << (defaultConfig.IsMassIndependent(0) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=1 mass independent: " << (defaultConfig.IsMassIndependent(1) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=2 mass independent: " << (defaultConfig.IsMassIndependent(2) ? "YES" : "NO") << std::endl;
    
    auto l2_waves = defaultConfig.GetWavesForL(2);
    std::cout << "   - Waves for L=2: ";
    for (const auto& wave : l2_waves) {
        std::cout << wave << " ";
    }
    std::cout << std::endl;
    
    // Test 2: Multiple a2 configuration
    std::cout << "\n2. Testing Multiple a2 Configuration:" << std::endl;
    auto multiA2Config = MassDependentFitter::CreateMultipleA2Config();
    
    auto l2_waves_multi = multiA2Config.GetWavesForL(2);
    std::cout << "   - Waves for L=2: ";
    for (const auto& wave : l2_waves_multi) {
        std::cout << wave << " ";
    }
    std::cout << "(" << l2_waves_multi.size() << " waves)" << std::endl;
    
    // Test 3: L=0 + L=2 configuration
    std::cout << "\n3. Testing L=0 + L=2 Configuration:" << std::endl;
    auto l0l2Config = MassDependentFitter::CreateL0L2Config();
    
    std::cout << "   - L=0 mass dependent: " << (l0l2Config.IsMassDependent(0) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=1 mass dependent: " << (l0l2Config.IsMassDependent(1) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=2 mass dependent: " << (l0l2Config.IsMassDependent(2) ? "YES" : "NO") << std::endl;
    
    auto l0_waves = l0l2Config.GetWavesForL(0);
    std::cout << "   - Waves for L=0: ";
    for (const auto& wave : l0_waves) {
        std::cout << wave << " ";
    }
    std::cout << std::endl;
    
    auto l2_waves_l0l2 = l0l2Config.GetWavesForL(2);
    std::cout << "   - Waves for L=2: ";
    for (const auto& wave : l2_waves_l0l2) {
        std::cout << wave << " ";
    }
    std::cout << std::endl;
    
    // Test 4: Custom configuration
    std::cout << "\n4. Testing Custom Configuration:" << std::endl;
    std::map<int, std::vector<TString>> customMassDep;
    customMassDep[0] = {"a0_980"};
    customMassDep[1] = {};
    customMassDep[2] = {"a2_1320", "a2_1700"};
    
    std::vector<int> customMassIndep = {}; // All are mass dependent
    
    auto customConfig = MassDependentFitter::CreateCustomConfig(customMassDep, customMassIndep);
    
    for (int l = 0; l <= 2; ++l) {
        std::cout << "   - L=" << l << " mass dependent: " << (customConfig.IsMassDependent(l) ? "YES" : "NO");
        if (customConfig.IsMassDependent(l)) {
            auto waves = customConfig.GetWavesForL(l);
            std::cout << " (waves: ";
            for (const auto& wave : waves) {
                std::cout << wave << " ";
            }
            std::cout << ")";
        }
        std::cout << std::endl;
    }
    
    std::cout << "\n=== Testing ParameterManager with Configurations ===" << std::endl;
    
    // Test 5: ParameterManager with different configurations
    std::vector<double> testMassBins = {1.0, 1.1, 1.2, 1.3, 1.4};
    std::vector<TString> testParNames = {
        "a_0_0", "b_0_0",
        "a_1_-1", "a_1_0", "a_1_1", 
        "b_1_-1", "b_1_0", "b_1_1",
        "a_2_-2", "a_2_-1", "a_2_0", "a_2_1", "a_2_2",
        "b_2_-2", "b_2_-1", "b_2_0", "b_2_1", "b_2_2"
    };
    
    std::cout << "\n5. Testing ParameterManager with Default Config:" << std::endl;
    MassDependentFitter::ParameterManager pmDefault;
    pmDefault.AddMassIndependentParameters(testMassBins, testParNames, 42, defaultConfig);
    pmDefault.AddMassDependentParameters(testParNames, 42, defaultConfig);
    
    std::cout << "   - Total parameters (default config): " << pmDefault.totalNpars << std::endl;
    
    std::cout << "\n6. Testing ParameterManager with Multiple a2 Config:" << std::endl;
    MassDependentFitter::ParameterManager pmMultiA2;
    pmMultiA2.AddMassIndependentParameters(testMassBins, testParNames, 42, multiA2Config);
    pmMultiA2.AddMassDependentParameters(testParNames, 42, multiA2Config);
    
    std::cout << "   - Total parameters (multi-a2 config): " << pmMultiA2.totalNpars << std::endl;
    
    std::cout << "\n7. Testing ParameterManager with L=0+L=2 Config:" << std::endl;
    MassDependentFitter::ParameterManager pmL0L2;
    pmL0L2.AddMassIndependentParameters(testMassBins, testParNames, 42, l0l2Config);
    pmL0L2.AddMassDependentParameters(testParNames, 42, l0l2Config);
    
    std::cout << "   - Total parameters (L=0+L=2 config): " << pmL0L2.totalNpars << std::endl;
    
    std::cout << "\n8. Testing ParameterManager with Custom Config:" << std::endl;
    MassDependentFitter::ParameterManager pmCustom;
    pmCustom.AddMassIndependentParameters(testMassBins, testParNames, 42, customConfig);
    pmCustom.AddMassDependentParameters(testParNames, 42, customConfig);
    
    std::cout << "   - Total parameters (custom config): " << pmCustom.totalNpars << std::endl;
    
    // Test 6: Show some parameter names for the custom configuration
    std::cout << "\n9. Sample parameter names from custom config:" << std::endl;
    for (int i = 0; i < std::min(10, pmCustom.totalNpars); ++i) {
        std::cout << "   [" << i << "] " << pmCustom.parIndexNames[i] << std::endl;
    }
    if (pmCustom.totalNpars > 10) {
        std::cout << "   ... (" << (pmCustom.totalNpars - 10) << " more parameters)" << std::endl;
    }
    
    std::cout << "\n=== Testing H Moments Factory Methods ===" << std::endl;
    
    // Test the three specifically requested H moment configurations
    std::cout << "\n10. Testing H Moments L=4 Only Configuration:" << std::endl;
    auto hL4OnlyConfig = MassDependentFitter::CreateL4OnlyConfig();
    std::cout << "   - H moments included: ";
    for (int L = 0; L <= 6; ++L) {
        if (hL4OnlyConfig.ShouldIncludeL(L)) {
            std::cout << "L=" << L << " ";
        }
    }
    std::cout << std::endl;
    
    std::cout << "\n11. Testing H Moments L=2,4 Only Configuration:" << std::endl;
    auto hL2L4OnlyConfig = MassDependentFitter::CreateL2L4OnlyConfig();
    std::cout << "   - H moments included: ";
    for (int L = 0; L <= 6; ++L) {
        if (hL2L4OnlyConfig.ShouldIncludeL(L)) {
            std::cout << "L=" << L << " ";
        }
    }
    std::cout << std::endl;
    
    std::cout << "\n12. Testing H Moments Include All Configuration:" << std::endl;
    auto hIncludeAllConfig = MassDependentFitter::CreateIncludeAllConfig();
    std::cout << "   - H moments included: ";
    for (int L = 0; L <= 6; ++L) {
        if (hIncludeAllConfig.ShouldIncludeL(L)) {
            std::cout << "L=" << L << " ";
        }
    }
    std::cout << std::endl;
    
    // Show detailed breakdown for each configuration
    std::cout << "\n13. Detailed H Moment Configuration Breakdown:" << std::endl;
    std::cout << "   Config Name     | L=0 | L=1 | L=2 | L=3 | L=4 | L=5 | L=6 |" << std::endl;
    std::cout << "   --------------- | --- | --- | --- | --- | --- | --- | --- |" << std::endl;
    
    std::vector<std::pair<TString, MassDependentFitter::HMomentsConfig>> hConfigs = {
        {"Include All", MassDependentFitter::CreateIncludeAllConfig()},
        {"L=4 Only", MassDependentFitter::CreateL4OnlyConfig()},
        {"L=2,4 Only", MassDependentFitter::CreateL2L4OnlyConfig()}
    };
    
    for (const auto& [name, config] : hConfigs) {
        std::cout << "   " << std::left << std::setw(15) << name << " |";
        for (int L = 0; L <= 6; ++L) {
            std::cout << "  " << (config.ShouldIncludeL(L) ? "Y" : "N") << "  |";
        }
        std::cout << std::endl;
    }
    
    // Example usage
    std::cout << "\n14. Usage Examples:" << std::endl;
    std::cout << "   Example 1: Fit only L=4 moments" << std::endl;
    std::cout << "   auto hConfig = MassDependentFitter::CreateL4OnlyConfig();" << std::endl;
    std::cout << "   auto fitter = MassDependentFitter(setup, massBins, 2, 0.0, {}, massDepConfig, hConfig);" << std::endl;
    
    std::cout << "\n   Example 2: Fit only L=2,4 moments" << std::endl;
    std::cout << "   auto hConfig = MassDependentFitter::CreateL2L4OnlyConfig();" << std::endl;
    std::cout << "   auto fitter = MassDependentFitter(setup, massBins, 2, 0.0, {}, massDepConfig, hConfig);" << std::endl;
    
    std::cout << "\n   Example 3: Include all H moments (default behavior)" << std::endl;
    std::cout << "   auto hConfig = MassDependentFitter::CreateIncludeAllConfig();" << std::endl;
    std::cout << "   auto fitter = MassDependentFitter(setup, massBins, 2, 0.0, {}, massDepConfig, hConfig);" << std::endl;
    
    std::cout << "\n=== Test Complete ===" << std::endl;
}
