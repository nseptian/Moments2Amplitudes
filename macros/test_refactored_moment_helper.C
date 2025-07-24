#include "MomentHelper.h"
#include <iostream>

void test_refactored_moment_helper() {
    std::cout << "Testing refactored MomentHelper class..." << std::endl;
    
    // Test basic functionality
    MomentHelper helper;
    
    // Test direct setting
    helper.Set("H_1_0_0", 1.5);
    helper.Set("H_0_1_0", 0.8);
    helper.Set("H_0_0_1", 0.3);
    
    std::cout << "Direct setting test:" << std::endl;
    std::cout << "H_1_0_0 = " << helper.GetVal("H_1_0_0") << std::endl;
    std::cout << "H_0_1_0 = " << helper.GetVal("H_0_1_0") << std::endl;
    std::cout << "H_0_0_1 = " << helper.GetVal("H_0_0_1") << std::endl;
    std::cout << "Non-existent = " << helper.GetVal("H_999_999_999") << std::endl;
    
    // Test PrintVals
    std::cout << "\nAll moments:" << std::endl;
    helper.PrintVals();
    
    std::cout << "\nMoments with prefix 'H_1':" << std::endl;
    helper.PrintVals("H_1");
    
    std::cout << "\nRefactored MomentHelper test completed successfully!" << std::endl;
}
