#include "MomentHelper.h"

void test_simple() {
    MomentHelper helper;
    helper.Set("H_0_0_0", 2.0);
    std::cout << "Value: " << helper.GetVal("H_0_0_0") << std::endl;
    std::cout << "Test completed successfully!" << std::endl;
}
