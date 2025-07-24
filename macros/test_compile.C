void test_compile() {
    cout << "Testing compilation..." << endl;
    
    // Test if the headers can be included and basic classes work
    try {
        std::vector<TString> testNames = {"test1", "test2"};
        std::map<TString, double> testMap;
        testMap["test"] = 1.0;
        
        cout << "Basic STL containers work" << endl;
        cout << "Test vector size: " << testNames.size() << endl;
        cout << "Test map value: " << testMap["test"] << endl;
        cout << "Headers included successfully" << endl;
    } catch (...) {
        cout << "Error in headers" << endl;
    }
    
    cout << "Compilation test complete" << endl;
}
