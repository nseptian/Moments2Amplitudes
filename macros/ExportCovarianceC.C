#include <TMatrixDSym.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

void ExportCovarianceC(TString dir, TString matrixFilename, TString outputFilename) {
    // int N = 15;
    // TMatrixDSym myMatrix(N);

    // TString dir = "/d/home/septian/Moments/plotsPhotoProdEtaPi0_m080240/merged/t075100/All.maxL_4/";
    // TString matrixFilename = dir + "covariance_matrix_Mpi0eta_1.300.root";

    // Open the ROOT file and read the covariance matrix
    cout << "Opening file: " << matrixFilename << std::endl;
    TFile* file = TFile::Open(matrixFilename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << matrixFilename << std::endl;
        return;
    }

    TObjString* paramNamesObj = dynamic_cast<TObjString*>(file->Get("paramNames"));
    if (!paramNamesObj) {
        std::cerr << "Error: paramNames not found in file." << std::endl;
        file->Close();
        return;
    }

    std::cout << "Parameter names from the file:" << std::endl;
    std::vector<std::string> paramNames;
    std::stringstream namesStream(paramNamesObj->GetString().Data());
    std::string paramName;
    while (std::getline(namesStream, paramName, ',')) {
        if (!paramName.empty()) {
            paramNames.push_back(paramName);
            std::cout << paramNames.back() << std::endl;
        }
    }

    int N = paramNames.size();

    TMatrixDSym* myMatrixPointer = (TMatrixDSym*)file->Get("MyCovMatrix");
    if (!myMatrixPointer) {
        std::cerr << "Error: MyCovMatrix not found in file." << std::endl;
        file->Close();
        return;
    }

    // Print the covariance matrix for verification
    std::cout << "Covariance matrix from the file:" << std::endl;
    myMatrixPointer->Print();

    TMatrixDSym& myMatrix = *myMatrixPointer;

    // TString outputFilename = "brufit_covariance_matrix_Mpi0eta_1.300.root";
    HS::FIT::BruCovarianceReader::SaveThis(dir+outputFilename, "MyCovMatrix", myMatrix, paramNames);
    
    std::cout << "Successfully exported " << N << "x" << N << " covariance matrix to " << outputFilename << std::endl;
}
