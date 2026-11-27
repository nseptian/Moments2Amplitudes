#include <TFile.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TDecompChol.h>
#include <TObjString.h>
#include <TList.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory>

// ============================================================================
// RegularizeCovarianceSPD:
// Regularizes a covariance matrix to be strictly Symmetric Positive-Definite (SPD)
// suitable for use as an MCMC proposal distribution.
//
// 1. Decouples Yield off-diagonals (isolating yield variance from physics moments)
// 2. Extracts absolute scale factors (sigma_i = sqrt(|diag_i|)) to repair negative signs
// 3. Projects correlation matrix R onto positive-definite cone via eigenvalue floor
// 4. Reconstructs SPD covariance matrix and verifies Cholesky decomposition
// ============================================================================
TMatrixDSym RegularizeCovarianceSPD(const TMatrixDSym& inMat, const std::vector<std::string>& paramNames, double minEigenFloor = 1e-4) {
    int n = inMat.GetNrows();
    TMatrixDSym mat = inMat;

    // First check if the raw matrix is already positive definite with all positive diagonals
    bool alreadySPD = true;
    for (int i = 0; i < n; ++i) {
        if (mat(i, i) <= 0.0) {
            alreadySPD = false;
            break;
        }
    }
    if (alreadySPD) {
        TDecompChol testChol(mat);
        if (testChol.Decompose()) {
            std::cout << "  [RegularizeCovarianceSPD] Matrix is already strictly positive-definite. No modification needed." << std::endl;
            return mat;
        }
    }

    std::cout << "  [RegularizeCovarianceSPD] Matrix is not positive-definite (has non-positive diagonals or eigenvalues)." << std::endl;
    std::cout << "  [RegularizeCovarianceSPD] Applying correlation-space spectral SPD projection..." << std::endl;

    // 1. Identify Yield parameter index
    int yieldIdx = -1;
    for (size_t i = 0; i < paramNames.size(); ++i) {
        TString name = paramNames[i].c_str();
        if (name.BeginsWith("Yld_") || name.Contains("Yield")) {
            yieldIdx = static_cast<int>(i);
            break;
        }
    }
    if (yieldIdx == -1) yieldIdx = n - 1; // Default fallback to last parameter

    // Decouple Yield off-diagonals and ensure positive yield variance
    double yieldVar = std::abs(mat(yieldIdx, yieldIdx));
    if (yieldVar <= 1.0) yieldVar = 10000.0;
    for (int i = 0; i < n; ++i) {
        if (i != yieldIdx) {
            mat(yieldIdx, i) = 0.0;
            mat(i, yieldIdx) = 0.0;
        }
    }
    mat(yieldIdx, yieldIdx) = yieldVar;

    // 2. Extract physics block variances: sigma_i = sqrt(|mat(i,i)|)
    std::vector<double> sigma(n, 1.0);
    for (int i = 0; i < n; ++i) {
        if (i == yieldIdx) {
            sigma[i] = std::sqrt(yieldVar);
            continue;
        }
        double d = std::abs(mat(i, i));
        // If variance is 0 or unphysically tiny, floor to reasonable moment scale (~0.05)
        sigma[i] = (d > 1e-8) ? std::sqrt(d) : 0.05;
    }

    // 3. Form Correlation Matrix: R_ij = mat_ij / (sigma_i * sigma_j)
    TMatrixDSym R(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == yieldIdx || j == yieldIdx) {
                R(i, j) = (i == j) ? 1.0 : 0.0;
            } else {
                double r = mat(i, j) / (sigma[i] * sigma[j]);
                R(i, j) = r;
            }
        }
    }
    for (int i = 0; i < n; ++i) R(i, i) = 1.0; // Diagonals of correlation matrix are strictly 1.0

    // 4. Spectral Projection on Correlation Matrix: R = V * Lambda * V^T
    TMatrixDSymEigen eigen(R);
    TVectorD eigenVals = eigen.GetEigenValues();
    TMatrixD eigenVecs = eigen.GetEigenVectors();

    int clampedCount = 0;
    for (int i = 0; i < n; ++i) {
        if (eigenVals(i) < minEigenFloor) {
            eigenVals(i) = minEigenFloor;
            clampedCount++;
        }
    }

    // Reconstruct SPD Correlation Matrix
    TMatrixDSym spdR(n);
    spdR.Zero();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += eigenVecs(i, k) * eigenVals(k) * eigenVecs(j, k);
            }
            spdR(i, j) = sum;
        }
    }

    // Renormalize diagonals of R to strictly 1.0
    for (int i = 0; i < n; ++i) {
        double d = std::sqrt(std::max(spdR(i, i), 1e-12));
        for (int j = 0; j < n; ++j) {
            spdR(i, j) /= d;
            spdR(j, i) = spdR(i, j);
        }
    }

    // 5. Rescale back to Covariance: Cov_ij = R_ij * sigma_i * sigma_j
    TMatrixDSym spdCov(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            spdCov(i, j) = spdR(i, j) * sigma[i] * sigma[j];
        }
    }

    // 6. Verify Cholesky decomposition
    TDecompChol chol(spdCov);
    bool cholOk = chol.Decompose();
    if (!cholOk) {
        std::cout << "  [RegularizeCovarianceSPD] Adding small diagonal jitter (1e-4) to guarantee Cholesky..." << std::endl;
        for (int i = 0; i < n; ++i) {
            spdCov(i, i) += sigma[i] * sigma[i] * 1e-4;
        }
        TDecompChol cholRetry(spdCov);
        cholOk = cholRetry.Decompose();
    }

    std::cout << "  [RegularizeCovarianceSPD] Alignment complete: " << clampedCount << "/" << n 
              << " eigenvalues regularized. Cholesky: " << (cholOk ? "PASSED (SPD Guaranteed!)" : "FAILED") << std::endl;

    return spdCov;
}

// Helper to save matrix and parameter names in BruFit-compatible format
bool SaveBruFitMatrix(const TString& outPath, const TString& matrixName, const TMatrixDSym& mat, const std::vector<std::string>& params) {
    std::unique_ptr<TFile> outFile(TFile::Open(outPath, "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outPath << std::endl;
        return false;
    }

    mat.Write(matrixName);

    TList nameList;
    nameList.SetName(matrixName + "_Names");
    for (const auto& name : params) {
        nameList.Add(new TObjString(name.c_str()));
    }
    nameList.Write(nameList.GetName(), TObject::kSingleKey);

    outFile->Close();
    return true;
}

void ExportCovarianceC(TString dir, TString matrixFilename, TString outputFilename) {
    // Open the ROOT file and read the covariance matrix
    std::cout << "Opening file: " << matrixFilename << std::endl;
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
        }
    }

    int N = paramNames.size();
    std::cout << "Loaded " << N << " parameter names." << std::endl;

    TMatrixDSym* myMatrixPointer = (TMatrixDSym*)file->Get("MyCovMatrix");
    if (!myMatrixPointer) {
        std::cerr << "Error: MyCovMatrix not found in file." << std::endl;
        file->Close();
        return;
    }

    std::cout << "Raw covariance matrix dimensions: " << myMatrixPointer->GetNrows() << "x" << myMatrixPointer->GetNcols() << std::endl;

    // Regularize the matrix into a strictly Symmetric Positive-Definite (SPD) matrix
    TMatrixDSym alignedMatrix = RegularizeCovarianceSPD(*myMatrixPointer, paramNames);

    TString outPath = dir;
    if (!outPath.EndsWith("/")) outPath += "/";
    outPath += outputFilename;

    if (SaveBruFitMatrix(outPath, "MyCovMatrix", alignedMatrix, paramNames)) {
        std::cout << "Successfully exported " << N << "x" << N << " aligned SPD covariance matrix to " << outPath << std::endl;
    } else {
        std::cerr << "Failed to export covariance matrix to " << outPath << std::endl;
    }

    file->Close();
}

// Directory-level batch export: regularizes all covariance_matrix_*.root files in dir
void ExportCovarianceC(TString dir) {
    if (!dir.EndsWith("/")) dir += "/";
    void* dirp = gSystem->OpenDirectory(dir);
    if (!dirp) {
        std::cerr << "Error opening directory: " << dir << std::endl;
        return;
    }
    const char* entry;
    std::vector<TString> files;
    while ((entry = gSystem->GetDirEntry(dirp))) {
        TString fn(entry);
        if (fn.BeginsWith("covariance_matrix_") && fn.EndsWith(".root")) {
            files.push_back(fn);
        }
    }
    gSystem->FreeDirectory(dirp);
    std::sort(files.begin(), files.end());

    std::cout << "ExportCovarianceC: Found " << files.size() << " covariance matrix file(s) in " << dir << std::endl;
    for (const auto& fn : files) {
        TString outFn = fn;
        outFn.ReplaceAll("covariance_matrix_", "BruCovMatrix_");
        ExportCovarianceC(dir, dir + fn, outFn);
    }
}
