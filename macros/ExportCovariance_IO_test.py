
import ROOT
import os
import subprocess

import sys
# Moments_dir = os.getenv("Moments")
Moments_dir = "/d/grid17/septian/Moments-dev/Moments.Nizar/"
sys.path.append(Moments_dir)

# Moments_config_dir = os.getenv("MOM2AMPS")
Mom2Amp_dir = "/d/home/septian/Moments2Amplitudes/"
print("Using Moment2Amplitude directory:", Mom2Amp_dir)
Moments_config_dir = Mom2Amp_dir + "/Moments.config/" 
sys.path.append(Moments_config_dir)

import MomentCalculator
from MomentCalculator import (
  MomentResult,
  MomentResultsKinematicBinning,
  MomentValue,
  QnMomentIndex,
)

BRUFIT_DIR = "/d/grid17/septian/brufit_dev/brufit_cov/"
os.environ["BRUFIT"] = BRUFIT_DIR

# ROOT.gROOT.ProcessLine(f".x {BRUFIT_DIR}/macros/LoadBru.C")
# ROOT.gROOT.ProcessLine(f".L ExportCovarianceC.C")

# call ExportCovarianceC function for each maxL and binIndex

maxLs = (
   2,
   # 4,
   # 6,
   # 8,
)

base_dir = "/d/grid17/septian/Moments2Amplitudes/Moments.out/IO_test_m128132_Lmax4_ZeroPWaves/"
maxL_dir = {maxL: f"{base_dir}.maxL_{maxL}/" for maxL in maxLs}

print("maxL directories:")
for maxL, dir in maxL_dir.items():
   print(f"maxL = {maxL}: {dir}")

print("Loading moment results...")
momentResults_maxL: dict[int, MomentResultsKinematicBinning] = {
   maxL: MomentResultsKinematicBinning.loadPickle(f"{dir}/unnorm_moments_phys.pkl") for maxL, dir in maxL_dir.items()
}

for maxL, momentResults in momentResults_maxL.items():
   print(f"Loaded {len(momentResults)} moment results for maxL = {maxL} from {maxL_dir[maxL]}/unnorm_moments_phys.pkl")

   # Get number of moments for this maxL
   paramNames: list[str] = []
   binIndex, momentResult = next(enumerate(momentResults))
   print(f"maxL = {maxL}, binIndex = {binIndex}, momentResult length: {len(momentResult)}")

   counter = 0

   for momentValue in momentResult:
      qnIndex: QnMomentIndex = momentValue.qn
      momentIndex = qnIndex.momentIndex
      L = qnIndex.L
      M = qnIndex.M
      # label = momentValue.label
      # print(f"qnNumber: momentIndex = {momentIndex}, L = {L}, M = {M}")
      if (momentIndex == 0 and L == 0 and M == 0):
         continue  # skip H_0_0_0, which is not included in paramNames

      label = f"H_{momentIndex}_{L}_{M}"
      paramNames.append(label)
      # print(counter, label)
      counter += 1
      if counter == (len(momentResult) - 1):
         break

   paramNames.append("Yld_Moments")  # add Yld_Moments at the end

   # print(f"Appended paramNames for maxL = {maxL}, number of moments: {len(paramNames)}")

   nmbMoments: int = len(paramNames)

# Load the unnormalized moment results from the pickle file
# momentResults: MomentResultsKinematicBinning = MomentResultsKinematicBinning.loadPickle(f"{dir}/unnorm_moments_phys.pkl")

# print(f"Loaded {len(momentResults)} moment results from {dir}/unnorm_moments_phys.pkl")
# print(f"Moment results length: {len(momentResults)}")

# paramNames: list[str] = [
#    "H_0_1_0","H_0_1_1","H_0_2_0","H_0_2_1","H_0_2_2",
#    "H_0_3_0","H_0_3_1","H_0_3_2","H_0_3_3","H_0_4_0",
#    "H_0_4_1","H_0_4_2","H_0_4_3","H_0_4_4","H_1_0_0",
#    "H_1_1_0","H_1_1_1","H_1_2_0","H_1_2_1","H_1_2_2",
#    "H_1_3_0","H_1_3_1","H_1_3_2","H_1_3_3","H_1_4_0",
#    "H_1_4_1","H_1_4_2","H_1_4_3","H_1_4_4","H_2_1_1",
#    "H_2_2_1","H_2_2_2","H_2_3_1","H_2_3_2","H_2_3_3",
#    "H_2_4_1","H_2_4_2","H_2_4_3","H_2_4_4","Yld_Moments"
# ]

# nmbMoments: int = len(paramNames)

# print(f"Number of moments: {nmbMoments}")


   for binIndex, momentResult in enumerate(momentResults):
      # Get unnormalized H_0_0_0 in this bin
      H_0_0_0_qn = QnMomentIndex(momentIndex = 0, L = 0, M = 0)
      H_0_0_0_index: int = momentResult.indices[H_0_0_0_qn]
      H_0_0_0_value: MomentValue = momentResult[H_0_0_0_index]
      H_0_0_0_real, H_0_0_0_real_uncert = H_0_0_0_value.real

   #  print(momentResult.binCenters)

      momentResult.normalize()

      massCenter = next((center for var, center in momentResult.binCenters.items() if var.name == "Mpi0eta"), None)
      if massCenter is None:
         print(f"Bin {binIndex}: Re[H_0_0_0] = {H_0_0_0_real} +- {H_0_0_0_real_uncert}")
      else:
         print(f"Bin {binIndex} (Mpi0eta = {massCenter}): Re[H_0_0_0] = {H_0_0_0_real} +- {H_0_0_0_real_uncert}")

      myMatrix: ROOT.TMatrixDSym = ROOT.TMatrixDSym(nmbMoments)  # covariance matrix physical moments

      flatIndices: list[int] = []  # indices of moments in momentResult corresponding to the entry in paramNames
      useRealParts: list[bool] = []  # True for Re part, False for Im part for each entry in paramNames
      h1_norm_moments: dict[ROOT.TH1D] = {}
      for label in paramNames:
         if label.startswith("H_"):
            quantumNumbers: list[str] = label.split("_")  # naming pattern: H_<moment index>_<L>_<M>
            momentIndex = int(quantumNumbers[1])
            qnIndex: QnMomentIndex = QnMomentIndex(
               momentIndex = momentIndex,
               L           = int(quantumNumbers[2]),
               M           = int(quantumNumbers[3]),
            )
            flatIndices.append(momentResult.indices[qnIndex])
            # Convention: momentIndex 0 and 1 use Re parts, momentIndex 2 uses Im parts.
            useRealParts.append(momentIndex in (0, 1))
            if (useRealParts[-1]):
               moment_value, moment_error = momentResult[flatIndices[-1]].real
            else:
               moment_value, moment_error = momentResult[flatIndices[-1]].imag

            # update the mass binning
            h1_norm_moments[label] = ROOT.TH1D(f"h1_norm_{label}", label, 1, 1.28, 1.32)
            moment_value *= 2.0
            moment_error *= 2.0
            print(f"{label}: {moment_value} +- {moment_error}")
            h1_norm_moments[label].SetBinContent(1, moment_value)
            h1_norm_moments[label].SetBinError(1, moment_error)

      flatIndices.append(H_0_0_0_index)
      useRealParts.append(True)  # Yld_Moments is built from Re[H_0_0_0]

      scales = [H_0_0_0_real if name == "Yld_Moments" else 2.0 for name in paramNames]

      for i in range(nmbMoments):
         for j in range(i, nmbMoments):
            covEst = momentResult.covariance((flatIndices[i], flatIndices[j]), (useRealParts[i], useRealParts[j]))
            orig = covEst[0, 0] if i == j else covEst[0, 1]
            scaled = orig * scales[i] * scales[j]
            myMatrix[i][j] = scaled
            myMatrix[j][i] = scaled

      binLabel = MomentCalculator.binLabel(momentResult)
      dir = maxL_dir[maxL]
      savePath = f"{dir}/covariance_matrix_{binLabel}.root"
      myMatrix.Print()

      print(f"Saving the moments value and covariance matrix for bin {binIndex} to {savePath}")

      rootFile = ROOT.TFile(savePath, "RECREATE")
      myMatrix.Write("MyCovMatrix")

      for label in paramNames:
         if label != "Yld_Moments":
            h1_norm_moments[label].Write()

      # Save parameter names in each ROOT file for downstream lookup.
      paramNamesStr = ROOT.TObjString(",".join(paramNames))
      paramNamesStr.Write("paramNames")
         

      rootFile.Close()

      # Saving in brufit format with ExportCovarianceC
      
      # dir: ROOT.TString = maxL_dir[maxL]
      # savePath: ROOT.TString = f"{dir}/BruCovMatrix_{binLabel}.root"
      # print(f"Exporting covariance matrix for bin {binIndex} to {savePath} in brufit format")
      # ROOT.ExportCovarianceC(dir, savePath, "MyCovMatrix")
      
      outputFilename = f"BruCovMatrix_{binLabel}.root"
      bash_script = "/d/grid17/septian/Moments2Amplitudes/macros/ExportCovariance.sh"

      subprocess.run([
          "bash",
          bash_script,
          dir,
          savePath,
          outputFilename
      ])