#!/d/grid17/septian/Moments.Nizar/.pixi/envs/default/bin/python
"""
This module performs the moment analysis of unpolarized and polarized
pi+ pi- photoproduction data. The calculated moments are written to
files to be read by the plotting script `photoProdPlotMoments.py`.

Usage: Run this module as a script to perform the moment calculations
and to generate the output files.
"""


from __future__ import annotations

from collections.abc import Iterable
from copy import deepcopy
import functools
import multiprocessing
import numpy as np
import os
import threadpoolctl

import ROOT
from wurlitzer import pipes, STDOUT

import sys
# Moments_dir = os.getenv("Moments")
Moments_dir = "/d/grid17/septian/Moments-dev/Moments.Nizar/"
sys.path.append(Moments_dir)

# Moments_config_dir = os.getenv("MOM2AMPS")
Mom2Amp_dir = "/d/home/septian/Moments2Amplitudes/"
print("Using Moment2Amplitude directory:", Mom2Amp_dir)
Moments_config_dir = Mom2Amp_dir + "/Moments.config/" 
sys.path.append(Moments_config_dir)

from Config_IO_test import (
  AnalysisConfig,
  BEAM_POL_INFOS,
  CFG_IO_TEST_LMAX4_ZEROPWAVES,
  CFG_IO_TEST_LMAX4_WITHPWAVES
)
from MomentCalculator import (
  DataSet,
  MomentCalculator,
  MomentCalculatorsKinematicBinning,
)
from PlottingUtilities import HistAxisBinning
import RootUtilities  # importing initializes OpenMP and loads `basisFunctions.C`
import Utilities


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def calculateAllMoments(
  cfg:                            AnalysisConfig,
  timer:                          Utilities.Timer        = Utilities.Timer(),
  forceIntegralMatrixCalculation: bool                   = True,  # if `True` integral matrices are recalculated even if pickled versions exist
  limitToDataEntryRange:          tuple[int, int] | None = None,  # for debugging: limit analysis to entry range [begin, end) of real-data data tree
  additionalCuts:                 Iterable[str] | None   = None,  # optional additional cuts to be applied to real and accepted phase-space data
  additionalColumnDefs:           dict[str, str]         = {},    # additional columns to define
) -> None:
  """Performs the moment analysis for the given configuration"""
  # setup MomentCalculators for all data samples and mass bins
  momentCalculators: dict[str | None, MomentCalculatorsKinematicBinning] = {  # dict: key = data sample label
    None : MomentCalculatorsKinematicBinning([]),
  } if cfg.method == AnalysisConfig.MethodType.LIN_ALG_BG_SUBTR_NEG_WEIGHTS else {
    "Sig" : MomentCalculatorsKinematicBinning([]),
    "Bkg" : MomentCalculatorsKinematicBinning([]),
  }
  assert len(cfg.massBinning) > 0, f"Need at least one mass bin, but found {len(cfg.massBinning)}"
  with timer.timeThis(f"Time to load data and setup MomentCalculators for {len(cfg.massBinning)} bins"):
    dataSamples: dict[str | None, ROOT.RDataFrame] = {  # dict: key = data sample label
      None : cfg.loadData(AnalysisConfig.DataType.REAL_DATA, additionalCuts, additionalColumnDefs),
    } if cfg.method == AnalysisConfig.MethodType.LIN_ALG_BG_SUBTR_NEG_WEIGHTS else {
      "Sig" : cfg.loadData(AnalysisConfig.DataType.REAL_DATA_SIGNAL, additionalCuts, additionalColumnDefs),    # select events in signal region (mixture of signal and background)
      "Bkg" : cfg.loadData(AnalysisConfig.DataType.REAL_DATA_SIDEBAND, additionalCuts, additionalColumnDefs),  # select events in sideband regions (ideally, pure background)
    }
    #TODO: implement event range limit
    # if limitToDataEntryRange is not None:
    #   print(f"Limiting analysis to entry range [{limitToDataEntryRange[0]}, {limitToDataEntryRange[1]}) of real data")
    #   data = data.Range(*limitToDataEntryRange)
    # use same MC events for all real-data samples
    dataPsAcc = cfg.loadData(AnalysisConfig.DataType.ACCEPTED_PHASE_SPACE, additionalCuts, additionalColumnDefs)
    if cfg.limitNmbPsAccEvents > 0 and dataPsAcc is not None:
      dataPsAcc = dataPsAcc.Range(cfg.limitNmbPsAccEvents)  #!Caution! Range() switches to single-threaded mode
    dataPsGen = cfg.loadData(AnalysisConfig.DataType.GENERATED_PHASE_SPACE, additionalColumnDefs = additionalColumnDefs)
    for labelDataSample, dataSample in dataSamples.items():
      if dataSample is None:
        print(f"No real-data events for type '{labelDataSample}'. Skipping.")
        continue
      print(f"Loaded {dataSample.Count().GetValue()} real-data events" + (f" of type '{labelDataSample}'" if labelDataSample is not None else ""))
      for massBinIndex, massBinCenter in enumerate(cfg.massBinning):
        massBinRange = cfg.massBinning.binValueRange(massBinIndex)
        print(f"Preparing {cfg.binVarMass.name} bin [{massBinIndex + 1} of {len(cfg.massBinning)}] at {massBinCenter} {cfg.binVarMass.unit} with range {massBinRange} {cfg.binVarMass.unit}")
        massBinFilter = cfg.massBinning.binFilter(massBinIndex)
        print(f"Applying filter '{massBinFilter}' to select kinematic bin")
        dataInBin = dataSample.Filter(massBinFilter)
        print(f"Loaded {dataInBin.Count().GetValue()} real-data events for kinematic bin")
        if "weightASBS" in dataInBin.GetColumnNames():
          print(f"Real-data sample contains {dataInBin.Sum('weightASBS').GetValue()} events after applying event weights")
        dataPsAccInBin = None if dataPsAcc is None else dataPsAcc.Filter(massBinFilter)
        nmbPsAccEvents = None if dataPsAcc is None else dataPsAccInBin.Count().GetValue()
        nmbPsGenEvents = None
        if nmbPsAccEvents is not None:
          if dataPsGen is not None:
            dataPsGenInBin = dataPsGen.Filter(massBinFilter)
            nmbPsGenEvents = dataPsGenInBin.Count().GetValue()
            print(f"Loaded phase-space events: number generated = {nmbPsGenEvents}; "
                  f"number accepted = {nmbPsAccEvents} "
                  f"-> efficiency = {nmbPsAccEvents / nmbPsGenEvents:.3f}")
          else:
            print(f"Loaded phase-space events: number accepted = {nmbPsAccEvents}")
        # setup moment calculators for data
        dataSet = DataSet(
          data           = dataInBin,
          phaseSpaceData = dataPsAccInBin,
          nmbGenEvents   = nmbPsGenEvents or nmbPsAccEvents or 0,
          polarization   = cfg.polarization,
        )
        momentCalculators[labelDataSample].append(
          MomentCalculator(
            indicesMeas          = cfg.momentIndicesMeas,
            indicesPhys          = cfg.momentIndicesPhys,
            dataSet              = dataSet,
            binCenters           = {cfg.binVarMass : massBinCenter},
            integralFileBaseName = f"{cfg.outFileDirName}/integralMatrix",
            flipSignOfWeights    = (labelDataSample == "Bkg"),  # flip sign of weights for background events
          )
        )

  if cfg.method in {AnalysisConfig.MethodType.LIN_ALG_BG_SUBTR_NEG_WEIGHTS, AnalysisConfig.MethodType.LIN_ALG_BG_SUBTR_MOMENTS}:
    # calculate integral matrix for all mass bins in all data sets
    nmbOpenMpThreads = ROOT.getNmbOpenMpThreads()
    # since all data sets use identical MC events, calculate integral matrix only for the first real-data sample
    labelFirstSample, momentCalculatorsFirstSample = next(iter(momentCalculators.items()))  # get first entry in momentCalculators
    with timer.timeThis(f"Time to calculate integral matrices for {len(momentCalculatorsFirstSample)} bins using {nmbOpenMpThreads} OpenMP threads"):
      print(f"Calculating acceptance integral matrices for {len(momentCalculatorsFirstSample)} bins using {nmbOpenMpThreads} OpenMP threads")
      momentCalculatorsFirstSample.calculateIntegralMatrices(forceIntegralMatrixCalculation)
      print(f"Acceptance integral matrix for first bin at {cfg.massBinning[0]} {cfg.binVarMass.unit}:\n{momentCalculatorsFirstSample[0].integralMatrix}")
      # eigenVals, _ = momentCalculatorsFirstSample[0].integralMatrix.eigenDecomp
      # print(f"Sorted eigenvalues of acceptance integral matrix for first bin at {cfg.massBinning[0]} {cfg.binVarMass.unit}:\n{np.sort(eigenVals)}")
    # assign integral matrices to `MomentCalculators` for other data samples
    labelsOtherSamples = (label for label in momentCalculators.keys() if label is not labelFirstSample)
    for labelDataSample in labelsOtherSamples:
      for massBinIndex in range(len(cfg.massBinning)):
        if len(momentCalculators[labelDataSample]) > 0:  # if data sample exists
          momentCalculators[labelDataSample][massBinIndex]._integralMatrix = momentCalculatorsFirstSample[massBinIndex].integralMatrix

  momentResultsFileBaseName = f"{cfg.outFileDirName}/{cfg.outFileNamePrefix}_moments"
  if cfg.calcAccPsMoments:
    # calculate moments of accepted phase-space data
    if cfg.method in {AnalysisConfig.MethodType.LIN_ALG_BG_SUBTR_NEG_WEIGHTS, AnalysisConfig.MethodType.LIN_ALG_BG_SUBTR_MOMENTS}:
      with timer.timeThis(f"Time to calculate moments of phase-space MC data using the linear-algebra method with {nmbOpenMpThreads} OpenMP threads"):
        # since all data samples use identical MC events, calculate moments only for first accepted phase-space sample
        momentCalculatorsFirstSample = next(iter(momentCalculators.values()))  # get first entry in momentCalculators
        print(f"Calculating moments of phase-space MC data for {len(momentCalculatorsFirstSample)} bins using {nmbOpenMpThreads} OpenMP threads")
        momentCalculatorsFirstSample.calculateMoments(
          dataSourceType = MomentCalculator.DataSourceType.ACCEPTED_PHASE_SPACE,
          normalize      = cfg.normalizeMoments,
        )
        momentCalculatorsFirstSample.momentResultsMeas.savePickle(f"{momentResultsFileBaseName}_accPs_meas.pkl")
        momentCalculatorsFirstSample.momentResultsPhys.savePickle(f"{momentResultsFileBaseName}_accPs_phys.pkl")
    else:
      #TODO calculate also using FIT method
      raise NotImplementedError(f"Calculation of accepted phase-space moments using method {cfg.method} is not implemented yet")

  if cfg.method == AnalysisConfig.MethodType.LIN_ALG_BG_SUBTR_NEG_WEIGHTS:
    #TODO calculate normalized and unnormalized moments
    with timer.timeThis(f"Time to calculate moments of real data for {len(momentCalculators[None])} bins using the linear-algebra method with {nmbOpenMpThreads} OpenMP threads and subtracting background using negative weights"):
      print(f"Calculating moments of real data for {len(momentCalculators[None])} bins using the linear-algebra method with {nmbOpenMpThreads} OpenMP threads and subtracting background using negative weights")
      momentCalculators[None].calculateMoments(
        normalize           = cfg.normalizeMoments,
        nmbBootstrapSamples = cfg.nmbBootstrapSamples,
      )
      momentCalculators[None].momentResultsMeas.savePickle(f"{momentResultsFileBaseName}_meas.pkl")
      momentCalculators[None].momentResultsPhys.savePickle(f"{momentResultsFileBaseName}_phys.pkl")

  elif cfg.method == AnalysisConfig.MethodType.LIN_ALG_BG_SUBTR_MOMENTS:
    for labelDataSample, momentCalculatorsDataSample in momentCalculators.items():
      with timer.timeThis(f"Time to calculate moments of '{labelDataSample}' real data with {len(momentCalculatorsDataSample)} bins using the linear-algebra method with {nmbOpenMpThreads} OpenMP threads and subtracting background at moment level"):
        print(f"Calculating moments of '{labelDataSample}' real data with {len(momentCalculatorsDataSample)} bins using the linear-algebra method with {nmbOpenMpThreads} OpenMP threads and subtracting background at moment level")
        momentCalculatorsDataSample.calculateMoments(
          normalize           = cfg.normalizeMoments,
          nmbBootstrapSamples = cfg.nmbBootstrapSamples,
        )
        momentCalculatorsDataSample.momentResultsMeas.savePickle(f"{momentResultsFileBaseName}_meas.pkl")
        momentCalculatorsDataSample.momentResultsPhys.savePickle(f"{momentResultsFileBaseName}_phys.pkl")
    if len(momentCalculators["Bkg"]) > 0:
      # perform background subtraction at moment level
      #!NOTE! This gives wrong estimates for the uncertainty (and covariances) of
      # H_0(0, 0) because the number of signal and background events are treated
      # as being constant. The uncertainty estimates for the other moments seem
      # to be unaffected.
      (momentCalculators["Sig"].momentResultsMeas - momentCalculators["Bkg"].momentResultsMeas).savePickle(f"{momentResultsFileBaseName}_meas.pkl")
      (momentCalculators["Sig"].momentResultsPhys - momentCalculators["Bkg"].momentResultsPhys).savePickle(f"{momentResultsFileBaseName}_phys.pkl")
    else:
      # background sample was empty
      momentCalculators["Sig"].momentResultsMeas.savePickle(f"{momentResultsFileBaseName}_meas.pkl")
      momentCalculators["Sig"].momentResultsPhys.savePickle(f"{momentResultsFileBaseName}_phys.pkl")

  elif cfg.method == AnalysisConfig.MethodType.MAX_LIKELIHOOD_FIT:
    # nmbFitAttempts          = 100
    # nmbParallelFitProcesses = 100
    nmbFitAttempts          = 10
    nmbParallelFitProcesses = 10
    startValueRandomSeed    = 123456789
    for labelDataSample, momentCalculatorsDataSample in momentCalculators.items():
      if len(momentCalculatorsDataSample) > 0:
        with timer.timeThis(f"Time to fit moments to '{labelDataSample}' real data with {len(momentCalculatorsDataSample)} bins running {nmbFitAttempts} fit attempts in {nmbParallelFitProcesses} processes and subtracting background at moment level"):
          print(f"Fitting moments to '{labelDataSample}' real data with {len(momentCalculatorsDataSample)} bins running {nmbFitAttempts} fit attempts in {nmbParallelFitProcesses} processes and subtracting background at moment level")
          momentCalculatorsDataSample.fitMomentsMultipleAttempts(
            nmbFitAttempts          = nmbFitAttempts,
            nmbParallelFitProcesses = nmbParallelFitProcesses,
            startValueRandomSeed    = startValueRandomSeed,
          )
    if len(momentCalculators["Bkg"]) > 0:
      # perform background subtraction at moment level
      (momentCalculators["Sig"].momentResultsPhys - momentCalculators["Bkg"].momentResultsPhys).savePickle(f"{momentResultsFileBaseName}_phys.pkl")
    else:
      # background sample was empty
      momentCalculators["Sig"].momentResultsPhys.savePickle(f"{momentResultsFileBaseName}_phys.pkl")

  else:
    raise ValueError(f"Unknown method {cfg.method}")


def _runAnalysisConfigurationInSeparateProcess(
  cfg:                            AnalysisConfig,
  forceIntegralMatrixCalculation: bool,
  dataDirBaseName:                str,
  subsystemLabel:                 str,
  outFileDirBaseNameCommon:       str,
  dataPeriod:                     str,
  tBinLabel:                      str,
  beamPolLabel:                   str,
  maxL:                           int,
) -> None:
  cfg = deepcopy(cfg)
  cfg.dataFileName       = f"{dataDirBaseName}/{tBinLabel}/{cfg.dataFileName}"
  cfg.psAccFileName      = f"{dataDirBaseName}/{tBinLabel}/{cfg.psAccFileName}"
  cfg.psGenFileName      = f"{dataDirBaseName}/{tBinLabel}/{cfg.psGenFileName}"
  # cfg.outFileDirBaseName = f"{outFileDirBaseNameCommon}/{tBinLabel}/{beamPolLabel}"
  # cfg.polarization       = "Pol" if BEAM_POL_INFOS[dataPeriod[:7]][beamPolLabel] is not None else None  #TODO why is this needed? shouldn't this be set in AnalysisConfig already?
  cfg.maxL = maxL
  cfg.init(createOutFileDir = True)
  thisSourceFileName = os.path.basename(__file__)
  logFileName = f"{cfg.outFileDirName}/{os.path.splitext(thisSourceFileName)[0]}_{cfg.outFileNamePrefix}.log"
  print(f"Writing output to log file '{logFileName}'")
  with open(logFileName, "w") as logFile, pipes(stdout = logFile, stderr = STDOUT):  # redirect all output into log file
    Utilities.printGitInfo()
    ROOT.gROOT.SetBatch(True)
    timer = Utilities.Timer()
    timer.start("Total execution time")
    threadController = threadpoolctl.ThreadpoolController()  # at this point all multi-threading libraries must be loaded
    print(f"Initial state of ThreadpoolController before setting number of threads:\n{threadController.info()}")
    with threadController.limit(limits = 4):
      print(f"State of ThreadpoolController after setting number of threads:\n{threadController.info()}")
      print(f"Using configuration:\n{cfg}")
      calculateAllMoments(
        cfg                            = cfg,
        timer                          = timer,
        forceIntegralMatrixCalculation = forceIntegralMatrixCalculation,
        additionalCuts                 = None,
        additionalColumnDefs           = {},
      )
    timer.stop("Total execution time")
    print(timer.summary)


# CPP_CODE_IS_IN_EFFICIENCY_HOLES = """
#   bool
#   isInEfficiencyHoles(
#     const double cosTheta,
#     const double phi
#   ) {
#     static std::vector<TEllipse*> holes = {
#       new TEllipse(-0.55,    0, 0.35, 30),
#       new TEllipse( 0.55, -180, 0.35, 30),
#       new TEllipse( 0.55, +180, 0.35, 30)
#     };
#     for (const auto& hole : holes) {
#       if (hole->IsInside(cosTheta, phi)) {
#         return true;
#       }
#     }
#     return false;
#   }
# """

CPP_CODE_IS_IN_EFFICIENCY_HOLES = """
  bool
  isInEfficiencyHoles(
    const double cosTheta,
    const double phi
  ) {
    // Define hole parameters: {x_center, y_center, x_radius, y_radius}
    static const double holes[3][4] = {
      {-0.55,    0.0, 0.35, 30.0},
      { 0.55, -180.0, 0.35, 30.0},
      { 0.55,  180.0, 0.35, 30.0}
    };

    for (int i = 0; i < 3; ++i) {
      double dx = cosTheta - holes[i][0];
      double dy = phi - holes[i][1];
      double rx = holes[i][2];
      double ry = holes[i][3];

      // Ellipse equation: (dx^2 / rx^2) + (dy^2 / ry^2) <= 1.0
      if (((dx * dx) / (rx * rx)) + ((dy * dy) / (ry * ry)) <= 1.0) {
        return true;
      }
    }
    return false;
  }
"""


if __name__ == "__main__":
  cfg = deepcopy(CFG_IO_TEST_LMAX4_WITHPWAVES)  # perform analysis of Nizar's polarized eta pi0 data
  subsystemLabel = "EtaPi0"
  os.chdir(Mom2Amp_dir)  # Change to the Moments directory
  dataDirBaseName = cfg.dataDirBaseName
  print("Using data directory:", dataDirBaseName)
  dataPeriods = ("merged", )
  tBinLabels = (
    "t010020",
  )
  beamPolLabels = ("All", )
  maxLs = (
    4,
    6,
    8,
  )
  forceIntegralMatrixCalculation = True  # if `True` integral matrices are recalculated even if pickled versions exist
  # forceIntegralMatrixCalculation = False

  ROOT.gInterpreter.Declare(CPP_CODE_IS_IN_EFFICIENCY_HOLES)

  outFileDirBaseNameCommon = cfg.outFileDirBaseName
  if os.path.isdir(outFileDirBaseNameCommon):
    print(f"Using output directory: {outFileDirBaseNameCommon}")
  else:
    print(f"Creating output directory: {outFileDirBaseNameCommon}")
    os.makedirs(outFileDirBaseNameCommon)

  
  workerContext = multiprocessing.get_context("fork")
  workerProcesses: list[multiprocessing.Process] = []
  for dataPeriod in dataPeriods:
    for tBinLabel in tBinLabels:
      for beamPolLabel in beamPolLabels:
        for maxL in maxLs:
          print(f"Performing moment analysis for data period '{dataPeriod}', t bin '{tBinLabel}', beam-polarization orientation '{beamPolLabel}', and L_max = {maxL}")
          workerProcess = workerContext.Process(
            target = _runAnalysisConfigurationInSeparateProcess,
            args = (
              cfg,
              forceIntegralMatrixCalculation,
              dataDirBaseName,
              subsystemLabel,
              outFileDirBaseNameCommon,
              dataPeriod,
              tBinLabel,
              beamPolLabel,
              maxL,
            ),
          )
          workerProcess.start()
          workerProcesses.append(workerProcess)

  for workerProcess in workerProcesses:
    workerProcess.join()
    if workerProcess.exitcode != 0:
      raise RuntimeError(f"analysis subprocess exited with code {workerProcess.exitcode}")
