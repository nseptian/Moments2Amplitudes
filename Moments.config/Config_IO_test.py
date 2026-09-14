"""
This module defines the input parameters for the various moment analyses
"""

from __future__ import annotations

from collections.abc import (
  Callable,
  Iterable,
)
from dataclasses import (
  dataclass,
  field,
)
from enum import Enum

import ROOT
import os
import sys

# Moments_dir = os.getenv("Moments")
Moments_dir = "/d/grid17/septian/Moments-dev/Moments.Nizar/"
sys.path.append(Moments_dir)

Mom2Amp_dir = "/d/home/septian/Moments2Amplitudes/"

from MomentCalculator import (
  KinematicBinningVariable,
  MomentIndices,
)
from PlottingUtilities import HistAxisBinning
import Utilities


def defineOverwriteRDataFrame(
  df:         ROOT.RDataFrame,
  colName:    str,
  newFormula: str,
) -> ROOT.RDataFrame:
  """Defines column `colName` using formula `newFormula`; overwriting column definition if it already exists"""
  if not df.HasColumn(colName):
    # print(f"Defining column '{colName}' = '{newFormula}'")
    df = df.Define(colName, newFormula)
  else:
    print(f"Redefining column '{colName}' to '{newFormula}'")
    df = df.Redefine(colName, newFormula)
  return df


@dataclass
class BeamPolInfo:
  """Stores information about beam polarization for a specific orientation"""
  pol:    float | str  # photon-beam polarization magnitude value or column name
  PhiLab: float | str  # azimuthal angle of photon beam polarization in lab frame [deg] value or column name

  def __str__(self) -> str:
    result = "pol = "
    if isinstance(self.pol, float):
      result += f"{self.pol:.4f}"
    else:
      result += f"'{self.pol}'"
    result += ", PhiLab = "
    if isinstance(self.PhiLab, float):
      result += f"{self.PhiLab:.1f} deg"
    else:
      result += f"'{self.PhiLab}'"
    return result


#TODO maybe using a function would be more flexible here
BEAM_POL_INFOS: dict[str, dict[str, BeamPolInfo | None]] = {  # <data-period label> : {<beam-polarization label> : BeamPolInfo(...)}; `None` means unpolarized
  "merged" : {  # several merged data periods with different polarization values
    "All" : BeamPolInfo(  # read polarization values from the given column names
      pol    = "beamPol",
      PhiLab = "beamPolPhiLabDeg",
    ),
    # approximate polarization values across all data-taking periods; don't use for serious analyses
    "PARA_0" : BeamPolInfo(
      pol    = 0.35,
      PhiLab = 0,
    ),
    "PERP_45" : BeamPolInfo(
      pol    = 0.35,
      PhiLab = 45,
    ),
    "PERP_90" : BeamPolInfo(
      pol    = 0.35,
      PhiLab = 90,
    ),
    "PARA_135" : BeamPolInfo(
      pol    = 0.35,
      PhiLab = -45,
    ),
    "AMO"   : None,
    "Unpol" : None,
  },
  # polarization values from Version 9 of `makePolVals` tool from https://halldweb.jlab.org/wiki-private/index.php/TPOL_Polarization
  # beam polarization angles in lab frame taken from `Lab Phi` column of tables 2 to 5 in GlueX-doc-3977
  "2017_01" : {  # polarization magnitudes obtained by running `.x makePolVals.C(17, 1, 0, 75)` in ROOT shell
    "PARA_0" : BeamPolInfo(
      pol    = 0.3537,
      PhiLab = 1.8,
    ),
    "PERP_45" : BeamPolInfo(
      pol    = 0.3484,
      PhiLab = 47.9,
    ),
    "PERP_90" : BeamPolInfo(
      pol    = 0.3472,
      PhiLab = 94.5,
    ),
    "PARA_135" : BeamPolInfo(
      pol    = 0.3512,
      PhiLab = -41.6,
    ),
    "AMO"   : None,
    "Unpol" : None,
  },
  "2018_01" : {  # polarization magnitudes obtained by running `.x makePolVals.C(18, 1, 0, 75)` in ROOT shell
    "PARA_0" : BeamPolInfo(
      pol    = 0.3420,
      PhiLab = 4.1,
    ),
    "PERP_45" : BeamPolInfo(
      pol    = 0.3474,
      PhiLab = 48.5,
    ),
    "PERP_90" : BeamPolInfo(
      pol    = 0.3478,
      PhiLab = 94.2,
    ),
    "PARA_135" : BeamPolInfo(
      pol    = 0.3517,
      PhiLab = -42.4,
    ),
    "AMO"   : None,
    "Unpol" : None,
  },
  "2018_08" : {  # polarization magnitudes obtained by running `.x makePolVals.C(18, 2, 0, 75)` in ROOT shell
    "PARA_0" : BeamPolInfo(
      pol    = 0.3563,
      PhiLab = 3.3,
    ),
    "PERP_45" : BeamPolInfo(
      pol    = 0.3403,
      PhiLab = 48.3,
    ),
    "PERP_90" : BeamPolInfo(
      pol    = 0.3430,
      PhiLab = 92.9,
    ),
    "PARA_135" : BeamPolInfo(
      pol    = 0.3523,
      PhiLab = -42.1,
    ),
    "AMO"   : None,
    "Unpol" : None,
  },
  "2019_11" : {  # Spring 2020 polarization magnitudes obtained by running `root makePolVals2020.C'(0)'`  #TODO did Alex already produce PhiLab values?
    "PARA_0" : BeamPolInfo(
      pol    = 0.3525,
      PhiLab = 0,
    ),
    "PERP_45" : BeamPolInfo(
      pol    = 0.3535,
      PhiLab = 45,
    ),
    "PERP_90" : BeamPolInfo(
      pol    = 0.3536,
      PhiLab = 90,
    ),
    "PARA_135" : BeamPolInfo(
      pol    = 0.3721,
      PhiLab = -45,
    ),
    "AMO"   : None,
    "Unpol" : None,
  },
}


@dataclass
class SubsystemInfo:
  """Stores information about the (A, B) two-body subsystem, for which moments are calculated, for the reaction beam + target -> A + B + recoil"""
  lvALabel:          str  # label of Lorentz-vector of daughter A (analyzer)
  lvBLabel:          str  # label of Lorentz-vector of daughter B
  lvRecoilLabel:     str  # label of Lorentz-vector of recoil particle
  pairLabel:         str  # label for particle pair (e.g. "PiPi" for pi+ pi- pair)
  ATLatexLabel:      str = ""  # optional LaTeX label for particle A (analyzer)
  BTLatexLabel:      str = ""  # optional LaTeX label for particle B
  recoilTLatexLabel: str = ""  # optional LaTeX label for recoil particle
  pairTLatexLabel:   str = ""  # optional LaTeX label for particle pair (e.g. "#pi#pi" for pi+ pi- pair)


@dataclass
class AnalysisConfig:
  """Stores configuration parameters for the moment analysis; defaults are for unpolarized pi+pi- analysis"""

  class MethodType(Enum):
    """Enumerates methods used for moment analysis"""
    LIN_ALG_BG_SUBTR_NEG_WEIGHTS = 0  # linear-algebra method for acceptance correction + background subtraction by negative event weights
    LIN_ALG_BG_SUBTR_MOMENTS     = 1  # linear-algebra method for acceptance correction + background subtraction by subtracting moments
    MAX_LIKELIHOOD_FIT           = 2  # maximum-likelihood fit for acceptance correction + background subtraction by subtracting moments

  class CoordSysType(Enum):
    """Enumerates coordinate systems in which moments can be calculated"""
    HF = 0  # helicity frame
    GJ = 1  # Gottfried-Jackson frame

  class DataType(Enum):
    """Enumerates used input data types"""
    REAL_DATA             = 0
    REAL_DATA_SIGNAL      = 1  # real data in signal region
    REAL_DATA_SIDEBAND    = 2  # real data in sideband region(s)
    GENERATED_PHASE_SPACE = 3
    ACCEPTED_PHASE_SPACE  = 4

  class DataFormat(Enum):
    """Enumerates formats of input data files"""
    ALEX            = 0  # Alex' data format  #TODO improve naming
    AMPTOOLS        = 1  # AmpTools format
    JPAC_MC         = 2  # MC truth data in JPAC text format
    TLORENTZVECTORS = 3  # TLorentzVector for each particle
    FSROOT          = 4  # FSROOT format

  # defaults are for unpolarized gamma p -> (pi+ pi-) p data
  method:                   AnalysisConfig.MethodType  = MethodType.LIN_ALG_BG_SUBTR_NEG_WEIGHTS  # method used to estimate moments from data
  frame:                    CoordSysType               = CoordSysType.HF  # coordinate system, in which moments are calculated
  subsystem:                SubsystemInfo              = SubsystemInfo(  # pi+pi- subsystem to analyze; pi+ is the analyzer
    lvALabel          = "pip",     # label of pi+ Lorentz-vector
    lvBLabel          = "pim",     # label of pi- Lorentz-vector
    lvRecoilLabel     = "recoil",  # label of recoil-proton Lorentz-vector
    pairLabel         = "PiPi",    #TODO treeName is identical
    ATLatexLabel      = "#it{#pi}^{#plus}",
    BTLatexLabel      = "#it{#pi}^{#minus}",
    recoilTLatexLabel = "#it{p}",
    pairTLatexLabel   = "#it{#pi}^{#plus}#it{#pi}^{#minus}",
  )
  dataDirBaseName:          str                        = "./dataPhotoProdPiPi/unpolarized"  # base directory for input data  #TODO turn into property and determine value depending on other member variables
  dataPeriods:              tuple[str, ...]            = (  # labels of data periods to process
    "2017_01",
    "2018_08",
  )
  tBinLabels:               tuple[str, ...]            = (  # labels of t bins to process
    "tbin_0.4_0.5",
  )
  beamPolLabels:            tuple[str, ...]            = (  # labels of beam polarizations to process, e.g. "PARA_0"; use "AMO" or "Unpol" for unpolarized data  #TODO store list of BeamPolInfos and generate labels with a function
    "Unpol",
  )
  inputDataFormats:         dict[DataType, DataFormat] = field(default_factory = lambda: {  # data formats for each input data type
    AnalysisConfig.DataType.REAL_DATA             : AnalysisConfig.DataFormat.ALEX,
    AnalysisConfig.DataType.ACCEPTED_PHASE_SPACE  : AnalysisConfig.DataFormat.ALEX,
    AnalysisConfig.DataType.GENERATED_PHASE_SPACE : AnalysisConfig.DataFormat.AMPTOOLS,
  })
  inputTreeName:            str                        = "kin"  # name of tree to read from input data and MC files

  treeName:                 str                        = "PiPi"  # name of tree to read from converted data and MC files  #TODO fix name

  # Spring 2017 low-energy data
  dataFileName:             str                        = "./dataPhotoProdPiPi/unpolarized/2017_01/tbin_0.4_0.5/PiPi/data_flat.root"  # file with real data to analyze  #TODO remove?
  psAccFileName:            str | None                 = "./dataPhotoProdPiPi/unpolarized/2017_01/tbin_0.4_0.5/PiPi/phaseSpace_acc_flat.root"  # file with accepted phase-space MC  #TODO remove?
  psGenFileName:            str | None                 = "./dataPhotoProdPiPi/unpolarized/2017_01/tbin_0.4_0.5/PiPi/phaseSpace_gen_flat.root"  # file with generated phase-space MC  #TODO remove?
  polarization:             float | str | None         = None  # photon-beam polarization; None = unpolarized photoproduction; polarized photoproduction: either polarization value or name of polarization column  #TODO use BeamPolInfo
  maxL:                     int | tuple[int, int]      = 8  # if int: maximum L of physical and measured moments; if tuple: (max L of physical moments, max L of measured moments)
  outFileDirBaseName:       str                        = "./plotsPhotoProdPiPiUnpolCLAS"  # base name of directory into which all output will be written
  # normalizeMoments:         bool                       = True
  normalizeMoments:         bool                       = False
  nmbBootstrapSamples:      int                        = 0
  # nmbBootstrapSamples:      int                        = 10000
  # plotAngularDistributions: bool                       = True
  plotAngularDistributions: bool                       = False
  # plotAccIntegralMatrices:  bool                       = True
  plotAccIntegralMatrices:  bool                       = False
  # calcAccPsMoments:         bool                       = True
  calcAccPsMoments:         bool                       = False
  # plotAccPsMoments:         bool                       = True
  plotAccPsMoments:         bool                       = False
  # plotMomentsInBins:        bool                       = True
  plotMomentsInBins:        bool                       = False
  # plotMeasuredMoments:      bool                       = True
  plotMeasuredMoments:      bool                       = False
  plotCovarianceMatrices:   bool                       = True
  # plotCovarianceMatrices:   bool                       = False
  limitNmbPsAccEvents:      int                        = 0
  # limitNmbPsAccEvents:      int                        = 100000
  binVarMass:               KinematicBinningVariable   = field(
    default_factory = lambda: KinematicBinningVariable(
      name      = "mass",
      label     = "#it{m}_{#it{#pi}^{#plus}#it{#pi}^{#minus}}",
      unit      = "GeV/#it{c}^{2}",
      nmbDigits = 3,
    )
  )
  massBinning:              HistAxisBinning            = field(
    default_factory = lambda: HistAxisBinning(  # same binning as used by CLAS
      nmbBins = 100,
      minVal  = 0.4,   # [GeV]
      maxVal  = 1.4,   # [GeV]
      _var    = None,  # set by init()  #TODO why keep a separate binVarMass member instead of using binning.var?
    )
  )

  def inputDataDirBasePath(
    self,
    dataPeriod: str,
    tBinLabel:  str,
  ) -> str:
    """Generates base path of directory with input data"""
    return f"{self.dataDirBaseName}/{dataPeriod}/{tBinLabel}/input"

  @staticmethod
  def _default_inputFilePaths(
    cfg:          AnalysisConfig,
    dataType:     AnalysisConfig.DataType,
    dataPeriod:   str,
    tBinLabel:    str,
    beamPolLabel: str
  ) -> tuple[str, ...]:
    """Default function that generates input file paths based on data type, data period, t bin label, and beam polarization label; can be overwritten by providing a different function to the `_inputFilePaths` member"""
    # one input file for each data type
    if dataType == AnalysisConfig.DataType.REAL_DATA:
      return (f"{cfg.inputDataDirBasePath(dataPeriod, tBinLabel)}/tree_data_{beamPolLabel}.root", )
    elif dataType == AnalysisConfig.DataType.ACCEPTED_PHASE_SPACE:
      return (f"{cfg.inputDataDirBasePath(dataPeriod, tBinLabel)}/tree_accepted_{beamPolLabel}.root", )
    elif dataType == AnalysisConfig.DataType.GENERATED_PHASE_SPACE:
      return (f"{cfg.inputDataDirBasePath(dataPeriod, tBinLabel)}/tree_thrown_{beamPolLabel}.root", )
    else:
      raise ValueError(f"Unknown data type: {dataType}")

  _inputFilePaths: Callable[[AnalysisConfig, AnalysisConfig.DataType, str, str, str], tuple[str, ...]] = field(
        default_factory = lambda: AnalysisConfig._default_inputFilePaths,
        repr            = False, 
        compare         = False,
    )

  def inputFilePaths(
    self,
    dataType:     AnalysisConfig.DataType,
    dataPeriod:   str,
    tBinLabel:    str,
    beamPolLabel: str,
  ) -> tuple[str, ...]:
    return self._inputFilePaths(self, dataType, dataPeriod, tBinLabel, beamPolLabel)

  def outputDataDirBasePath(
    self,
    dataPeriod: str,
    tBinLabel:  str,
  ) -> str:
    """Generates base path of directory with converted output data"""
    return f"{self.dataDirBaseName}/{dataPeriod}/{tBinLabel}/{self.subsystem.pairLabel}"

  def outputFilePath(
    self,
    dataType:     AnalysisConfig.DataType,
    dataPeriod:   str,
    tBinLabel:    str,
    beamPolLabel: str,
  ) -> str:
    if dataType == AnalysisConfig.DataType.REAL_DATA:
      return f"{self.outputDataDirBasePath(dataPeriod, tBinLabel)}/data_flat_{beamPolLabel}.root"
    elif dataType == AnalysisConfig.DataType.ACCEPTED_PHASE_SPACE:
      return f"{self.outputDataDirBasePath(dataPeriod, tBinLabel)}/phaseSpace_acc_flat_{beamPolLabel}.root"
    elif dataType == AnalysisConfig.DataType.GENERATED_PHASE_SPACE:
      return f"{self.outputDataDirBasePath(dataPeriod, tBinLabel)}/phaseSpace_gen_flat_{beamPolLabel}.root"
    else:
      raise ValueError(f"Unknown data type: {dataType}")

  @property
  def maxLPhys(self) -> int:
    """Returns maximum L of physical moments"""
    return self.maxL if isinstance(self.maxL, int) else self.maxL[0]

  @property
  def maxLMeas(self) -> int:
    """Returns maximum L of measured moments"""
    return self.maxL if isinstance(self.maxL, int) else self.maxL[1]

  @property
  def momentIndicesPhys(self) -> MomentIndices:
    """Returns moment indices for physical moments"""
    return MomentIndices(maxL = self.maxLPhys, polarized = (self.polarization is not None))

  @property
  def momentIndicesMeas(self) -> MomentIndices:
    """Returns moment indices for measured moments"""
    return MomentIndices(maxL = self.maxLMeas, polarized = (self.polarization is not None))

  @property
  def momentIndices(self) -> tuple[MomentIndices, MomentIndices]:
    """Returns moment indices for physical and measured moments"""
    return (self.momentIndicesPhys, self.momentIndicesMeas)

  @property
  def outFileDirName(self) -> str:
    """Returns name of directory into which all output will be written"""
    return f"{self.outFileDirBaseName}.maxL_{self.maxL if isinstance(self.maxL, int) else f'{self.maxL[0]}_{self.maxL[1]}'}"

  @property
  def outFileNamePrefix(self) -> str:
    """Returns name prefix prepended to output file names"""
    return "norm" if self.normalizeMoments else "unnorm"

  def init(
    self,
    createOutFileDir: bool = False,
  ) -> None:
    """Creates output directory and initializes member variables; needs to be called before passing the config to any consumer"""
    if isinstance(self.maxL, tuple):
      assert self.maxL[0] <= self.maxL[1], f"Maximum L for physical moments {self.maxL[0]=} must be smaller than or equal to maximum L for measured moments {self.maxL[1]=}"
    self.massBinning.var = self.binVarMass
    if createOutFileDir:
      Utilities.makeDirPath(self.outFileDirName)

  def loadData(
    self,
    dataType:             AnalysisConfig.DataType,
    additionalCuts:       Iterable[str] | None = None,  # optional additional cuts to be applied to loaded data
    additionalColumnDefs: dict[str, str]       = {},    # additional columns to define
  ) -> ROOT.RDataFrame | None:
    """Returns a ROOT RDataFrame with given data type, applies optional additional cuts, and defines optional additional columns"""
    df = None
    if dataType == AnalysisConfig.DataType.REAL_DATA:
      print(f"Loading real data from tree '{self.treeName}' in file '{self.dataFileName}'")
      df = ROOT.RDataFrame(self.treeName, self.dataFileName)
    elif dataType == AnalysisConfig.DataType.REAL_DATA_SIGNAL:
      print(f"Loading real-data signal events with weight = 1 from tree '{self.treeName}' in file '{self.dataFileName}'")
      df = ROOT.RDataFrame(self.treeName, self.dataFileName)
      if "Weight_SBsubtracted" in df.GetColumnNames():
        df = df.Filter("Weight_SBsubtracted == 1")
      # if there is no `Weight_SBsubtracted` column assume all events have weight 1
    elif dataType == AnalysisConfig.DataType.REAL_DATA_SIDEBAND:
      df = ROOT.RDataFrame(self.treeName, self.dataFileName)
      if "Weight_SBsubtracted" in df.GetColumnNames():
        df = df.Filter("Weight_SBsubtracted < 0")
        if df.Count().GetValue() > 0:
          print(f"Loading real-data sideband events with weight < 0 from tree '{self.treeName}' in file '{self.dataFileName}'")
        else:
          df = None  # no events with weight < 0
      else:
        # if there is no `Weight_SBsubtracted` column assume all events have weight 1 and hence there are no sideband events
        df = None
    elif dataType == AnalysisConfig.DataType.GENERATED_PHASE_SPACE:
      if self.psGenFileName is None:
        print("??? Warning: File name for generated phase-space data was not provided. Acceptance may not be calculated correctly.")
        df = None
      else:
        print(f"Loading generated phase-space data from tree '{self.treeName}' in file '{self.psGenFileName}'")
        df = ROOT.RDataFrame(self.treeName, self.psGenFileName)
    elif dataType == AnalysisConfig.DataType.ACCEPTED_PHASE_SPACE:
      if self.psAccFileName is None:
        print("??? Warning: File name for accepted phase-space data was not provided. Assuming perfect acceptance.")
        df = None
      else:
        print(f"Loading accepted phase-space data from tree '{self.treeName}' in file '{self.psAccFileName}'")
        df = ROOT.RDataFrame(self.treeName, self.psAccFileName)
    else:
      raise ValueError(f"Unknown data type: {dataType}")
    # add additional columns if requested
    if df is not None and additionalColumnDefs:
      for columnName, columnFormula in additionalColumnDefs.items():
        print(f"Defining additional column '{columnName}' = '{columnFormula}'")
        df = defineOverwriteRDataFrame(df, columnName, columnFormula)
    # apply additional cuts
    if df is not None and additionalCuts is not None:
      for cut in additionalCuts:
        print(f"Applying additional cut: '{cut}'")
        df = df.Filter(cut)
    return df


# configuration for Nizar's polarized gamma p -> (eta pi0) p data, with eta -> gamma gamma
CFG_IO_TEST_LMAX4_ZEROPWAVES = AnalysisConfig(
  frame              = AnalysisConfig.CoordSysType.GJ,
  subsystem          = SubsystemInfo(  # eta pi0 subsystem to analyze; eta is the analyzer
    lvALabel          = "eta",     # label of eta Lorentz-vector
    lvBLabel          = "pi0",     # label of pi0 Lorentz-vector
    lvRecoilLabel     = "recoil",  # label of recoil-proton Lorentz-vector
    pairLabel         = "EtaPi0",
    ATLatexLabel      = "#it{#eta}",
    BTLatexLabel      = "#it{#pi}^{0}",
    recoilTLatexLabel = "#it{p}",
    pairTLatexLabel   = "#it{#eta}#it{#pi}^{0}",
  ),
  dataDirBaseName    = "./data/IO_test_Lmax4_ZeroPWaves/",
  dataPeriods        = ("merged", ),  # merged Phase-I + Spring 2020 data with different polarization values
  tBinLabels         = (
    "t010020",
  ),
  beamPolLabels      = ("All", ),  # input files contain all beam polarization orientations merged together
  inputDataFormats   = {
    AnalysisConfig.DataType.REAL_DATA             : AnalysisConfig.DataFormat.AMPTOOLS,
    AnalysisConfig.DataType.ACCEPTED_PHASE_SPACE  : AnalysisConfig.DataFormat.AMPTOOLS,
    AnalysisConfig.DataType.GENERATED_PHASE_SPACE : AnalysisConfig.DataFormat.AMPTOOLS,
  },
  # treeName           = "EtaPi0",
  treeName           = "kin",  # tree name in input files; same for all data types 
  dataFileName       = "data_flat_All.root",
  psAccFileName      = "phaseSpace_acc_flat_All.root",
  psGenFileName      = "phaseSpace_gen_flat_All.root",
  polarization       = "Pol",  # read polarization from tree column
  # maxL               = 4,
  outFileDirBaseName = "./Moments.out/IO_test_m128132_Lmax4_ZeroPWaves/",
  binVarMass         = KinematicBinningVariable(
    name      = "Mpi0eta",
    label     = "#it{m}_{#it{#eta}#it{#pi}^{0}}",
    unit      = "GeV/#it{c}^{2}",
    nmbDigits = 3,
  ),
  massBinning        = HistAxisBinning(nmbBins = 1, minVal = 1.28, maxVal = 1.32),  # 40 MeV wide bins
)