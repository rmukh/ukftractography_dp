#include "tractography.h"
#include "utilities.h"

#include <itkMacro.h> // needed for ITK_VERSION_MAJOR
#if ITK_VERSION_MAJOR < 5
#include "itkMultiThreader.h"
#else
#include "itkMultiThreaderBase.h"
#endif

#include "UKFTractographyCLP.h"

// TODO make configurable?
static const bool verbose = true;

namespace
{
void ukf_setAndTell(ukfPrecisionType &x, const ukfPrecisionType y, const std::string &name)
{
  if (verbose)
  {
    x = y;
    std::cout << "- " << name << ": " << y << std::endl;
  }
}

void ukf_tell(const ukfPrecisionType &x, const std::string &name)
{
  if (verbose)
  {
    std::cout << "* " << name << ": " << x << std::endl;
  }
}
}; // anonymous namespace

int ukf_parse_cli(int argc, char **argv, UKFSettings &s)
{
  PARSE_ARGS;

  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " --dwiFile dMRI file --maskFile mask file --tracts output fiber tracts file name" << std::endl;
    std::cout << "To see other parameters call " << argv[0] << " --help" << std::endl;
    return 1;
  }

  ukfPrecisionType l_stoppingThreshold = stoppingThreshold;
  ukfPrecisionType l_stepLength = stepLength;
  ukfPrecisionType l_recordLength = recordLength;
  ukfPrecisionType l_maxHalfFiberLength = maxHalfFiberLength;
  ukfPrecisionType l_seedingThreshold = seedingThreshold;
  ukfPrecisionType l_Qm = Qm;
  ukfPrecisionType l_Ql = Ql;
  ukfPrecisionType l_Qw = Qw;
  ukfPrecisionType l_Qt = Qt;
  ukfPrecisionType l_Qwiso = Qwiso;
  ukfPrecisionType l_Rs = Rs;

  ukfPrecisionType l_minRTOP1stop = minRTOP1stop;
  ukfPrecisionType l_maxNMSE = maxNMSE;
  ukfPrecisionType l_maxUKFIterations = maxUKFIterations;
  ukfPrecisionType l_maxODFthresh = maxODFthresh;
  ukfPrecisionType l_FWthresh = FWthresh;

  // If sigmaSignal is not set minimum of voxel size is used for interpolation
  ukfPrecisionType SIGMA_SIGNAL = sigmaSignal;

  // HANDLE ERRORNOUS INPUT
  if (dwiFile.empty() || maskFile.empty() || tracts.empty())
  {
    std::cout << "Error! Must indicate DWI data, mask and tracts output files!" << std::endl
              << std::endl;
    return 1; //This is to indicate that the module returns with error
  }

  if (l_maxHalfFiberLength <= 0)
  {
    std::cout << "Invalid maximum half fiber length!" << std::endl;
    return 1;
  }

  if (std::ceil(l_maxHalfFiberLength / l_stepLength) <= 1)
  {
    std::cout << "Too large step length or too small fiber cutoff limit!" << std::endl;
    return 1;
  }

  if (l_recordLength < l_stepLength)
  {
    std::cout << "recordLength should be greater than stepLength" << std::endl;
    return 1;
  }

  // SETTING THE DEFAULT PARAMETERS
  std::cout << "Using the Diffustion Propagator Bi-exponential Spherical Ridgelets model. Setting the default parameters accordingly:\n";
  std::cout << "\"*\": set by user\n";
  std::cout << "\"-\": default setting\n";

  if (labels.size() == 0)
  {
    labels.push_back(1); //Default to use label 1
  }

  if (l_seedingThreshold == 0.18)
  {
    ukf_setAndTell(l_seedingThreshold, FULL_BRAIN_MEAN_SIGNAL_MIN, "seedingThreshold"); 
  }
  else
  {
    ukf_tell(l_seedingThreshold, "seedingThreshold");
  }

  if (l_Qm == 0.0)
  {
    ukf_setAndTell(l_Qm, 0.001, "Qm");
  }
  else
  {
    ukf_tell(l_Qm, "Qm");
  }

  if (l_Ql == 0.0)
  {
    ukf_setAndTell(l_Ql, 150.0, "Ql");
  }

    if (l_Qt == 0.0)
    {
      ukf_setAndTell(l_Qt, 50.0, "Qt");
    }
    else
    {
      ukf_tell(l_Qt, "Qt");
    }

    if (l_minRTOP1stop == 500.0)
    {
      ukf_setAndTell(l_minRTOP1stop, 500.0, "minRTOP1stop");
    }
    else
    {
      ukf_tell(l_minRTOP1stop, "minRTOP1stop");
    }

    if (l_maxNMSE == 0.15)
    {
      ukf_setAndTell(l_maxNMSE, 0.15, "maxNMSE");
    }
    else
    {
      ukf_tell(l_maxNMSE, "maxNMSE");
    }

    if (l_maxUKFIterations < 0.0)
    {
      std::cout << "Error: maxUKFIterations cannot be negative. Exiting" << std::endl;
      exit(1);
    }
    if (l_maxUKFIterations == 5)
    {
      ukf_setAndTell(l_maxUKFIterations, 5, "maxUKFIterations");
    }
    else
    {
      ukf_tell(l_maxUKFIterations, "maxUKFIterations");
    }

    if (l_maxODFthresh == 0.7)
    {
      ukf_setAndTell(l_maxODFthresh, 0.7, "maxODFthresh");
    }
    else
    {
      ukf_tell(l_maxODFthresh, "maxODFthresh");
    }

    if (l_FWthresh == 0.65)
    {
      ukf_setAndTell(l_FWthresh, 0.65, "FWthresh");
    }
    else
    {
      ukf_tell(l_FWthresh, "FWthresh");
    }
  

  if (l_Rs == 0.0)
  {
      ukf_setAndTell(l_Rs, 0.015, "Rs");
  }
  else
  {
    ukf_tell(l_Rs, "Rs");
  }

  if (l_stepLength == 0.3)
  {
    ukf_setAndTell(l_stepLength, 0.3, "stepLength");
  }
  else
  {
    ukf_tell(l_stepLength, "stepLength");
  }

  if (l_recordLength == 0.9)
  {
    ukf_setAndTell(l_recordLength, 0.9, "recordLength");
  }
  else
  {
    ukf_tell(l_recordLength, "recordLength");
  }

    if (l_Qw == 0.0)
    {
        ukf_setAndTell(l_Qw, 0.002, "Qw");
    }
    else
    {
      ukf_tell(l_Qw, "Qw");
    }
  
    ukf_setAndTell(l_Qwiso, 0.002, "Qwiso");
  

  if (l_stoppingThreshold == 0.1)
  {
    ukf_setAndTell(l_stoppingThreshold, 0.1, "stoppingThreshold");
  }
  else
  {
    ukf_tell(l_stoppingThreshold, "stoppingThreshold");
  }

  if (seedsPerVoxel == 1)
  {
    std::cout << "- seedsPerVoxel: " << seedsPerVoxel << std::endl;
  }
  else
  {
    std::cout << "* seedsPerVoxel: " << seedsPerVoxel << std::endl;
  }

  if (seedsPerVoxel != static_cast<int>(seedsPerVoxel) && seedsPerVoxel > 1.0)
  {
    cout << "You set seedsPerVoxel = " << seedsPerVoxel << " which is real number that > 1, "
                                                           "so it will be rounded to nearest integer = "
         << static_cast<int>(seedsPerVoxel) << endl;
    cout << "If you want to specify fraction of seed voxels than restart program with value in the range (0,1)" << endl;
    seedsPerVoxel = static_cast<int>(seedsPerVoxel);
  }

  // initializing settings
  //UKFSettings& s -- from function argument
  {
    s.record_nmse = recordNMSE;
    s.record_trace = recordTrace;
    s.record_state = recordState;
    s.record_cov = recordCovariance;
    s.record_free_water = recordFreeWater;
    s.record_tensors = recordTensors;
    s.record_weights = recordWeights;
    s.record_uncertainties = recordUncertainties;
    s.record_rtop = recordRTOP;
    s.transform_position = true; // TODO hard-coded :/
    s.store_glyphs = storeGlyphs;
    s.mean_signal_min = l_stoppingThreshold;
    s.seeding_threshold = l_seedingThreshold;
    s.seeds_per_voxel = seedsPerVoxel;
    s.stepLength = l_stepLength;
    s.recordLength = l_recordLength;
    s.maxHalfFiberLength = maxHalfFiberLength;
    s.labels = labels;
    s.num_threads = numThreads;

    s.Qm = l_Qm;
    s.Ql = l_Ql;
    s.Qw = l_Qw;
    s.Qt = l_Qt;
    s.Qwiso = l_Qwiso;
    s.Rs = l_Rs;

    s.rtop1_min_stop = l_minRTOP1stop;
    s.max_nmse = l_maxNMSE;
    s.maxUKFIterations = l_maxUKFIterations;
    s.max_odf_threshold = l_maxODFthresh;
    s.fw_thresh = l_FWthresh;

    // TODO these should be header-initialized once we use C++11
    s.p0 = P0;
    s.sigma_signal = SIGMA_SIGNAL;
    s.sigma_mask = SIGMA_MASK;
    s.min_radius = MIN_RADIUS;

    s.output_file = tracts;
    s.dwiFile = dwiFile;
    s.seedsFile = seedsFile;
    s.maskFile = maskFile;
    s.csfFile = csfFile;
    s.writeAsciiTracts = writeAsciiTracts;
    s.writeUncompressedTracts = writeUncompressedTracts;
  }

  return EXIT_SUCCESS;
}
