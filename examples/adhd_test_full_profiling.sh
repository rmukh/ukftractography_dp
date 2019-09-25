start=`date +%s`

SRC="../UKFTractography"

#logFile Name
logFile="log.txt"

# BINARY
BINARY='../build/UKFTractography-build/UKFTractography/bin/UKFTractography'

# VOLUME
dwi_path="/home/rinat/Desktop/ukftests/dwi-Aligned-Ed-Bet-Merged-unring-Epi.nhdr"

# MASK
mask_path="/home/rinat/Desktop/ukftests/epi_corrected_tensormask.nrrd"

seeds_path="/home/rinat/Desktop/ukftests/Segmentation-label_363x.nrrd"

# OUTPUT FIBER
output_path='/home/rinat/Desktop/ukftests/seeds_tc_profiling.vtk'

# --seedsFile $seeds_path \
valgrind --tool=callgrind \
 ../build/UKFTractography-build/UKFTractography/bin/UKFTractography \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --tracts $output_path \
 --seedsFile $seeds_path \
 --seedsPerVoxel 1 \
 --diffusionPropagator \
 --minRTOP1stop 600 \
 --recordNMSE \
 --recordWeights \
 --recordRTOP \
 --seedingThreshold 0.18 \
 --numTensor 3 | tee $logFile

