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

# SEEDS
seeds_path="/home/rinat/Desktop/ukftests/Segmentation-label_363x.nrrd"

# OUTPUT FIBER
output_path='/home/rinat/Desktop/ukftests/adhd363/seeds_tc_363x_multi_speed.vtk'

# --seedsFile $seeds_path \
valgrind -v $BINARY \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --tracts $output_path \
 --stepLength 0.5 \
 --seedsFile $seeds_path

