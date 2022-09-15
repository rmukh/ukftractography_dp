start=`date +%s`

SRC="../UKFTractography"

#logFile Name
logFile="log.txt"

# BINARY
BINARY='../build/UKFTractography-build/UKFTractography/bin/UKFTractography'

# VOLUME
dwi_path="$SRC/Data/Input/dwi.nhdr"

# MASK
mask_path="$SRC/Data/Input/dwi-mask.nhdr"

# SEEDS
seeds_path="$SRC/Data/Input/seeds_tc.nhdr"

# OUTPUT FIBER
output_path='./fibers_tc.vtk'

# --seedsFile $seeds_path \
eval $BINARY \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --seedsFile $seeds_path \
 --tracts $output_path \
 --stepLength 0.2 \
 --seedsPerVoxel 5 \
 --recordRTOP \
 --recordUncertainties | tee $logFile
 end=`date +%s`

runtime=$((end - start))
echo "Output file name $output_path" | tee -a $logFile
echo "CPU Runtime was $runtime sec."  | tee -a $logFile
