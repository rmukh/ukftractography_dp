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
output_path='./seeds_tc.vtk'

eval $BINARY \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --diffusionPropagator \
 --tracts $output_path \
 --seedsFile $seeds_path \
 --minRTOP 20 \
 --seedsPerVoxel 5 \
 --numTensor 3 | tee -a $logFile
 end=`date +%s`

runtime=$(python -c "print(${end} - ${start})")
echo "Output file name $output_path" | tee -a $logFile
echo "CPU Runtime was $runtime"  | tee -a $logFile

#diff $output_path $SRC/ukf/Data/Baseline/seeds_tc.vtk