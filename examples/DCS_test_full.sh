start=`date +%s`

SRC="../UKFTractography"

#logFile Name
logFile="log.txt"

# BINARY
BINARY='../build/UKFTractography-build/UKFTractography/bin/UKFTractography'

# VOLUME
dwi_path="/home/rinat/Desktop/ukftests/dwi_0008_ed.nhdr"

# MASK
mask_path="/home/rinat/Desktop/ukftests/dwi_0008_ed_mask_edited.nhdr"

# OUTPUT FIBER
output_path='/home/rinat/Desktop/ukftests/dcs_full_post2.vtk'

csf_path='/home/rinat/Desktop/ukftests/dwi_0008_CSF.nrrd'
wm_path='/home/rinat/Desktop/ukftests/dwi_0008_WM.nrrd'

#seeds_path='/home/rinat/Desktop/ukftests/dwi_0008_seeds.nrrd'

# --seedsFile $seeds_path \
eval $BINARY \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --csfFile $csf_path \
 --wmFile $wm_path \
 --tracts $output_path \
 --seedsPerVoxel 0.01 \
 --diffusionPropagator \
 --minRTOP 20 \
 --recordNMSE \
 --recordWeights \
 --recordRTOP \
 --seedingThreshold 0.3 \
 --numTensor 3 | tee $logFile
 end=`date +%s`

runtime=$((end - start))
echo "Output file name $output_path" | tee -a $logFile
echo "CPU Runtime was $runtime sec."  | tee -a $logFile

#diff $output_path $SRC/ukf/Data/Baseline/seeds_tc.vtk

