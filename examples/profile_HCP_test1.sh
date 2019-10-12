start=`date +%s`

SRC="../UKFTractography"

#logFile Name
logFile="log.txt"

# BINARY
BINARY='../build/UKFTractography-build/UKFTractography/bin/UKFTractography'

# VOLUME
dwi_path="/home/rinat/Desktop/HCP_subj/100307-dwi-new-convert.nhdr"

# MASK
mask_path="/home/rinat/Desktop/HCP_subj/100307_diff_mask.nrrd"

# OUTPUT FIBER
output_path="/home/rinat/Desktop/ukftests/HCP/HCP_demo_output.vtk"

#csf_path="/data/pnl/home/rz892/HCP/spm/csf_mask.nrrd"
#wm_path="/data/pnl/home/rz892/HCP/spm/wm_mask.nrrd"
seeds_path="/home/rinat/Desktop/HCP_subj/seeds_stem.nrrd"
echo "Start time $(date)" | tee $logFile

# --seedsFile $seeds_path \
sudo perf record -g $BINARY \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --tracts $output_path \
 --seedsFile $seeds_path \
 --seedsPerVoxel 3 \
 --minRTOP1stop 500 \
 --Qm 0.0001 \
 --stepLength 0.5 \
 --maxODFthresh 0.3 \
 --recordNMSE \
 --recordWeights \
 --recordRTOP

echo "Output file name $output_path" | tee -a $logFile

#diff $output_path $SRC/ukf/Data/Baseline/seeds_tc.vtk

