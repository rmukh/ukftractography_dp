start=`date +%s`

SRC="../UKFTractography"

#logFile Name
logFile="output/log1.txt"

# BINARY
BINARY='../build/UKFTractography-build/UKFTractography/bin/UKFTractography'

# VOLUME
dwi_path="/data/pnl/home/rz892/HCP/100307-dwi.nhdr"

# MASK
mask_path="/data/pnl/home/rz892/HCP/nodif_brain_mask.nrrd"

# OUTPUT FIBER
output_path="/data/pnl/home/rz892/HCP/results/NO_wm_mask_NO_csf_mask_rtop1_500_sl_0.3_qm_0.005_FW_0.65_GA_0.1_odf_0.3.vtk"

#csf_path="/data/pnl/home/rz892/HCP/spm/csf_mask.nrrd"
#wm_path="/data/pnl/home/rz892/HCP/spm/wm_mask.nrrd"
echo "Start time $(date)" | tee $logFile

# --seedsFile $seeds_path \
eval $BINARY \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --tracts $output_path \
 --seedsPerVoxel 1 \
 --diffusionPropagator \
 --Qm 0.005 \
 --stepLength 0.3 \
 --maxODFthresh 0.3 \
 --recordNMSE \
 --recordWeights \
 --recordRTOP \
 --numTensor 3 | tee -a $logFile
 end=`date +%s`

runtime=$((end - start))
echo "Output file name $output_path" | tee -a $logFile
echo "CPU Runtime was $runtime sec."  | tee -a $logFile

#diff $output_path $SRC/ukf/Data/Baseline/seeds_tc.vtk

