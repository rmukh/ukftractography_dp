start=`date +%s`

SRC="../UKFTractography"

#logFile Name
logFile="log.txt"

# BINARY
BINARY='../build/UKFTractography-build/UKFTractography/bin/UKFTractography'

# VOLUME
dwi_path="/rfanfs/pnl-zorro/home/rinat/HCP/100307-dwi-new-convert.nhdr"

# MASK
mask_path="/rfanfs/pnl-zorro/home/rinat/HCP/100307_diff_mask.nrrd"

# OUTPUT FIBER
output_path="/rfanfs/pnl-zorro/home/rinat/results/fast_rtop1_500_sl_0.5_qm_0.0001_FW_0.65_odf_0.3_seeds_3.vtk"

#csf_path="/data/pnl/home/rz892/HCP/spm/csf_mask.nrrd"
wm_path="/rfanfs/pnl-zorro/home/rinat/HCP/spm/wm_mask.nrrd"
seeds_path="/home/rinat/Desktop/HCP_subj/seeds_stem.nrrd"
echo "Start time $(date)" | tee $logFile

# --seedsFile $seeds_path \
eval $BINARY \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --tracts $output_path \
 --wmFile $wm_path \
 --seedsPerVoxel 3 \
 --diffusionPropagator \
 --minRTOP1stop 500 \
 --Qm 0.0001 \
 --stepLength 0.5 \
 --maxODFthresh 0.3 \
 --recordNMSE \
 --recordWeights \
 --recordRTOP \
 --recordState \
 --recordUncertainties \
 --numTensor 3 | tee -a $logFile
 end=`date +%s`

runtime=$((end - start))
echo "Output file name $output_path" | tee -a $logFile
echo "CPU Runtime was $runtime sec."  | tee -a $logFile

#diff $output_path $SRC/ukf/Data/Baseline/seeds_tc.vtk

