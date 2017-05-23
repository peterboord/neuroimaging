#!/bin/bash
# usage:
# epi_reg_withB0.sh funcParFile b0parFile examCardSenseFactor meanFuncFile outFile_in_T1 T1 T1_brain
# e.g.
# funcParFile="*_WIP_ME-RS_SENSE_3_1.PAR"
# b0parFile="*_WIP_B0_*.PAR"
# examCardSenseFactor=2.5
# meanFuncFile=raw_mean_func
# outFile_in_T1=raw_func_to_T1
# T1=T1
# T1_brain=T1_brain
funcParFile=$1
b0parFile=$2
examCardSenseFactor=$3
meanFuncFile=$4
outFile_in_T1=$5
T1=$6
T1_brain=$7

b0par2nii.sh $b0parFile
mv -f *B0*x1.nii.gz B0_Mag.nii.gz
mv -f *B0*x2.nii.gz B0_Phase.nii.gz
# B0_Mag_brain_mask
echo bet
bet B0_Mag.nii.gz B0_Mag_brain.nii.gz -R -m
# B0_Mag_brain_restore
echo fast
fast -B -t 2 -o B0_Mag_brain B0_Mag_brain
rm B0_Mag_brain_pve*.nii.gz B0_Mag_brain_seg*.nii.gz B0_Mag_brain_mixeltype.nii.gz
# B0_phase_rescaled
echo fslmaths
fslmaths -dt float B0_Phase -div 500 -mul 3.14 B0_phase_rescaled -odt float
# B0_Mag_brain_mask_ero1
fslmaths B0_Mag_brain_mask -ero B0_Mag_brain_mask_ero1
# B0_phase_rescaled_unwrapped
echo prelude
prelude -p B0_phase_rescaled -a B0_Mag -o B0_phase_rescaled_unwrapped -m B0_Mag_brain_mask_ero1
# B0_Phase_rads
echo fslmaths
fslmaths B0_phase_rescaled_unwrapped -mul 1000 B0_Phase_rads
rm B0_phase_rescaled.nii.gz B0_phase_rescaled_unwrapped.nii.gz
echo epi_reg
epi_reg --epi=$meanFuncFile --t1=$T1 --t1brain=$T1_brain --fmap=B0_Phase_rads --fmapmag=B0_Mag --fmapmagbrain=B0_Mag_brain_restore \
        --echospacing=`GetEchoSpacing_pb.sh $funcParFile $examCardSenseFactor` --pedir=y --noclean --out=$outFile_in_T1

