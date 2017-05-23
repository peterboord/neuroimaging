restDir=/projects2/udall/standardFC
sessions=$(notdir $(wildcard $(restDir)/RC4[1-2][0-9][0-9]-[1-2]))

sessions=RC4101-1

prefiltered_func_data_mcf=$(sessions:%=%/prefiltered_func_data_mcf.mat)
rest_brain_res_reg=$(sessions:%=%/rest_brain_res_reg.nii.gz)
rawDir=/projects2/udall/rest/raw
betDir=/projects2/udall/rest/bet

all: $(rest_brain_res_reg)

RC4%/prefiltered_func_data_mcf.mat:
	echo $@ ;\
	session=RC4$* ;\
	mkdir -p $$session ;\
	cd $$session ;\
	mcflirt -stats -mats -plots prefiltered_func_data_mcf.par -in $(rawDir)/$${session}_rest.nii.gz -out prefiltered_func_data_mcf
RC4%/prefiltered_func_data_mcf.mat/inv_MAT_0000: RC4%/prefiltered_func_data_mcf.mat
	cd $$(dirname $@) ;\
	for mat in $$(ls MAT_????); do \
		convert_xfm -omat inv_$$mat -inverse $$mat ;\
	done
# problem: mean across vols influenced by global signal of all time points, so intensity of inv_ref_mean could be larger at one timepoint due to large global means elsewhere
# solution: normalize each volume with global signal, calc ref_mean, inv xfm ref_mean, then scale back with original global mean
# for each vol
# mask, divide by spatial mean, xfm to ref
# av across all vols for ref_mean
# inv xfm from ref
# scale back 
#
# mask func with bet mask
# global signal = mean across non-zero voxels
# loop1: divide (norm) by global signal - reduces the bias from images with high mean signal
# xfm normed images to ref using mcflirt xfms - use nearestneighbour to avoid blurring gm/wm edge
# merge images and calc *median* of registered norm images - median *classifies* voxel as gm/wm
# loop2: inv xfm mean reg norm images back to original image space - blurs gm/wm edge
# scale back to original range with global signal
# residual image = original image - xfm'ed ref
# xfm res to ref using mcflirt xfms
# merge images and calc mean of registered norm images
#
# Order of volreg & tshift?
# Subtraction of gm/wm estimate is best before tshift as gm/wm independent of slice timing. 
# Also, gm/wm misreg likely to cause larger transients than slice timing lags, so removing gm/wm contrast will help tshift filtering.
#
# Order of volreg & retroicor & tshift
# Frequencies above Nyquist can't be properly interpolated (ref 3dTshift), so minimize physio (retroicor) before tshift.

RC4%/rest_brain_res_reg.nii.gz: RC4%/prefiltered_func_data_mcf.mat/inv_MAT_0000
	echo $@ ;\
	session=RC4$* ;\
	cd $$session/prefiltered_func_data_mcf.mat ;\
	fslmaths $(rawDir)/$${session}_rest -mas $(betDir)/$$session/rest_brain_mask rest_brain_masked ;\
	fslstats -t rest_brain_masked -M > global_signal ;\
	fslroi rest_brain_masked rest_brain_masked_z0 0 -1 0 -1 0 1 0 -1 ;\
	fslmerge -z rest_brain_masked_pad0 rest_brain_masked_z0 rest_brain_masked ;\
	fslroi rest_brain_masked rest_brain_masked_zend 0 -1 0 -1 $$(($$(fslval rest_brain_masked dim3)-1)) 1 0 -1 ;\
	fslmerge -z rest_brain_masked_padded rest_brain_masked rest_brain_masked_zend ;\
	fslsplit rest_brain_masked_padded rest_brain_ ;\
	volNrUnpadded=0 ;\
	for global_signal in $$(cat global_signal); do \
		volNr=$$(printf "%04d" $$volNrUnpadded) ;\
		mat=MAT_$$volNr ;\
		fslmaths rest_brain_$${volNr}.nii.gz -div $$global_signal rest_brain_$${volNr}_norm ;\
		flirt -in rest_brain_$${volNr}_norm -out rest_brain_$${volNr}_norm_reg -ref rest_brain_$${volNr}_norm -applyxfm -init $$mat ;\
		volNrUnpadded=$$(($$volNrUnpadded + 1)) ;\
	done ;\
	fslmerge -t rest_brain_norm_reg $$(ls rest_brain_*_norm_reg.nii.gz) ;\
	fslmaths rest_brain_norm_reg -Tmedian rest_brain_norm_reg_median ;\
	volNrUnpadded=0 ;\
	for global_signal in $$(cat global_signal); do \
		volNr=$$(printf "%04d" $$volNrUnpadded) ;\
		mat=MAT_$$volNr ;\
		flirt -in rest_brain_norm_reg_median -out inv_rest_brain_norm_reg_median_$$mat -ref rest_brain_norm_reg_median -applyxfm -init inv_$$mat ;\
		fslmaths inv_rest_brain_norm_reg_median_$$mat -mul $$global_signal inv_rest_brain_norm_reg_median_$$mat ;\
		fslmaths rest_brain_$${volNr}.nii.gz -sub inv_rest_brain_norm_reg_median_$$mat rest_brain_$${volNr}_res ;\
		flirt -in rest_brain_$${volNr}_res -out rest_brain_$${volNr}_res_reg -ref rest_brain_$${volNr}_res -applyxfm -init $$mat ;\
		volNrUnpadded=$$(($$volNrUnpadded + 1)) ;\
	done ;\
	fslmerge -t ../rest_brain_res_reg $$(ls rest_brain_*_res_reg.nii.gz) ;\
	fslroi ../rest_brain_res_reg ../rest_brain_res_reg 0 -1 0 -1 1 $$(fslval rest_brain_masked dim3) 0 -1 ;\
	fslmaths ../rest_brain_res_reg -Tmean ../rest_brain_res_reg_mean
		
