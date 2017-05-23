#!/bin/bash
MNI2mm=/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz
subjDir=/projects2/act-plus/subjects
subjPaths=$(wildcard $(subjDir)/109[0-9][0-9][0-9])
physioDir=/projects2/act-plus/physio
physioFiles=rest_cardio rest_resp
.PHONY: all

#subjPaths=$(subjDir)/109001
#subjPaths=$(subjDir)/109055
#physioFiles=rest_resp

all: $(subjPaths)
	
$(subjDir)/%: FORCE
	mkdir -p $* ;\
	cd $* ;\
	ln -sf $(subjDir)/$*/session1/rest/rest_e001.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/rest/rest_e002.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/rest/rest_e003.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/rest/rest_e00213_medn.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/rest/rest_e00213_tsoc.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/rest/meica.rest_e00213/eBvrmask.nii.gz mask.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/memprage/T1.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/memprage/T1_brain.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/fast_segmentation/T1_brain_pve_0.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/fast_segmentation/T1_brain_pve_1.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/fast_segmentation/T1_brain_pve_2.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/fast_segmentation/T1_brain_mixeltype.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/fast_segmentation/T1_brain_pveseg.nii.gz ;\
	ln -sf $(subjDir)/$*/session1/fieldmap/B0_mag_fMRI.nii.gz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/wm.mgz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/aseg.mgz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/aparc+aseg.mgz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/orig.mgz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/rawavg.mgz ;\
	for file in $(physioFiles); do \
		ln -sf $(physioDir)/$*_$${file}_ds.txt $${file}_ds.txt ;\
		ln -sf $(physioDir)/$*_$${file}_peaks.txt $${file}_peaks.txt ;\
	done ;\
	fslmaths fcDiff_rest_e00213_medn -abs -add mean_func add_abs_fcDiff_rest_e00213_medn ;\
	fslmaths pic_rest_e00213_medn.nii.gz -mul -0.5 -add 0.5 dis_pic_rest_e00213_medn.nii.gz ;\
	cat ../design.fsf | sed "s/SUBJECT/$*/g" > ../fsf/$*.fsf
temp:
	#mri_convert orig.mgz orig.nii.gz ;\
	#fslmaths rest_e00213_tsoc.nii.gz -Tmean mean_func ;\
	#run_pic $$MCRROOT . rest_e00213_medn.nii.gz mask.nii.gz
	flirt -in mean_func -ref T1_brain -omat rest_to_T1.mat -out rest_to_T1 ;\
	flirt -in T1_brain -ref $(MNI2mm) -omat T1_to_std.mat -out T1_to_std ;\
	convert_xfm -omat rest_to_std.mat -concat T1_to_std.mat rest_to_T1.mat ;\
	convert_xfm -omat std_to_rest.mat -inverse rest_to_std.mat ;\
	flirt -in corrAbsFcDiffPic -ref $(MNI2mm) -applyxfm -init rest_to_std.mat -out mni_corrAbsFcDiffPic

FORCE:
	
