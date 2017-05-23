subjectDir=/projects2/act-plus/subjects
subjects=$(wildcard $(subjectDir)/109[0-9][0-9][0-9])
physioDir=/projects2/act-plus/physio
physioFiles=rest_cardio rest_resp
.PHONY: all

#subjects=109003
#physioFiles=rest_resp

all: $(subjects)

$(subjectDir)/%: FORCE
	mkdir -p $* ;\
	cd $* ;\
	ln -sf $(subjectDir)/$*/session1/rest/rest_e001.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/rest/rest_e002.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/rest/rest_e003.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/memprage/T1.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/memprage/T1_brain.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/fast_segmentation/T1_brain_pve_0.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/fast_segmentation/T1_brain_pve_1.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/fast_segmentation/T1_brain_pve_2.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/fast_segmentation/T1_brain_mixeltype.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/fast_segmentation/T1_brain_pveseg.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/fieldmap/B0_mag_fMRI.nii.gz ;\
	ln -sf $(subjectDir)/$*/session1/fieldmap/B0_phase_fMRI.nii.gz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/wm.mgz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/aseg.mgz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/aparc+aseg.mgz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/orig.mgz ;\
	ln -sf /projects2/act-plus/freesurfer/$*/mri/rawavg.mgz ;\
	for file in $(physioFiles); do \
		ln -sf $(physioDir)/$*_$${file}_ds.txt $${file}_ds.txt ;\
		ln -sf $(physioDir)/$*_$${file}_peaks.txt $${file}_peaks.txt ;\
	done

FORCE:
	
