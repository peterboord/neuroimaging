subjectsDir=/NAS_II/Projects/Udall/subjects
#/RC4101/session1/fMRIb0_mag_brain.nii.gz
sessions=$(notdir $(wildcard /projects2/udall/folders/sessions/RC4[1-2][0-9][0-9]-[1-2]))
EF_UD_warp=$(sessions:%=%/b0unwarp.feat/unwarp/EF_UD_warp.nii.gz)

all: $(EF_UD_warp)

RC4%/b0unwarp.feat/unwarp/EF_UD_warp.nii.gz:
	session=RC4$* ;\
	subject=$$(echo $$session|cut -d'-' -f 1) ;\
	sessionNr=$$(echo $$session|cut -d'-' -f 2) ;\
	\cp --preserve=timestamps $(subjectsDir)/$$subject/session$$sessionNr/fMRIb0_mag_brain.nii.gz $$session/. ;\
	\cp --preserve=timestamps $(subjectsDir)/$$subject/session$$sessionNr/fMRIb0_mag.nii.gz $$session/. ;\
	\cp --preserve=timestamps $(subjectsDir)/$$subject/session$$sessionNr/fMRIb0_phase.nii.gz $$session/. ;\
	cat b0unwarp.fsf |sed 's/SESSION/RC4$*/g' > fsf/RC4$*_b0unwarp.fsf
