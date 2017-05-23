
FSdir=/projects2/udall/freesurferS1S2
taskDir=/projects2/udall/task/fnirt
RC4files=$(notdir $(wildcard $(taskDir)/RC4[1-2][0-9][0-9]))
FSfiles=$(RC4files:RC4%=RC4%_brain.nii.gz)

all: $(FSfiles)

%_brain.nii.gz:
	mri_convert $(FSdir)/$*.s1/mri/brain.mgz $@
