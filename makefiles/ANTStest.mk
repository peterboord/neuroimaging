
subjDir=../task/RC4101/session1
all: \
$(subjDir)/T1_brain_in_template.nii.gz \
$(subjDir)/T1_brain_in_template.nii.gz \
$(subjDir)/T1_brain_in_MNI_WtAtWaAa.nii.gz \
$(subjDir)/T1_brain_in_MNI_WaAaWtAt.nii.gz \
$(subjDir)/T1_brain_in_MNI_AtWtAaWa.nii.gz \
$(subjDir)/T1_brain_in_MNI_AaWaAtWt.nii.gz \
$(subjDir)/T1_brain_in_MNI_AaAtWaWt.nii.gz \
$(subjDir)/T1_brain_in_MNI_WtWaAtAa.nii.gz \
$(subjDir)/T1_brain_in_MNI_WaWtAaAt.nii.gz \
$(subjDir)/T1_brain_in_MNI_AtAaWtWa.nii.gz

# two-step, step 1
$(subjDir)/T1_brain_in_template.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_template.nii.gz \
	-R udall_template.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt

# two-step, step 2
$(subjDir)/T1_brain_in_template.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain_in_template.nii.gz \
	$(subjDir)/T1_brain_in_MNI.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	ANTStemplate_to_MNI_Warp.nii.gz \
	ANTStemplate_to_MNI_Affine.txt

# WtAtWaAa
$(subjDir)/T1_brain_in_MNI_WtAtWaAa.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_MNI_WtAtWaAa.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt \
	ANTStemplate_to_MNI_Warp.nii.gz \
	ANTStemplate_to_MNI_Affine.txt

# WaAaWtAt - CORRECT!!!!!!!!!!!!
$(subjDir)/T1_brain_in_MNI_WaAaWtAt.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_MNI_WaAaWtAt.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	ANTStemplate_to_MNI_Warp.nii.gz \
	ANTStemplate_to_MNI_Affine.txt \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt
	
# AtWtAaWa
$(subjDir)/T1_brain_in_MNI_AtWtAaWa.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_MNI_AtWtAaWa.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz \
	ANTStemplate_to_MNI_Affine.txt \
	ANTStemplate_to_MNI_Warp.nii.gz

# AaWaAtWt
$(subjDir)/T1_brain_in_MNI_AaWaAtWt.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_MNI_AaWaAtWt.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	ANTStemplate_to_MNI_Affine.txt \
	ANTStemplate_to_MNI_Warp.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz

# AaAtWaWt
$(subjDir)/T1_brain_in_MNI_AaAtWaWt.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_MNI_AaAtWaWt.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	ANTStemplate_to_MNI_Affine.txt \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt \
	ANTStemplate_to_MNI_Warp.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz
	
# WtWaAtAa
$(subjDir)/T1_brain_in_MNI_WtWaAtAa.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_MNI_WtWaAtAa.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz \
	ANTStemplate_to_MNI_Warp.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt \
	ANTStemplate_to_MNI_Affine.txt

# WaWtAaAt
$(subjDir)/T1_brain_in_MNI_WaWtAaAt.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_MNI_WaWtAaAt.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz \
	ANTStemplate_to_MNI_Warp.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt \
	ANTStemplate_to_MNI_Affine.txt

# AtAaWtWa
$(subjDir)/T1_brain_in_MNI_AtAaWtWa.nii.gz:
	${ANTSPATH}WarpImageMultiTransform 3 \
	$(subjDir)/T1_brain.nii.gz \
	$(subjDir)/T1_brain_in_MNI_AtAaWtWa.nii.gz \
	-R ${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz \
	$(subjDir)/T1_to_ANTStemplate_Affine.txt \
	ANTStemplate_to_MNI_Affine.txt \
	$(subjDir)/T1_to_ANTStemplate_Warp.nii.gz \
	ANTStemplate_to_MNI_Warp.nii.gz
