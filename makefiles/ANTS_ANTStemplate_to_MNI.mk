
# Set the default number of threads used by ANTS to 1 to parallelize on 
# gridengine
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

all: ANTStemplate_to_MNI_deformed.nii.gz

ANTStemplate_to_MNI_deformed.nii.gz:
	${ANTSPATH}antsIntroduction.sh -d 3 -i udall_template.nii.gz -m 30x90x20 -o ANTStemplate_to_MNI_ -s CC -r /usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz -t GR
