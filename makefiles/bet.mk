rawDir=/projects2/udall/rest/raw
restDir=/projects2/udall/standardFC
sessions=$(notdir $(wildcard $(restDir)/RC4[1-2][0-9][0-9]-[1-2]))
rest_brain=$(sessions:%=%/rest_brain.nii.gz)
# Set open MP number of threads to be 1, so that we can parallelize
# using make.
export OMP_NUM_THREADS=1

all: $(rest_brain)

RC4%/rest_brain.nii.gz:
	session=RC4$* ;\
        echo $@ ;\
	mkdir -p $$session ;\
        cd $$session ;\
	bet $(rawDir)/$${session}_rest.nii.gz rest_brain.nii.gz -F -m



