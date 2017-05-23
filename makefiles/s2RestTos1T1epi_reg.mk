
sessions=$(wildcard RC4[0-9][0-9][0-9]-2)
subjects=$(subst -2,,$(sessions))
ANTSdir=/projects2/udall/ANTS
restDir=/projects2/udall/standardFC
rest_to_s1T1=$(sessions:%=%/rest_to_s1T1.mat)
rawDir=/projects2/udall/rest/raw
# Set open MP number of threads to be 1, so that we can parallelize
# using make.
export OMP_NUM_THREADS=1

all: $(rest_to_s1T1)

RC4%/rest_to_s1T1.mat:
	session=RC4$* ;\
	subject=$$(echo $$session|cut -d'-' -f 1) ;\
	echo $@ ;\
	mkdir -p $$session ;\
	cd $$session ;\
	fslroi $(rawDir)/$${session}_rest.nii.gz $${session}_rest_0.nii.gz $$(($$(fslnvols $(rawDir)/$${session}_rest.nii.gz) / 2)) 1 ;\
	epi_reg --epi=$${session}_rest_0.nii.gz --t1=$(ANTSdir)/$$subject/T1.nii.gz --t1brain=$(ANTSdir)/$$subject/T1_brain.nii.gz --out=rest_to_s1T1.mat
