SHELL=/bin/sh

udallDir=/projects2/udall
restDir=$(udallDir)/standardFC
subjectList=$(notdir $(wildcard $(restDir)/RC4[1-2][0-9][0-9]))

# Set open MP number of threads to be 1, so that we can parallelize
# using make.
export OMP_NUM_THREADS=1
subjectList=RC4103

rawDir=/projects2/udall/rest/raw
scriptDir=/projects2/udall/pboord/pic/scripts
.PHONY: all
all: $(subjectList)

RC4%:
	for sessionNr in 1 2; do \
		funcNii=$(rawDir)/$@-$${sessionNr}_rest.nii ;\
		funcNiiGz=$${funcNii}.gz ;\
		if [ ! -e $$funcNiiGz ]; then \
			if [ -e $$funcNii ]; then \
				gzip $$funcNiiGz ;\
			else \
				echo No data for $@-$${sessionNr} ;\
				exit ;\
			fi ;\
		fi ;\
	done ;\
	cmd="$(scriptDir)/preproc_pestica $@" ;\
	echo $$cmd ;\
	eval $$cmd ;\

