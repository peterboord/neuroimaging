
restDir=/projects2/udall/standardFC
sessions=$(notdir $(wildcard $(restDir)/RC4[1-2][0-9][0-9]-[1-2]))
prefiltered_func_data_mcf=$(sessions:%=%/prefiltered_func_data_mcf.par)
rawDir=/projects2/udall/rest/raw

all: $(prefiltered_func_data_mcf)

RC4%/prefiltered_func_data_mcf.par:
	echo $@ ;\
	session=RC4$* ;\
	mkdir -p $$session ;\
	cd $$session ;\
	mcflirt -nn_final -stats -mats -plots prefilted_func_data_mcf.par -in $(rawDir)/$${session}_rest.nii.gz -out prefiltered_func_data_mcf
