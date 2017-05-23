#!/bin/bash
subjDir=$(shell pwd)
subjects=$(wildcard $(subjDir)/3[0-9][0-9][0-9])
subj_raw=$(subjects:%=%~raw.hpf)
#subj_despike=$(subjects:%=%~despike.hpf)
#subj_tshift=$(subjects:%=%~tshift.hpf)
#subj_volreg=$(subjects:%=%~volreg.hpf)
#subj_regress=$(subjects:%=%~regress.hpf)
#subj_raw=$(subjects:%=%~raw)
#subj_despike=$(subjects:%=%~despike)
#subj_tshift=$(subjects:%=%~tshift)
#subj_volreg=$(subjects:%=%~volreg)
#subj_regress=$(subjects:%=%~regress)
.PHONY: all

#all: $(subj_raw) $(subj_despike) $(subj_tshift) $(subj_volreg) $(subj_regress)
all: $(subj_raw)

%: FORCE

$(subjDir)/%: FORCE
	session=$$(echo $* | cut -d'~' -f 1); \
	filePrefix=$$(echo $* | cut -d'~' -f 2); \
	proc_pic_iowa $$session

FORCE:
	
