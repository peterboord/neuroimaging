#!/bin/bash
subjDir=$(shell pwd)
subjects=$(wildcard $(subjDir)/3[0-9][0-9][0-9])
subj_raw=$(subjects:%=%~cov_raw.hpf)
subj_despike=$(subjects:%=%~cov_despike.hpf)
subj_tshift=$(subjects:%=%~cov_tshift.hpf)
subj_volreg=$(subjects:%=%~cov_volreg.hpf)
subj_regress=$(subjects:%=%~cov_regress.hpf)
.PHONY: all

all: $(subj_raw) $(subj_despike) $(subj_tshift) $(subj_volreg) $(subj_regress)

%: FORCE

$(subjDir)/%: FORCE
	session=$$(echo $* | cut -d'~' -f 1); \
	filePrefix=$$(echo $* | cut -d'~' -f 2); \
	post_proc_pic_iowa $$session $$filePrefix

FORCE:
	
