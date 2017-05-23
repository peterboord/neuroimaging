#!/bin/bash
cwd=$(shell pwd)
dataDir=/project_space/pboord/PIC/iowa/raw
compareDir=$(dataDir)/compare/noRetroicor
#seeds=PCC WM
seeds=PCC
#procSteps=$(wildcard $(compareDir)/*)
procSteps=$(compareDir)/raw.hpf_despike.hpf
.PHONY: all

all: $(procSteps)

%: FORCE

$(compareDir)/%: FORCE
	stepA=$$(echo $* | cut -d'_' -f 1); \
	stepB=$$(echo $* | cut -d'_' -f 2); \
	for seed in $(seeds); do \
		mkdir -p $@/pic/$$seed; \
		fslmerge -t $@/pic/$$seed/4dcope.nii $$(ls $(dataDir)/3???/3???.results/$$stepA/pic/$$seed/0.nii) $$(ls $(dataDir)/3???/3???.results/$$stepB/pic/$$seed/0.nii); \
	done

FORCE:
	
